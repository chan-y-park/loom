import logging
import numpy
import sympy

from cmath import exp, pi
from math import floor
from scipy import interpolate

from geometry import find_xs_at_z_0, align_sheets_for_e_6_ffr
from misc import (
    cpow, remove_duplicate, ctor2, r2toc, delete_duplicates, is_root
)

x, z = sympy.symbols('x z')

# TODO: Generalize this.
# Number of x's at a fixed z
NUM_ODE_XS_OVER_Z = 2

# FIXME: Use the configuration data 'branch_point_cutoff'
# within a disc of such radius from any branch point,
# the intersection of a S-wall originating from there
# with the corresponding cut, will be ignored.
# BRANCH_POINT_RADIUS = 0.01 

# Desired precision on the phase of seeds
# Warning: setting it too small will bring the seeding point
# too close to a branch point.
SEED_PHASE_PRECISION = 0.001
SEED_PRECISION_MAX_DEPTH = 5

# FIXME: Use the configuration data 'puncture_cutoff'
# cut S-walls if they get too close to punctures
# PUNCTURE_RADIUS = 0.001


class Joint:
    def __init__(self, z=None, s_wall_1=None, s_wall_2=None,
                 t_1=None, t_2=None, sw_data=None, logger_name='loom',):
        self.z = None
        self.M = None
        self.parents = None
        self.label = None
        self.root = None
        self.ode_xs = None

        if z is None:
            # Return without setting attributes.
            return

        self.z = z
        self.M = s_wall_1.M[t_1] + s_wall_2.M[t_2]
        self.parents = [s_wall_1.label, s_wall_2.label]
        self.label = [s_wall_1.label, s_wall_2.label]

        alpha_1 = s_wall_1.get_root_at_t(t_1)
        alpha_2 = s_wall_2.get_root_at_t(t_2) 
        self.root = alpha_1 + alpha_2

        # FIXME: The following, including self.ode_xs, can be removed
        # once the seeding of an S-wall is done by using a root.
        ffr_xs_at_z = sw_data.get_sheets_at_z(z, ffr=True).values()
        ffr_new_wall_weight_pairs = (
            sw_data.g_data.ordered_weight_pairs(self.root, ffr=True)
        )
        ffr_w_p_0 = ffr_new_wall_weight_pairs[0]
        ode_x1 = ffr_xs_at_z[ffr_w_p_0[0]]
        ode_x2 = ffr_xs_at_z[ffr_w_p_0[1]]
        self.ode_xs = [ode_x1, ode_x2]

    def __eq__(self, other):
        return self.label == other.label

    def get_json_data(self):
        json_data = {
            'z': ctor2(self.z),
            'M': ctor2(self.M),
            'parents': [parent for parent in self.parents],
            'label': self.label,
            'root': self.root.tolist(),
            'ode_xs': [ctor2(x) for x in self.ode_xs],
        }
        return json_data

    def set_from_json_data(self, json_data):
        self.z = r2toc(json_data['z'])
        self.M = r2toc(json_data['M'])
        self.parents = [parent for parent in json_data['parents']]
        self.label = json_data['label']
        self.root = numpy.array(json_data['root'])
        self.ode_xs = [r2toc(x) for x in json_data['ode_xs']]

    def is_equal_to(self, other, accuracy):
        if(abs(self.z - other.z) > accuracy):
            return False

        if numpy.array_equal(self.root, other.root) is not True:
            return False

        if(abs(self.M - other.M) > accuracy):
            return False

        # FIXME: Probably the following comparison is not needed 
        # after comparing the roots of two S-walls.
        # Decide whether to do the comparison or not.
        for i in range(NUM_ODE_XS_OVER_Z):
            if abs(self.ode_xs[i] - other.ode_xs[i]) > accuracy:
                return False

        return True


# TODO: Instead of seeding with x_0, it would make 
# more sense to give the initial root-type.
# From that, and the trivialization module, we can
# extract the value of x_0, in principle.


class SWall(object):
    def __init__(self, z_0=None, x_0=None, M_0=None, parents=None,
                 label=None, n_steps=None, logger_name='loom'):
        """
        SWall.z is a NumPy array of length n_steps+1,
        where z[t] is the base coordinate.

        SWall.x is a Numpy array of the fiber coordinates at t, i.e.
        SWall.x[t] = [x[t][0], x[t][1], ...]. 
        (NOTE: for now, we just use two sheets, they serve only for 
        numerical ode evolution. To get the sheets use the method
        SWall.get_sheets_at_t() instead.)
        """
        self.logger_name = logger_name
        # FIXME: self.zs & self.xs instead of z & x?
        if n_steps is None:
            self.z = []
            self.x = []
            self.M = []
        else:
            self.z = numpy.empty(n_steps + 1, complex)
            self.x = numpy.empty(n_steps + 1, (complex, NUM_ODE_XS_OVER_Z))
            self.M = numpy.empty(n_steps + 1, complex)
            self.z[0] = z_0
            self.x[0] = x_0
            self.M[0] = M_0
        self.parents = parents
        self.label = label

        # cuts_intersections = [[b_pt_idx, i, '(c)cw'], ...]
        self.cuts_intersections = []
        self.root_basepoint = []
        self.local_roots = []
        # local_weight_pairs is a list of pair of intgers.
        self.local_weight_pairs = []        

    def __setitem__(self, t, data):
        """
        Set the data of the S-wall at t, where
            data = [z[t], x[t][0], x[t][1], ..., M[t]]
        """
        self.z[t] = data[0]
        self.x[t] = data[1:NUM_ODE_XS_OVER_Z + 1]
        self.M[t] = data[NUM_ODE_XS_OVER_Z + 1]

    def __getitem__(self, t):
        """
        Get the data of the S-wall at t, where
            data = [z[t], x[t][0], x[t][1], ..., M[t]]
        """
        return numpy.concatenate([[self.z[t]], self.x[t], [self.M[t]]])

    def resize(self, size):
        """
        Resize z & x arrays to discard garbage data.
        """
        self.z.resize((size))
        self.x.resize((size, NUM_ODE_XS_OVER_Z))
        self.M.resize((size))

    def get_json_data(self):
        json_data = {
            'z': numpy.array([self.z.real, self.z.imag]).T.tolist(),
            'M': numpy.array([self.M.real, self.M.imag]).T.tolist(),
            'x': numpy.rollaxis(
                numpy.array([self.x.real, self.x.imag]), 0, 3
            ).tolist(),
            'parents': [parent for parent in self.parents],
            'label': self.label,
            'cuts_intersections': [
                [br_loc.label, t, d]
                for br_loc, t, d in self.cuts_intersections
            ],
            # 'root_basepoint': ctor2(self.root_basepoint),
            'local_roots': [root.tolist() for root in self.local_roots],
            'local_weight_pairs': self.local_weight_pairs,
        }
        return json_data

    def set_from_json_data(self, json_data, branch_loci):
        self.z = numpy.array([r2toc(z_t) for z_t in json_data['z']])
        self.M = numpy.array([r2toc(M_t) for M_t in json_data['M']])
        self.x = numpy.array(
            [[r2toc(x_i) for x_i in x_t] for x_t in json_data['x']]
        )
        self.parents = [parent for parent in json_data['parents']]
        self.label = json_data['label']
        self.cuts_intersections = []
        for br_loc_label, t, d in json_data['cuts_intersections']:
            for br_loc in branch_loci:
                if br_loc_label == br_loc.label:
                    self.cuts_intersections.append([br_loc, t, d])
        self.local_roots = numpy.array(json_data['local_roots'])
        self.local_weight_pairs = json_data['local_weight_pairs']
        # self.root_basepoint = r2toc(json_data['root_basepoint'])

    def get_splits(self, endpoints=False):
        splits = [t for bp, t, d in self.cuts_intersections]
        if endpoints is True:
            return [0] + splits + [len(self.z) - 1]
        else:
            return splits

    def grow(
        self,
        ode,
        branch_point_zs,
        puncture_point_zs,
        config,
    ):
        bpzs = branch_point_zs
        ppzs = puncture_point_zs
        num_of_steps = config['num_of_steps']
        size_of_small_step = config['size_of_small_step']
        size_of_large_step = config['size_of_large_step']
        size_of_bp_neighborhood = config['size_of_bp_neighborhood']
        size_of_puncture_cutoff = config['size_of_puncture_cutoff']
        mass_limit = config['mass_limit']

        step = 0
        z_i = self.z[0]
        M_i = self.M[0]
        # Prepare a 1-dim array for ode
        y_i = self[0]
        ode.set_initial_value(y_i)

        while ode.successful() and step < num_of_steps:
            step += 1
            # Stop if z is inside a cutoff of a puncture.
            if len(ppzs) > 0:
                min_d = min([abs(z_i - ppz) for ppz in ppzs])
                if min_d < size_of_puncture_cutoff:
                    self.resize(step)
                    break

            # Stop if M exceeds mass limit.
            if mass_limit is not None:
                if M_i > mass_limit:
                    self.resize(step)
                    break

            # Adjust the step size if z is near a branch point.
            if (len(bpzs) > 0 and
                (min([abs(z_i - bpz) for bpz in bpzs])
                 < size_of_bp_neighborhood)):
                dt = size_of_small_step * min([1.0, abs(y_i[1] - y_i[2])])
            else:
                dt = size_of_large_step * min([1.0, abs(y_i[1] - y_i[2])])

            y_i = ode.integrate(ode.t + dt)
            z_i = y_i[0]
            M_i = y_i[NUM_ODE_XS_OVER_Z + 1]
            self[step] = y_i

    def determine_root_types(self, sw_data, cutoff_radius=0,):
        """
        Determine at which points the wall crosses a cut, 
        for instance [55, 107, 231] would mean that 
        it changes root-type 3 times. 
        Then, pick a suitable point along the swall, away 
        from branch points or singularities, and determine 
        the root there. Finally, extend the determination 
        of the root type to other segments by following
        the wall across the various splits induced by cuts,
        both forward and backwards, using the Weyl monodromy.
        """
        logger = logging.getLogger(self.logger_name)
        g_data = sw_data.g_data
        # branching will occur at branch points or irregular singularities
        branch_loci = sw_data.branch_points + sw_data.irregular_singularities
        # Adding a minimal radius of 1.0 is necessary 
        # in case there is only a single branch point at z=0,
        # otherwise max_radius would be 0.
        max_radius = 2 * max(
            [abs(c_l.z) for c_l in branch_loci] + [1.0]
        )
        # Parametrize the z-coordinate of the k-wall's coordinates
        # as a function of the index.
        num_of_zs = len(self.z)
        traj_t = numpy.arange(num_of_zs)
        traj_z_r = self.z.real
        
        # Scan over branching loci cuts, see if path ever crosses one 
        # based on x-coordinates only, and create a list of intersections, 
        # each element being [br_loc, t, 'cw'|'ccw'].
        cuts_intersections = []
        for br_loc in branch_loci:
            br_loc_x = br_loc.z.real
            br_loc_y = br_loc.z.imag

            if num_of_zs > 3:
                # If the length of the S-wall's coordinates
                # is greater than 3, use the B-spline of SciPy
                # to find the intersections between cuts and the S-wall.
                g = interpolate.splrep(traj_t, traj_z_r - br_loc_x, s=0)

                # Produce a list of integers corresponding to indices of  
                # the S-wall's coordinate list that seem to cross branch-cuts
                # based on the z-coordinate's real part.
                t_zeros = interpolate.sproot(g).tolist()

            elif 1 < num_of_zs <= 3:
                # There are two or three data points on the S-wall.
                # Use a linear interpolation for every pair of data points.
                t_zeros = []
                for i in range(num_of_zs - 1):
                    y_a, y_b = (traj_z_r - br_loc_x)[i:i + 2]
                    if y_a * y_b < 0:
                        t_zeros.append(y_a / (y_a - y_b))
            else:
                raise RuntimeError('{} has only one data point.'
                                   .format(self.label))

            intersection_ts = []
            for t_zero in t_zeros:
                t = int(floor(t_zero))
                # Remove a duplicate, if any.
                if (t in intersection_ts):
                    continue
                # Enforce imaginary-part of z-coordinate intersection 
                # criterion: branch cuts extend vertically.
                if self.z[t].imag > br_loc_y:
                    intersection_ts.append(t)

            for t in intersection_ts:
                # Drop intersections of a primary S-wall with the 
                # branch cut emanating from its parent branch-point
                # if such intersections happens within a short 
                # distance from the starting point.
                if (
                    br_loc.__class__.__name__ == 'BranchPoint' and
                    br_loc.label == self.parents[0] and 
                    (abs(br_loc.z - self.z[t]) < cutoff_radius)
                ):
                    continue

                # Check that the intersection actually happens
                # and is not an artifact of the interpolation used above
                # which could become an extrapolation
                x_p = self.z[t - 1].real
                x_n = self.z[t + 1].real
                if not ((x_p < br_loc_x < x_n) or (x_n < br_loc_x < x_p)):
                    logger.info('Drop a fake cut intersection.')
                    continue

                # Add 
                # [the branch-point, t,
                #  the direction (either 'cw' or 'ccw')]
                # to each intersection.
                cuts_intersections.append(
                    [br_loc, t, clock(left_right(self.z, t))]
                )

        # Now sort intersections according to where they happen in proper 
        # time; recall that the elements of cuts_intersections are 
        # organized  as      [..., [branch_point_idx, t, 'ccw'] ,...]
        # where 't' is the integer of proper time at the intersection.
        self.cuts_intersections = sorted(
            cuts_intersections, cmp=lambda k1, k2: cmp(k1[1], k2[1])
        )
        logger.debug(
            'S-wall {} intersects the following cuts at the points\n{}.'
            .format(self.label, self.cuts_intersections)
        )

        # Add the actual intersection point to the S-wall
        # then update the attribute SWall.cuts_intersections accordingly    
        self.enhance_at_cuts(sw_data)

        # Choose a suitable point along the wall
        # we pick the one whose z coordinate's real part is 
        # farthest from critical loci of the fibration        
        # Ideally, we would like this basepoint to be not too far away 
        # on the C-plane, because getting close to infinity
        # means colliding sheets usually.
        # But some walls are born beyond the max_radius in general
        # in that case, we just choose the t=0 coordinate
        # Also we keep t < 100, to avoid numerical errors induced by 
        # numerics of ode_int, which tend to spoil the values of 
        # s_wall.x[t] far along the trajectory.
        t_0 = sorted(
            ([
                [t_i, min(
                    z_r_distance_from_ramification_loci(z_i, sw_data)
                )]
                for t_i, z_i in enumerate(self.z) if (
                    (abs(z_i) < max_radius or t_i == 0) and t_i < 100
                ) 
            ]), cmp=lambda a, b: cmp(a[1], b[1])
        )[-1][0]

        z_0 = self.z[t_0]
        xs_0 = self.x[t_0]
        self.root_basepoint = [t_0, z_0, xs_0]

        # Determine the initial root-type
        initial_root = get_s_wall_root(z_0, xs_0, sw_data,)

        if is_root(initial_root, sw_data.g_data) is False:
            logging.info(
                'Warning: could not assign a root to {}'
                .format(self.label)
            )

        # A list of ordered pairs [...[i, j]...]
        # such that weights[j] - weights[i] = root
        initial_weight_pairs = (
            sw_data.g_data.ordered_weight_pairs(initial_root,)
        )
            
        self.local_roots = [initial_root]
        self.local_weight_pairs = [initial_weight_pairs]
        
        if len(self.cuts_intersections) > 0:
            # Note that we reverse the time-ordering!
            intersections_before_t_0 = [
                i_p for i_p in self.cuts_intersections if i_p[1] < t_0
            ][::-1]

            intersections_after_t_0 = [
                i_p for i_p in self.cuts_intersections if i_p[1] > t_0
            ]

            # Fill in the root types that occur after the basepoint
            for k in range(len(intersections_after_t_0)):
                br_loc, t, direction = intersections_after_t_0[k]

                current_root = self.local_roots[-1]
                new_root = g_data.weyl_monodromy(
                    current_root, br_loc, direction
                )
                new_weight_pairs = g_data.ordered_weight_pairs(new_root)

                self.local_roots.append(new_root)
                self.local_weight_pairs.append(new_weight_pairs)

            # Fill in the root types that occur before the basepoint
            # recall that their time-ordering has already been reversed
            # so the first one in the list is the closest to t_0, and so on
            for k in range(len(intersections_before_t_0)):
                br_loc, t, direction = intersections_before_t_0[k]

                current_root = self.local_roots[0]
                new_root = g_data.weyl_monodromy(
                    current_root, br_loc, direction, reverse=True
                )
                new_weight_pairs = g_data.ordered_weight_pairs(new_root)

                self.local_roots.insert(0, new_root)
                self.local_weight_pairs.insert(0, new_weight_pairs)

    def get_root_at_t(self, t):
        """
        Given an integer t which parametrizes a point 
        of the trajectory in proper time, return the local 
        root at that point.
        """
        if t < 0 or t > (len(self.z) - 1):
            raise ValueError
        else:
            closed_splits = self.get_splits() + [len(self.z) - 1]
            for i, sp in enumerate(closed_splits):
                if t <= sp:
                    return self.local_roots[i]
                    break
                else:
                    pass

    def get_weight_pairs_at_t(self, t):
        """
        Given an integer t which parametrizes a point 
        of the trajectory in proper time, return the local 
        pair of weights at that point.
        """        
        if t < 0 or t > (len(self.z) - 1):
            raise ValueError
        else:
            closed_splits = self.get_splits() + [len(self.z) - 1]
            for i, sp in enumerate(closed_splits):
                if t <= sp:
                    return self.local_weight_pairs[i]
                    break
                else:
                    pass

    def get_sheets_at_t(self, t, sw_data):
        """
        Given an integer t which parametrizes a point 
        of the trajectory in proper time, return 
        a list of pairs of values.
        For sheet labels [... [i, j] ...] at given t,
        return [... [x_i, x_j] ...]
        """
        z = self.z[t]
        # the following is a dictionary
        xs_at_z = sw_data.get_sheets_at_z(z)
        weight_pairs = self.get_weight_pairs_at_t(t)
        return [[xs_at_z[w_p[0]], xs_at_z[w_p[1]]] for w_p in weight_pairs]

    def enhance_at_cuts(self, sw_data):
        """
        Add the intersection points of Swalls and branch cuts
        also update the intersection data accordingly
        """
        if len(self.cuts_intersections) == 0:
            return None

        wall_pieces_z = []
        wall_pieces_x = []
        wall_pieces_M = []

        # split the wall into pieces, at the end of each 
        # piece add the corresponding intersection point
        t_0 = 0
        for int_data in self.cuts_intersections:
            br_loc, t, chi = int_data

            z_1 = self.z[t]
            z_2 = self.z[t + 1]
            
            z_to_add = get_intermediate_z_point(z_1, z_2, br_loc.z)

            xs_1 = self.x[t]
            xs_2 = self.x[t + 1]
            xs_to_add = [
                get_intermediate_value(xs_1[i], xs_2[i], z_1, z_2, z_to_add)
                for i in range(NUM_ODE_XS_OVER_Z)
            ]

            M_to_add = get_intermediate_value(
                self.M[t], self.M[t + 1], z_1, z_2, z_to_add
            )

            z_piece = numpy.concatenate(
                (self.z[t_0:t + 1], numpy.array([z_to_add], dtype=complex))
            )
            x_piece = numpy.concatenate(
                (self.x[t_0:t + 1], numpy.array([xs_to_add], dtype=complex))
            )
            M_piece = numpy.concatenate(
                (self.M[t_0:t + 1], numpy.array([M_to_add], dtype=complex))
            )
            wall_pieces_z.append(z_piece)
            wall_pieces_x.append(x_piece)
            wall_pieces_M.append(M_piece)
            t_0 = t + 1

        # Get the last piece of the wall.
        wall_pieces_z.append(self.z[t_0:])
        wall_pieces_x.append(self.x[t_0:])
        wall_pieces_M.append(self.M[t_0:])

        # Make updates. 
        self.z = numpy.concatenate(wall_pieces_z)
        self.x = numpy.concatenate(wall_pieces_x)
        self.M = numpy.concatenate(wall_pieces_M)

        # Update the intersection data.
        new_cuts_intersections = []
        for i, int_data in enumerate(self.cuts_intersections):
            br_loc, t_old, chi = int_data
            t_new = t_old + i + 1
            new_cuts_intersections.append([br_loc, t_new, chi])
        
        self.cuts_intersections = new_cuts_intersections


def get_s_wall_root(z, ffr_xs, sw_data):
    x_i, x_j = ffr_xs
    
    # Recall that S-wall numerical evolution
    # is based on the first fundamental representation.
    # In particular, the xs above are values of sheets 
    # from the 1st fundamental rep cover, and should 
    # be compared with the corresponding trivialization.
    # This is taken care of by setting the key argument
    # ffr=True when calling get_sheets_at_z.
    # The following is a dictionary
    sheets_at_z = sw_data.get_sheets_at_z(z, ffr=True)
    xs_at_z = sheets_at_z.values()
    
    # Sheet matching x_i
    closest_to_x_i = sorted(xs_at_z, key=lambda x: abs(x - x_i))[0]
    i = [k for k, v in sheets_at_z.iteritems() if v == closest_to_x_i][0]

    # Sheet matching x_j
    closest_to_x_j = sorted(xs_at_z, key=lambda x: abs(x - x_j))[0]
    j = [k for k, v in sheets_at_z.iteritems() if v == closest_to_x_j][0]

    return sw_data.g_data.ffr_weights[j] - sw_data.g_data.ffr_weights[i]


def get_intermediate_z_point(z_1, z_2, bp_z_med):
    """
    get the intermediate point between z_1 and z_2
    in correspondence of the real part of z_med
    """
    # FIXME: division by zero may happend when x_1 and x_2 are too close.
    x_1 = z_1.real
    y_1 = z_1.imag
    x_2 = z_2.real
    y_2 = z_2.imag
    x_med = bp_z_med.real
    slope = (y_2 - y_1) / (x_2 - x_1)
    y_med = y_1 + slope * (x_med - x_1)
    return x_med + 1j * y_med


def get_intermediate_value(v_1, v_2, z_1, z_2, z_med):
    """
    Along an S-wall, get the intermediate value between v_1 and v_2
    in correspondence of z_med
    """
    # FIXME: division by zero may happend when z_1 and z_2 are too close.
    slope = (v_2 - v_1) / (z_2 - z_1)
    v_med = v_1 + slope * (z_med - z_1)
    return v_med


def get_s_wall_seeds(sw, theta, branch_point, config, logger_name):
    """
    S-walls are seeded from branch points.
    Each branch point has a number of ramification 
    points lying above it.
    Regardless of the representation, it is sufficient
    to consider one of these ramification points
    to extract the seed data.
    We thus stick to (any)one ramification point of the 
    fundamental representation to get the seeds.
    It is assumed that we are working in the first 
    fundamental representation of the Lie algebra.
    """
    logger = logging.getLogger(logger_name)

    # FIXME: reintroduce the handling of massless punctures
    # see previous versions of this function, left above in comment.

    accuracy = config['accuracy']
    initial_seed_size = config['size_of_small_step']
    seeds = []
    # max_r_i = max([rp.i for rp in branch_point.ffr_ramification_points])
    min_dt = 1.0

    for rp in branch_point.ffr_ramification_points:
        z_0 = rp.z
        x_0 = rp.x
        r_i = rp.i
        rp_type = rp.ramification_type
        # sw_diff_coeff = rp.sw_diff_coeff
        sw_diff_coeffs_a_b = rp.sw_diff_coeffs_a_b
        logger.debug('Analyze ramification point (z,x)={}'.format([z_0, x_0]))
        logger.debug('Ramification index = {}'.format(r_i))
        logger.debug('Ramification type = {}'.format(rp_type))
        logger.debug(
            'leading coefficients of SW diff: a = {}\t b={}\n'
            .format(sw_diff_coeffs_a_b[0], sw_diff_coeffs_a_b[1])
        )

        # Construct the seeding points for the branch point
        # by studying the type of ramification structure of the r.p.
        if rp_type == 'type_I':
            phases = [exp(2 * pi * 1j * float(i) / r_i) for i in range(r_i)]
            phi = [[p1 - p2 for p1 in phases] for p2 in phases]
            
            omega = exp(2.0 * pi * 1j * float(r_i) / float(r_i + 1))

            a, b = sw_diff_coeffs_a_b
            dz_phases = ([
                (1.0 / cpow((-1.0 * a / b), 1, r_i + 1)) *
                exp(1j * theta * float(r_i) / (r_i + 1)) *
                ((-1.0 / phi[i][j]) ** (float(r_i) / (r_i + 1))) * (omega ** s)
                for i in range(r_i) for j in range(r_i) 
                for s in range(r_i + 1) if i != j
            ])

            norm_dz_phases = [d / abs(d) for d in dz_phases]
            # these are the normalized phases of the seeds
            # with respect to the branch point:
            zetas = remove_duplicate(
                norm_dz_phases, lambda p1, p2: abs(p1 - p2) < (accuracy)
            )

            print 'a = {}\nb = {}'.format(a, b)
            print 'phases = {}'.format(zetas)
        
        elif rp_type == 'type_II':
            if r_i % 2 == 1:
                raise Exception('Cannot have a type II ramification point' +
                                'with odd ramification index.')
            # defining this object just for enhanced readability of code 
            # in comparing with notes on classification of ramifications
            r_k = r_i / 2
            phases = [
                exp(2 * pi * 1j * float(i) / (2.0 * r_k)) for i in range(r_k)
            ]
            phi = [[p1 - p2 for p1 in phases] for p2 in phases]
            psi = [[
                (phases[i] + phases[j]) * numpy.sign(i - j) for i in range(r_k)
            ] for j in range(r_k)]
            
            a, b = sw_diff_coeffs_a_b
            omega = exp(2.0 * pi * 1j * float(2 * r_k) / float(2 * r_k + 1))

            dz_phases = ([
                (1.0 / cpow((-1.0 * a / b), 1, 2 * r_k + 1)) *
                exp(1j * theta * float(2 * r_k) / (2 * r_k + 1)) *
                ((-1.0 / phi[i][j]) ** (float(2 * r_k) / (2 * r_k + 1))) 
                * (omega ** s)
                for i in range(r_k) for j in range(r_k) 
                for s in range(2 * r_k + 1) if i != j
            ] + [
                (1.0 / cpow((-1.0 * a / b), 1, 2 * r_k + 1)) *
                exp(1j * theta * float(2 * r_k) / (2 * r_k + 1)) *
                ((-1.0 / psi[i][j]) ** (float(2 * r_k) / (2 * r_k + 1))) 
                * (omega ** s)
                for i in range(r_k) for j in range(r_k) 
                for s in range(2 * r_k + 1) if i != j
            ])

            norm_dz_phases = [d / abs(d) for d in dz_phases]
            # these are the normalized phases of the seeds
            # with respect to the branch point:
            zetas = remove_duplicate(
                norm_dz_phases, lambda p1, p2: abs(p1 - p2) < accuracy
            )

        elif rp_type == 'type_III':
            if r_i % 2 == 1:
                raise Exception('Cannot have a type III ramification point' +
                                'with odd ramification index.')
            # defining this object just for enhanced readability of code 
            # in comparing with notes on classification of ramifications
            r_k = r_i / 2

            phases = [
                exp(2 * pi * 1j * float(i) / (2.0 * (r_k - 1))) 
                for i in range(r_k - 1)
            ] + [0.0]
            phi = [[p1 - p2 for p1 in phases] for p2 in phases]
            psi = [[
                (phases[i] + phases[j]) * numpy.sign(i - j) 
                for i in range(r_k)
            ] for j in range(r_k)]
            # print 'phi = {}'.format(phi)
            omega = exp(
                2.0 * pi * 1j * float(2 * r_k - 2) / float(2 * r_k - 1)
            )

            a, b = sw_diff_coeffs_a_b
            dz_phases = ([
                (1.0 / cpow((-1.0 * a / b), 1, 2 * r_k - 1)) *
                exp(1j * theta * float(2 * r_k - 2) / (2 * r_k - 1)) *
                ((-1.0 / phi[i][j]) ** (float(2 * r_k - 2) / (2 * r_k - 1))) * 
                (omega ** s)
                for i in range(r_k) for j in range(r_k) 
                for s in range(2 * r_k - 1) if i != j
            ] + [
                (1.0 / cpow((-1.0 * a / b), 1, 2 * r_k - 1)) *
                exp(1j * theta * float(2 * r_k - 2) / (2 * r_k - 1)) *
                ((-1.0 / psi[i][j]) ** (float(2 * r_k - 2) / (2 * r_k - 1))) * 
                (omega ** s)
                for i in range(r_k) for j in range(r_k) 
                for s in range(2 * r_k - 1) if i != j
            ])

            norm_dz_phases = [d / abs(d) for d in dz_phases]
            # these are the normalized phases of the seeds
            # with respect to the branch point:
            zetas = remove_duplicate(
                norm_dz_phases, lambda p1, p2: abs(p1 - p2) < accuracy
            )

        elif rp_type == 'type_IV':
            phases = [
                exp(2 * pi * 1j * float(i) / (12.0)) 
                for i in range(12)
            ] + [0.0]
            phi = [[
                phases[i] - phases[j] for i in range(13) 
            ] for j in range(13)]

            # print 'phi = {}'.format(phi)
            omega = exp(
                2.0 * pi * 1j / 13.0
            )
            a, b = sw_diff_coeffs_a_b
            dz_phases = ([
                (1.0 / cpow((-1.0 * a / b), 1, 13)) *
                exp(1j * theta * 12.0 / 13.0) *
                ((-1.0 / phi[i][j]) ** (12.0 / 13.0)) * 
                (omega ** s)
                for i in range(13) for j in range(13) 
                for s in range(13) if (i != j and e_6_compatible(i, j))
            ])
            norm_dz_phases = [d / abs(d) for d in dz_phases]
            # these are the normalized phases of the seeds
            # with respect to the branch point:
            zetas = remove_duplicate(
                norm_dz_phases, lambda p1, p2: abs(p1 - p2) < accuracy
            )

        # Now for each seeding point z_1 we identify two sheets
        # of the cover which match the phase of the displacement z_1-z_0
        for zeta in zetas:
            raise_precision = True
            precision_level = 0
            while (
                raise_precision is True and 
                precision_level <= SEED_PRECISION_MAX_DEPTH
            ):
                logger.debug(
                    'Seeding precision level = {}'.format(precision_level)
                )
                dt = initial_seed_size * 0.1 ** (precision_level)
                if dt < min_dt:
                    min_dt = dt
                # z_1 = z_0 + accuracy * zeta
                z_1 = z_0 + dt * zeta

                if rp_type == 'type_I':
                    x_s = find_xs_at_z_0(sw, z_1, x_0, r_i, ffr=True)
                    # print '\n\nat z_1={} the sheets are {}'.format(z_1, x_s)
                    # a list of the type
                    # [... [phase, [x_i, x_j]] ...]
                    x_i_x_j_phases = []
                    for i, x_i in enumerate(x_s): 
                        for j, x_j in enumerate(x_s):
                            if i != j:
                                # v_i = complex(
                                #     sw.diff.num_v.subs([(z, z_1), (x, x_i)])
                                # )
                                # v_j = complex(
                                #     sw.diff.num_v.subs([(z, z_1), (x, x_j)])
                                # )
                                # if (v_j - v_i) != 0:
                                #     ij_factor = (
                                #         -1.0 * exp(1j * theta) / (x_j - x_i)
                                #     )
                                #     x_i_x_j_phases.append(
                                #         [
                                #             (ij_factor) / abs(ij_factor), 
                                #             [x_i, x_j]
                                #         ]
                                #     )
                                ij_factor = (
                                    -1.0 * exp(1j * theta) / (x_j - x_i)
                                )
                                x_i_x_j_phases.append(
                                    [
                                        (ij_factor) / abs(ij_factor), 
                                        [x_i, x_j]
                                    ]
                                )

                elif rp_type == 'type_II' or rp_type == 'type_III':
                    # we assume that the ramification index is maximal
                    # therefore we ask for all the sheets at z_1.
                    x_s = find_xs_at_z_0(sw, z_1, ffr=True)

                    # order of magnitude of expected separation 
                    # of sheets at z_1
                    dx = abs(a / b) * (dt ** (1.0 / float(r_i)))
                    x_accuracy = min([accuracy, dx])
                    # OLD
                    # dx = abs(sw_diff_coeff) * (dt ** (1.0 / float(r_i)))
                    # x_accuracy = min([accuracy, dx])

                    # a list of the type
                    # [... [phase, [x_i, x_j]] ...]
                    # where we eclude x_i=x_j and x_i=-x_j 
                    # since in D-type there are no roots 
                    # between a weight v and -v.
                    x_i_x_j_phases = []
                    for i, x_i in enumerate(x_s):
                        for j, x_j in enumerate(x_s):
                            if (
                                abs(x_j - x_i) > x_accuracy 
                                and abs(x_j + x_i) > x_accuracy
                            ):
                                # v_i = complex(
                                #     sw.diff.num_v.subs([(z, z_1), (x, x_i)])
                                # )
                                # v_j = complex(
                                #     sw.diff.num_v.subs([(z, z_1), (x, x_j)])
                                # )
                                # ij_factor = (
                                #     -1.0 * exp(1j * theta) / (v_j - v_i)
                                # )
                                # # ij_factor = -1.0 * exp(1j*theta)/(x_j - x_i)
                                # x_i_x_j_phases.append(
                                #     [(ij_factor) / abs(ij_factor), [x_i, x_j]]
                                # )

                                ij_factor = (
                                    -1.0 * exp(1j * theta) / (x_j - x_i)
                                )
                                # ij_factor = -1.0 * exp(1j*theta)/(x_j - x_i)
                                x_i_x_j_phases.append(
                                    [(ij_factor) / abs(ij_factor), [x_i, x_j]]
                                )

                elif rp_type == 'type_IV':
                    x_s = find_xs_at_z_0(
                        sw, z_1, x_0, r_i, 
                        ffr=True, use_sage=True
                    )
                    # Note: although we are working near a degenerate locus,
                    # it is correct to enforce the option
                    # near_degenerate_branch_locus=False
                    aligned_xs = align_sheets_for_e_6_ffr(
                        x_s, 
                        sw.g_data.weights, 
                        near_degenerate_branch_locus=False
                    )
                    
                    # a list of the type
                    # [... [phase, [x_i, x_j]] ...]
                    x_i_x_j_phases = []
                    for i, x_i in enumerate(aligned_xs): 
                        w_i = sw.g_data.weights[i]
                        for j, x_j in enumerate(aligned_xs):
                            w_j = sw.g_data.weights[j]
                            if is_root(w_j - w_i, sw.g_data):
                                # v_i = complex(
                                #     sw.diff.num_v.subs([(z, z_1), (x, x_i)])
                                # )
                                # v_j = complex(
                                #     sw.diff.num_v.subs([(z, z_1), (x, x_j)])
                                # )
                                # # print '\n\ni ={}, j={}'.format(i,j)
                                # # print 'w_j - w_i = {}'.format(w_j -w_i)
                                # # print 'x_j - x_i = {}'.format(x_j -x_i)
                                # # print 'v_j - v_i = {}'.format(v_j -v_i)
                                # if (v_j - v_i) != 0:
                                #     ij_factor = (
                                #         -1.0 * exp(1j * theta) / (v_j - v_i)
                                #     )
                                #     x_i_x_j_phases.append(
                                #         [
                                #             (ij_factor) / abs(ij_factor), 
                                #             [x_i, x_j]
                                #         ]
                                #     )

                                ij_factor = (
                                    -1.0 * exp(1j * theta) / (x_j - x_i)
                                )
                                x_i_x_j_phases.append(
                                    [
                                        (ij_factor) / abs(ij_factor), 
                                        [x_i, x_j]
                                    ]
                                )

                closest_pair = sorted(
                    x_i_x_j_phases, key=lambda p: abs(p[0] - zeta)
                )[0][1]
                phase_mismatch = abs(sorted(
                    x_i_x_j_phases, key=lambda p: abs(p[0] - zeta)
                )[0][0] - zeta)
                if phase_mismatch < SEED_PHASE_PRECISION:
                    logger.debug(
                        'Reached desired precision for seed'
                        '\nThe Mismatch between its phase '
                        'and that of the displacement '
                        'is : {}'.format(phase_mismatch)
                    )
                    raise_precision = False
                else:
                    if precision_level < SEED_PRECISION_MAX_DEPTH:
                        precision_level += 1
                    else:
                        logger.warning(
                            'Could not get the desired precision '
                            'on the seed of an S-wall.\nThe Mismatch between '
                            'the phase of a seed and that of the displacement '
                            'is : {}'.format(phase_mismatch)
                        )
                        raise_precision = False
                        break
            
            M_0 = 0
            seeds.append([z_1, closest_pair, M_0])

    # for higher-index ramification points we need greater accuracy to 
    # keep all the correct seeds, since dt is also their displacement
    # |z_1-z_0| we cannot just use dt, but must choose a small 
    # fraction of it
    seeds = delete_duplicates(seeds, lambda s: s[0], accuracy=(min_dt / 100))
    logger.debug('Number of S-walls emanating = {}'.format(len(seeds)))
    logger.debug('these are the seeds {}\n'.format(seeds))
    branch_point.seeds = seeds
    return seeds


def z_r_distance_from_ramification_loci(z, sw_data):
    critical_loci = sw_data.branch_points + sw_data.irregular_singularities
    return [abs(z.real - c_l.z.real) for c_l in critical_loci]


def branch_locus_from_label(sw_data, br_loc_label):
    branch_loci = sw_data.branch_points + sw_data.irregular_singularities
    for br_loc in branch_loci:
        if br_loc.label == br_loc_label:
            return br_loc
    raise Exception(
        'Could not find any branching locus labeled {}'.format(br_loc_label)
    )


def clock(direction):
    if direction == 'left':
        return 'ccw'
    elif direction == 'right':
        return 'cw'
    else:
        raise ValueError('Cannot read direction!')


def left_right(l, point):
    """
    given the list 
    l = [..., z, ...]
    and a point in the list (specified by the corresponding integer),
    determines whether x increases or decreases at that point, 
    returning repsectively 'left' or 'right'
    """
    if point > len(l) - 1:
        raise ValueError(
            'Can\'t determine direction, point doesnt belong to list!'
        )
    elif point > 0:
        if l[point - 1].real < l[point].real:
            return 'right'
        else:
            return 'left'
    elif point == 0:
        if l[point].real < l[point + 1].real:
            return 'right'
        else:
            return 'left'


def e_6_compatible(i, j):
    """
    Given two integers i,j =0, ..., 12
    determine whether the corresponding phases
    in the E_6 seeds can be paired together.
    This accounts for which sheets actually differ 
    by a root, or not. The following prescription
    is derived from the Coxeter projection diagram 
    of E_6, and looking for roots that connect
    a given weight with other weights.
    """
    dist = abs(i - j) % 12
    if (
        (0 <= i <= 12) and (0 <= i <= 12) and (dist <= 3 or dist >= 9) or 
        (i == 12 or j == 12)
    ):
        return True
    else:
        return False

