
import logging
import numpy
import sympy
# import pdb

from cmath import exp, pi
from math import floor
from scipy import interpolate

from geometry import find_xs_at_z_0
from misc import (cpow, remove_duplicate, ctor2, r2toc, delete_duplicates,
                  clock, left_right,)

x, z = sympy.symbols('x z')

# TODO: Generalize this.
# Number of x's at a fixed z
NUM_ODE_XS_OVER_Z = 2

# within a disc of such radius from any branch point,
# the intersection of a S-wall oroginating from there
# with the corresponding cut, will be ignored.
BRANCH_POINT_RADIUS = 0.01 

# Desired precision on the phase of seeds
# Warning: setting it too small will bring the seeding point
# too close to a branch point.
SEED_PHASE_PRECISION = 0.001
SEED_PRECISION_MAX_DEPTH = 5


class Joint:
    def __init__(self, z=None, s_wall_1=None, s_wall_2=None,
                 t_1=None, t_2=None, label=None, sw_data=None):
        self.z = None
        self.M = None
        self.parents = None
        self.label = None
        self.root = None
        self.ode_xs = None
        self.combos = None

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
        # self.indices = [t_1, t_2]

        # FIXME: Erase the following if we decide not to use
        # x-coords in calculating joints of S-walls.
        # Probably this will not be additional information because
        # x's are calculated from the x's of the ffr cover.
        # xs_at_z = find_xs_at_z_0(sw_data, z)
        new_wall_weight_pairs = sw_data.g_data.ordered_weight_pairs(self.root)
        # ### w_p_0 = new_wall_weight_pairs[0]
        # x_i_s = [xs_at_z[w_p[0]] for w_p in new_wall_weight_pairs]
        # x_j_s = [xs_at_z[w_p[1]] for w_p in new_wall_weight_pairs]
        # self.x_s = [[x_i_s[k], x_j_s[k]] for k in range(len(x_i_s))]

        # FIXME: The following, including self.ode_xs, can be removed
        # once the seeding of an S-wall is done by using a root.
        # ffr_xs_at_z = find_xs_at_z_0(sw_data, z, ffr=True)
        ffr_xs_at_z = sw_data.get_sheets_at_z(z, ffr=True).values()
        ffr_new_wall_weight_pairs = (
            sw_data.g_data.ordered_weight_pairs(self.root, ffr=True)
        )
        ffr_w_p_0 = ffr_new_wall_weight_pairs[0]
        ode_x1 = ffr_xs_at_z[ffr_w_p_0[0]]
        ode_x2 = ffr_xs_at_z[ffr_w_p_0[1]]
        self.ode_xs = [ode_x1, ode_x2]

        # For each pair of pairs: [i, j] from the first wall
        # and [j, k] from the second wall, we build 
        # the list of 'combos'.
        # This is a list [..., [[i, j], [j, k]],...]
        combos = []
        weight_pairs_1 = s_wall_1.get_weight_pairs_at_t(t_1)
        weight_pairs_2 = s_wall_2.get_weight_pairs_at_t(t_2)
        for w_p in new_wall_weight_pairs:
            w_i = w_p[0]
            w_k = w_p[1]
            ij_in_1 = [pair for pair in weight_pairs_1 if (pair[0] == w_i)]
            ij_in_2 = [pair for pair in weight_pairs_2 if (pair[0] == w_i)]
            jk_in_1 = [pair for pair in weight_pairs_1 if (pair[1] == w_k)]
            jk_in_2 = [pair for pair in weight_pairs_2 if (pair[1] == w_k)]
            if len(ij_in_1) > 0 and len(jk_in_2) > 0:
                combos.append([ij_in_1[0], jk_in_2[0]])
            elif len(ij_in_2) > 0 and len(jk_in_1) > 0:
                combos.append([ij_in_2[0], jk_in_1[0]])
            else:
                raise ValueError(
                    'Cannot pick ij, jk pairs from colliding S-walls\n{}\n{}'
                    .format(weight_pairs_1, weight_pairs_2)
                )
        self.combos = combos

    def __eq__(self, other):
        return self.label == other.label

    def get_json_data(self):
        # FIXME: Determine what data to save.
        json_data = {
            'z': ctor2(self.z),
            'M': ctor2(self.M),
            'parents': [parent for parent in self.parents],
            'label': self.label,
            'root': self.root.tolist(),
            'ode_xs': [ctor2(x) for x in self.ode_xs],
            'combos': self.combos,
        }
        return json_data

    def set_json_data(self, json_data):
        # FIXME: Determine what data to load.
        self.z = r2toc(json_data['z'])
        self.M = r2toc(json_data['M'])
        self.parents = [parent for parent in json_data['parents']]
        self.label = json_data['label']
        self.root = numpy.array(json_data['root'])
        self.ode_xs = [r2toc(x) for x in json_data['ode_xs']]
        self.combos = json_data['combos']

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
        # FIXME: Erase the following if we decide not to carry
        # ``x_s`` as an attribute of ``Joint``.
        # if(len(self.x_s) != len(other.x_s)):
        #     return False
        # for i in range(len(self.x_s)):
        #     if ((abs(self.x_s[i][0] - other.x_s[i][0]) > accuracy)
        #         or (abs(self.x_s[i][1] - other.x_s[i][1]) > accuracy)):
        #         return False
        return True


# TODO: Instead of seeding with x_0, it would make 
# more sense to give the initial root-type.
# From that, and the trivialization module, we can
# extract the value of x_0, in principle.


class SWall(object):
    def __init__(self, z_0=None, x_0=None, M_0=None, parents=None,
                 label=None, n_steps=None,):
        """
        SWall.z is a NumPy array of length n_steps+1,
        where z[t] is the base coordinate.

        SWall.x is a Numpy array of the fiber coordinates at t, i.e.
        SWall.x[t] = [x[t][0], x[t][1], ...]. 
        (NOTE: for now, we just use two sheets, they serve only for 
        numerical ode evolution. To get the sheets use the method
        SWall.get_sheets_at_t() instead.)
        """
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
            'cuts_intersections': self.cuts_intersections,
            'local_roots': [root.tolist() for root in self.local_roots],
            'local_weight_pairs': self.local_weight_pairs,
        }
        return json_data

    def set_json_data(self, json_data):
        self.z = numpy.array([r2toc(z_t) for z_t in json_data['z']])
        self.M = numpy.array([r2toc(M_t) for M_t in json_data['M']])
        self.x = numpy.array(
            [[r2toc(x_i) for x_i in x_t] for x_t in json_data['x']]
        )
        self.parents = [parent for parent in json_data['parents']]
        self.label = json_data['label']
        self.cuts_intersections = json_data['cuts_intersections']
        self.local_roots = numpy.array(json_data['local_roots'])
        self.local_weight_pairs = json_data['local_weight_pairs']

    def get_splittings(self):
        return [t for bp, t, d in self.cuts_intersections]

    def get_turning_points(self):
        """
        Return a list of indices of turning points of SWall.z,
        i.e. dx/dy = 0 or dy/dx = 0, where x = z[t].real and y = z[t].imag.
        """
        tps = []

        if len(self.z) < 3:
            return tps

        x_0 = self.z[0].real
        y_0 = self.z[0].imag

        for t in range(1, len(self.z) - 1):
            x_1 = self.z[t].real
            y_1 = self.z[t].imag
            x_2 = self.z[t + 1].real
            y_2 = self.z[t + 1].imag
            if (
                (x_1 - x_0) * (x_2 - x_1) < 0 or (y_1 - y_0) * (y_2 - y_1) < 0
            ): 
                tps.append(t)
            x_0 = x_1
            y_0 = y_1
        return tps

    def grow(
        self,
        ode,
        ramification_point_zs,
        puncture_point_zs,
        config,
        clipping_radius=None,
    ):
        rpzs = ramification_point_zs
        ppzs = puncture_point_zs
        # z_range_limits = config['z_range_limits']
        num_of_steps = config['num_of_steps']
        size_of_small_step = config['size_of_small_step']
        size_of_large_step = config['size_of_large_step']
        size_of_neighborhood = config['size_of_neighborhood']
        size_of_puncture_cutoff = config['size_of_puncture_cutoff']
        mass_limit = config['mass_limit']

        step = 0
        z_i = self.z[0]
        M_i = self.M[0]
        # Prepare a 1-dim array for ode
        y_i = self[0]
        ode.set_initial_value(y_i)

        # if z_range_limits is not None:
        #     (
        #         [z_real_min, z_real_max], [z_imag_min, z_imag_max] 
        #         = z_range_limits
        #     )

        while ode.successful() and step < num_of_steps:
            step += 1
            # Stop if z is inside a cutoff of a puncture.
            if len(ppzs) > 0:
                min_d = min([abs(z_i - ppz) for ppz in ppzs])
                if min_d < size_of_puncture_cutoff:
                    self.resize(step)
                    break

            # Stop if z is ouside the range limit.
            # if z_range_limits is not None:
            #     if (z_i.real < z_real_min or
            #         z_i.real > z_real_max or
            #         z_i.imag < z_imag_min or
            #         z_i.imag > z_imag_max):
            #         self.resize(step)
            #         break

            # Stop if M exceeds mass limit.
            if mass_limit is not None:
                if M_i > mass_limit:
                    self.resize(step)
                    break

            # Adjust the step size if z is near a branch point.
            if (
                len(rpzs) > 0 and
                min([abs(z_i - rpz) for rpz in rpzs]) < size_of_neighborhood
            ):
                dt = size_of_small_step
            else:
                dt = size_of_large_step

            y_i = ode.integrate(ode.t + dt)
            z_i = y_i[0]
            M_i = y_i[NUM_ODE_XS_OVER_Z + 1]
            self[step] = y_i
        if clipping_radius is not None:
            i_clip = None
            for i, z_i in enumerate(self.z):
                if abs(z_i) > clipping_radius:
                    i_clip = i
                    break
            if i_clip is not None:
                self.z = self.z[:i_clip]
                self.x = self.x[:i_clip]


    def determine_root_types(self, sw_data):
        """
        Determine at which points the wall crosses a cut, 
        for instance [55, 107, 231] would mean that 
        it changes root-type 3 times. 
        Then, pick a suitable point along the swall, away 
        from branch points or singularities, and determine 
        the root there. Finally, extend the determination 
        of the root type to other segments by following
        the wall across the various splittings induced by cuts,
        both forward and backwards, using the Weyl monodromy.
        """
        g_data = sw_data.g_data
        branch_points = sw_data.branch_points
        irregular_singularities = sw_data.irregular_singularities
        # Adding a minimal radius of 1.0, this is necessary 
        # in case there is only a single branch point at z=0,
        # otherwise max_radius would be 0.
        max_radius = 2 * max(
            [abs(c_l.z) for c_l in branch_points + irregular_singularities] + 
            [1.0]
        )

        # branching will occur at branch points or irregular singularities
        branching_loci = branch_points + irregular_singularities
        br_loc_zs_r = [bl.z.real for bl in branching_loci]
        
        # If the length of the S-wall's coordinates
        # is 3 or less, do not check cuts.
        # Otherwise, the interpolation method would
        # raise an error.
        if len(self.z) > 3:
            # parametrizing the z-coordinate of the k-wall's coordinates
            # as a function of proper time
            traj_t = numpy.array(range(len(self.z)))
            traj_z_r = numpy.array([z.real for z in self.z])
            
            # Scan over branching loci cuts, see if path ever crosses one 
            # based on x-coordinates only
            _cuts_intersections = []
            for br_loc_idx, x_0 in list(enumerate(br_loc_zs_r)):
                g = interpolate.splrep(traj_t, traj_z_r - x_0, s=0)
                # now produce a list of integers corresponding to points in  
                # the S-wall's coordinate list that seem to cross branch-cuts
                # based on the z-coordinate's real part.
                # Now get the list of intersection points as integer values 
                # of the proper time, approximated by the floor function

                __intersection_ts = map(
                    int, map(floor, interpolate.sproot(g))
                )
                # Remove duplicates.
                _intersection_ts = delete_duplicates(__intersection_ts)
                # Enforce imaginary-part of z-coordinate intersection 
                # criterion: branch cuts extend vertically.
                y_0 = branching_loci[br_loc_idx].z.imag
                intersection_ts = [t for t in _intersection_ts
                                   if self.z[t].imag > y_0]
                
                intersections = []
                for t in intersection_ts:
                    branch_locus = branching_loci[br_loc_idx]
                    
                    if branch_locus.__class__.__name__ == 'BranchPoint':
                        # # Drop intersections of a primary S-wall with the 
                        # # branch cut emanating from its parent branch-point
                        # # if such intersections happens at t=0 or t=1.
                        # if (branch_locus.label == self.parents[0] 
                        #                             and (t == 0 or t == 1)):
                        #     continue

                        # Drop intersections of a primary S-wall with the 
                        # branch cut emanating from its parent branch-point
                        # if such intersections happens within a short 
                        # distance from the starting point
                        if (
                            branch_locus.label == self.parents[0] and 
                            # t<6):
                            abs(branch_locus.z - self.z[t]) 
                            < BRANCH_POINT_RADIUS
                        ):
                            # abs(branch_locus.z - self.z[t]) < accuracy):
                            continue

                    # Check that the intersection actually happens
                    # and is not an artifact of the interpolation used above
                    # which could become an extrapolation

                    if not (
                        (
                            self.z[t - 1].real < branch_locus.z.real and
                            branch_locus.z.real < self.z[t + 1].real
                        )
                        or 
                        (
                            self.z[t + 1].real < branch_locus.z.real and
                            branch_locus.z.real < self.z[t - 1].real
                        )
                    ):
                        logging.info('Dropping a fake cut intersection.')
                        continue

                    # Add 
                    # [the branch-point identifier(index), t,
                    #  the direction (either 'cw' or 'ccw')]
                    # to each intersection.
                    intersections.append(
                        [branch_locus.label, t, clock(left_right(self.z, t))]
                    )
                _cuts_intersections += intersections

            # Now sort intersections according to where they happen in proper 
            # time; recall that the elements of cuts_intersections are 
            # organized  as      [..., [branch_point_idx, t, 'ccw'] ,...]
            # where 't' is the integer of proper time at the intersection.
            self.cuts_intersections = sorted(
                _cuts_intersections, 
                cmp=lambda k1, k2: cmp(k1[1], k2[1])
            )
            logging.debug(
                'S-wall {} intersects the following '
                'cuts at the points\n{}.'.format(
                    self.label, self.cuts_intersections
                )
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
            t_0 = sorted(
                ([
                    [t_i, min(
                        z_r_distance_from_ramification_loci(z_i, sw_data)
                    )]
                    for t_i, z_i in enumerate(self.z) if (
                        abs(z_i) < max_radius or t_i == 0
                    )
                ]), cmp=lambda a, b: cmp(a[1], b[1])
            )[-1][0]
        
        else:
            self.cuts_intersections = []
            t_0 = 0

        z_0 = self.z[t_0]
        xs_0 = self.x[t_0]
        self.root_basepoint = [t_0, z_0, xs_0]

        # Determine the initial root-type
        initial_root = get_s_wall_root(z_0, xs_0, sw_data,)
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
                br_loc_label, t, direction = intersections_after_t_0[k]
                branch_locus = branch_locus_from_label(sw_data, br_loc_label)

                current_root = self.local_roots[-1]
                new_root = g_data.weyl_monodromy(
                    current_root, branch_locus, direction
                )
                new_weight_pairs = g_data.ordered_weight_pairs(new_root)

                self.local_roots.append(new_root)
                self.local_weight_pairs.append(new_weight_pairs)

            # Fill in the root types that occur before the basepoint
            # recall that their time-ordering has already been reversed
            # so the first one in the list is the closest to t_0, and so on
            for k in range(len(intersections_before_t_0)):
                br_loc_label, t, direction = intersections_before_t_0[k]
                branch_locus = branch_locus_from_label(sw_data, br_loc_label)

                current_root = self.local_roots[0]
                new_root = g_data.weyl_monodromy(
                    current_root, branch_locus, direction, reverse=True
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
            closed_splittings = self.get_splittings() + [len(self.z) - 1]
            for i, sp in enumerate(closed_splittings):
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
            closed_splittings = self.get_splittings() + [len(self.z) - 1]
            for i, sp in enumerate(closed_splittings):
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
            br_loc_label, t, chi = int_data
            br_loc = branch_locus_from_label(sw_data, br_loc_label)

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
            br_loc_label, t_old, chi = int_data
            t_new = t_old + i + 1
            new_cuts_intersections.append([br_loc_label, t_new, chi])
        
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


def get_s_wall_seeds(sw, theta, branch_point, config,):
    """
    S-walls are seeded from branch points.
    Each branch point has a number of ramification 
    points lying above it.
    Regardless of the representation, it is sufficient
    to consider one of these ramification points
    to extract the seed data.
    We thus stick to (any)one ramification point of the 
    fundamental representation to get the seeds.
    """

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
        sw_diff_coeff = rp.sw_diff_coeff
        logging.debug('Analyze ramification point (z,x)={}'.format([z_0, x_0]))
        logging.debug('Ramification index = {}'.format(r_i))
        logging.debug('Ramification type = {}'.format(rp_type))
        logging.debug(
            'leading coefficient of SW diff = {}\n'.format(sw_diff_coeff)
        )

        # Construct the seeding points for the branch point
        # by studying the type of ramification structure of the r.p.
        if rp_type == 'type_I':
            phases = [exp(2 * pi * 1j * float(i) / r_i) for i in range(r_i)]
            phi = [[p1 - p2 for p1 in phases] for p2 in phases]
            
            omega = exp(2.0 * pi * 1j * float(r_i) / float(r_i + 1))

            dz_phases = ([
                (1.0 / cpow(sw_diff_coeff, r_i, r_i + 1)) *
                exp(1j * theta * float(r_i) / (r_i + 1)) *
                ((1.0 / phi[i][j]) ** (float(r_i) / (r_i + 1))) * (omega ** s)
                for i in range(r_i) for j in range(r_i) 
                for s in range(r_i + 1) if i != j
            ])

            norm_dz_phases = [d / abs(d) for d in dz_phases]
            # these are the normalized phases of the seeds
            # with respect to the branch point:
            zetas = remove_duplicate(
                norm_dz_phases, lambda p1, p2: abs(p1 - p2) < (accuracy)
            )
        
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
            
            omega = exp(2.0 * pi * 1j * float(2 * r_k) / float(2 * r_k + 1))
            dz_phases = ([
                (1.0 / cpow(sw_diff_coeff, 2 * r_k, 2 * r_k + 1)) *
                exp(1j * theta * float(2 * r_k) / (2 * r_k + 1)) *
                ((1.0 / phi[i][j]) ** (float(2 * r_k) / (2 * r_k + 1))) 
                * (omega ** s)
                for i in range(r_k) for j in range(r_k) 
                for s in range(2 * r_k + 1) if i != j
            ] + [
                (1.0 / cpow(sw_diff_coeff, 2 * r_k, 2 * r_k + 1)) *
                exp(1j * theta * float(2 * r_k) / (2 * r_k + 1)) *
                ((1.0 / psi[i][j]) ** (float(2 * r_k) / (2 * r_k + 1))) 
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
            dz_phases = ([
                (1.0 / cpow(sw_diff_coeff, 2 * r_k - 2, 2 * r_k - 1)) *
                exp(1j * theta * float(2 * r_k - 2) / (2 * r_k - 1)) *
                ((1.0 / phi[i][j]) ** (float(2 * r_k - 2) / (2 * r_k - 1))) * 
                (omega ** s)
                for i in range(r_k) for j in range(r_k) 
                for s in range(2 * r_k - 1) if i != j
            ] + [
                (1.0 / cpow(sw_diff_coeff, 2 * r_k - 2, 2 * r_k - 1)) *
                exp(1j * theta * float(2 * r_k - 2) / (2 * r_k - 1)) *
                ((1.0 / psi[i][j]) ** (float(2 * r_k - 2) / (2 * r_k - 1))) * 
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
        
        # Now for each seeding point z_1 we identify two sheets
        # of the cover which match the phase of the displacement z_1-z_0

        for zeta in zetas:
            raise_precision = True
            precision_level = 0
            while (
                raise_precision is True and 
                precision_level <= SEED_PRECISION_MAX_DEPTH
            ):
                logging.debug(
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
                                v_i = complex(
                                    sw.diff.num_v.subs([(z, z_1), (x, x_i)])
                                )
                                v_j = complex(
                                    sw.diff.num_v.subs([(z, z_1), (x, x_j)])
                                )
                                ij_factor = (
                                    -1.0 * exp(1j * theta) / (v_j - v_i)
                                )
                                x_i_x_j_phases.append(
                                    [(ij_factor) / abs(ij_factor), [x_i, x_j]]
                                )

                elif rp_type == 'type_II' or rp_type == 'type_III':
                    # we assume that the ramification index is maximal
                    # therefore we ask for all the sheets at z_1.
                    x_s = find_xs_at_z_0(sw, z_1, ffr=True)

                    # order of magnitude of expected separation 
                    # of sheets at z_1
                    dx = abs(sw_diff_coeff) * (dt ** (1.0 / float(r_i)))
                    x_accuracy = min([accuracy, dx])

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
                                v_i = complex(
                                    sw.diff.num_v.subs([(z, z_1), (x, x_i)])
                                )
                                v_j = complex(
                                    sw.diff.num_v.subs([(z, z_1), (x, x_j)])
                                )
                                ij_factor = (
                                    -1.0 * exp(1j * theta) / (v_j - v_i)
                                )
                                # ij_factor = -1.0 * exp(1j*theta)/(x_j - x_i)
                                x_i_x_j_phases.append(
                                    [(ij_factor) / abs(ij_factor), [x_i, x_j]]
                                )

                closest_pair = sorted(
                    x_i_x_j_phases, key=lambda p: abs(p[0] - zeta)
                )[0][1]
                phase_mismatch = abs(sorted(
                    x_i_x_j_phases, key=lambda p: abs(p[0] - zeta)
                )[0][0] - zeta)
                if phase_mismatch < SEED_PHASE_PRECISION:
                    logging.debug(
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
                        logging.warning(
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
    logging.info('Number of S-walls emanating = {}'.format(len(seeds)))
    logging.debug('these are the seeds {}\n'.format(seeds))
    return seeds


def is_root(np_array, g_data):
    ans = False
    for rt in list(g_data.roots):
        if (np_array == rt).all():
            ans = True
            break
        else:
            pass
    return ans


def get_joint(z, s_wall_1, s_wall_2, t_1, t_2, sw_data=None, label=None):
    """
    Return a joint if formed, otherwise return None.
    """
    logging.debug('evaluating possible joint at z = {}'.format(z))
    alpha_1 = s_wall_1.get_root_at_t(t_1)
    alpha_2 = s_wall_2.get_root_at_t(t_2)

    if is_root(alpha_1 + alpha_2, sw_data.g_data):
        return Joint(z, s_wall_1, s_wall_2, t_1, t_2, label, sw_data)
    else:
        return None


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