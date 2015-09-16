import logging
import numpy
import sympy
import pdb

from cmath import exp, pi, phase
from math import floor
from scipy import interpolate

from geometry import get_local_sw_diff, find_xs_at_z_0
from misc import (gather, cpow, remove_duplicate, unravel, ctor2, r2toc,
                  GetSWallSeedsError, n_nearest_indices, delete_duplicates,
                  clock, left_right,)

x, z = sympy.symbols('x z')

# TODO: Generalize this.
# Number of x's at a fixed z
NUM_ODE_XS_OVER_Z = 2

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
        #self.indices = [t_1, t_2]

        # FIXME: Erase the following if we decide not to use
        # x-coords in calculating joints of S-walls.
        # Probably this will not be additional information because
        # x's are calculated from the x's of the ffr cover.
        #xs_at_z = find_xs_at_z_0(sw_data, z)
        new_wall_weight_pairs = sw_data.g_data.ordered_weight_pairs(self.root)
        ##w_p_0 = new_wall_weight_pairs[0]
        #x_i_s = [xs_at_z[w_p[0]] for w_p in new_wall_weight_pairs]
        #x_j_s = [xs_at_z[w_p[1]] for w_p in new_wall_weight_pairs]
        #self.x_s = [[x_i_s[k], x_j_s[k]] for k in range(len(x_i_s))]

        # FIXME: The following, including self.ode_xs, can be removed
        # once the seeding of an S-wall is done by using a root.
        ffr_xs_at_z = find_xs_at_z_0(sw_data, z, ffr=True)
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
            ij_in_1 = [pair for pair in weight_pairs_1 if (pair[0]==w_i)]
            ij_in_2 = [pair for pair in weight_pairs_2 if (pair[0]==w_i)]
            jk_in_1 = [pair for pair in weight_pairs_1 if (pair[1]==w_k)]
            jk_in_2 = [pair for pair in weight_pairs_2 if (pair[1]==w_k)]
            if len(ij_in_1)>0 and len(jk_in_2)>0:
                combos.append([ij_in_1[0], jk_in_2[0]])
            elif len(ij_in_2)>0 and len(jk_in_1)>0:
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
        #if(len(self.x_s) != len(other.x_s)):
        #    return False
        #for i in range(len(self.x_s)):
        #    if ((abs(self.x_s[i][0] - other.x_s[i][0]) > accuracy)
        #        or (abs(self.x_s[i][1] - other.x_s[i][1]) > accuracy)):
        #        return False
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
        ### FIXME: self.zs & self.xs instead of z & x?
        if n_steps is None:
            self.z = []
            self.x = []
            self.M = []
        else:
            self.z = numpy.empty(n_steps+1, complex)
            self.x = numpy.empty(n_steps+1, (complex, NUM_ODE_XS_OVER_Z))
            self.M = numpy.empty(n_steps+1, complex)
            self.z[0] = z_0
            self.x[0] = x_0
            self.M[0] = M_0
        self.parents = parents
        self.label = label

        # cuts_intersections = [[b_pt_idx, i, '(c)cw'], ...]
        self.cuts_intersections = []
        self.local_roots = []
        # local_weight_pairs is a list of pair of intgers.
        self.local_weight_pairs = []        


    def __setitem__(self, t, data):
        """
        Set the data of the S-wall at t, where
            data = [z[t], x[t][0], x[t][1], ..., M[t]]
        """
        self.z[t] = data[0]
        self.x[t] = data[1:NUM_ODE_XS_OVER_Z+1]
        self.M[t] = data[NUM_ODE_XS_OVER_Z+1]


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

        for t in range(1, len(self.z)-1):
            x_1 = self.z[t].real
            y_1 = self.z[t].imag
            x_2 = self.z[t+1].real
            y_2 = self.z[t+1].imag
            if ((x_1 - x_0) * (x_2 - x_1) < 0 or
                (y_1 - y_0) * (y_2 - y_1) < 0):
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
    ):
        rpzs = ramification_point_zs
        ppzs = puncture_point_zs
        #z_range_limits = config['z_range_limits']
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

        #if z_range_limits is not None:
        #    [z_real_min, z_real_max], [z_imag_min, z_imag_max] = z_range_limits

        while ode.successful() and step < num_of_steps:
            step += 1
            # Stop if z is inside a cutoff of a puncture.
            if len(ppzs) > 0:
                min_d = min([abs(z_i - ppz) for ppz in ppzs])
                if min_d < size_of_puncture_cutoff:
                    self.resize(step)
                    break

            # Stop if z is ouside the range limit.
            #if z_range_limits is not None:
            #    if (z_i.real < z_real_min or
            #        z_i.real > z_real_max or
            #        z_i.imag < z_imag_min or
            #        z_i.imag > z_imag_max):
            #        self.resize(step)
            #        break

            # Stop if M exceeds mass limit.
            if mass_limit is not None:
                if M_i > mass_limit:
                    self.resize(step)
                    break

            # Adjust the step size if z is near a branch point.
            if (len(rpzs) > 0 and
                min([abs(z_i - rpz) for rpz in rpzs]) < size_of_neighborhood):
                dt = size_of_small_step
            else:
                dt = size_of_large_step

            y_i = ode.integrate(ode.t + dt)
            z_i = y_i[0]
            M_i = y_i[NUM_ODE_XS_OVER_Z+1]
            self[step] = y_i


    def determine_root_types(self, sw_data):

        """
        Determine at which points the wall crosses a cut, 
        for instance [55, 107, 231] would mean that 
        it changes root-type 3 times. 
        """
        g_data = sw_data.g_data
        branch_points = sw_data.branch_points

        # Determine the initial root-type
        z_0 = self.z[0]
        xs_0 = self.x[0]
        initial_root = get_s_wall_root(z_0, xs_0, sw_data,)
        # A list of ordered pairs [...[i, j]...]
        # such that weights[j] - weights[i] = root
        initial_weight_pairs = (
            sw_data.g_data.ordered_weight_pairs(initial_root,)
        )
        
        self.local_roots = [initial_root]
        self.local_weight_pairs = [initial_weight_pairs]

        # If the length of the S-wall's coordinates
        # is 3 or less, do not check cuts.
        # Otherwise, the interpolation methood would
        # raise an error.
        if len(self.z) <= 3:
            return None
        bpzs_r = [bp.z.real for bp in branch_points]
        
        # parametrizing the z-coordinate of the k-wall's coordinates
        # as a function of proper time
        traj_t = numpy.array(range(len(self.z)))
        traj_z_r = numpy.array([z.real for z in self.z])
        
        # Scan over branch cuts, see if path ever crosses one 
        # based on x-coordinates only
        for b_pt_idx, x_0 in list(enumerate(bpzs_r)):
            g = interpolate.splrep(traj_t, traj_z_r - x_0, s=0)
            # now produce a list of integers corresponding to points in the 
            # S-wall's coordinate list that seem to cross branch-cuts
            # based on the z-coordinate's real part.
            # Now get the list of intersection points as integer values 
            # of the proper time, approximated by the floor function

            __intersection_ts = map(int, map(floor, interpolate.sproot(g)))
            # Remove duplicates.
            _intersection_ts = delete_duplicates(__intersection_ts)
            # Enforce imaginary-part of z-coordinate intersection criterion:
            # branch cuts extend vertically.
            y_0 = branch_points[b_pt_idx].z.imag
            intersection_ts = [t for t in _intersection_ts
                               if self.z[t].imag > y_0]
            
            intersections = []
            for t in intersection_ts:
                bp = branch_points[b_pt_idx]
                # Drop intersections of a primary S-wall with the 
                # branch cut emanating from its parent branch-point
                # if such intersections happens at t=0 or t=1.
                if bp.label == self.parents[0] or t == 0 or t == 1:
                    continue
                # Add 
                # [the branch-point identifier(index), t,
                #  the direction (either 'cw' or 'ccw')]
                # to each intersection.
                intersections.append(
                    [b_pt_idx, t, clock(left_right(self.z, t))]
                )
            self.cuts_intersections += intersections

        # TODO: Might be worth implementing an algorithm for handling 
        # overlapping branch cuts: e.g. the one with a lower starting point 
        # will be taken to be on the left, or a similar criterion.
        # Still, there will be other sorts of problems, it is necessary
        # to just rotate the z-plane and avoid such situations.

        # Now sort intersections according to where they happen in proper 
        # time; recall that the elements of cuts_intersections are organized 
        # as      [..., [branch_point_idx, t, 'ccw'] ,...]
        # where 't' is the integer of proper time at the intersection.
        self.cuts_intersections = sorted(
            self.cuts_intersections, 
            cmp = lambda k1, k2: cmp(k1[1], k2[1])
        )
        logging.debug(
            'S-wall {} intersects the following'
            'cuts at the points\n{}.'.format(self.label, intersections)
        )

        if len(self.cuts_intersections) > 0:
            # Add the actual intersection point to the S-wall
            # then update the attribute SWall.cuts_intersections accordingly
            self.enhance_at_cuts(branch_points)
            
            for k in range(len(self.cuts_intersections)):
                bp_idx, t, direction = self.cuts_intersections[k]
                branch_point = branch_points[bp_idx]

                current_root = self.local_roots[-1]
                new_root = g_data.weyl_monodromy(
                    current_root, branch_point, direction
                )
                new_weight_pairs = g_data.ordered_weight_pairs(new_root)

                self.local_roots.append(new_root)
                self.local_weight_pairs.append(new_weight_pairs)



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


    def enhance_at_cuts(self, branch_points):
        # Add the intersection points of Swalls and branch cuts
        # also update the intersection data accordingly
        wall_pieces_z = []
        wall_pieces_x = []
        wall_pieces_M = []

        # split the wall into piececs, at the end of each 
        # piece add the corresponding intersection point
        t_0 = 0
        for int_data in self.cuts_intersections:
            bp_idx, t, chi = int_data
            bp = branch_points[bp_idx]

            z_1 = self.z[t]
            z_2 = self.z[t+1]
            z_to_add = get_intermediate_z_point(z_1, z_2, bp.z)

            xs_1 = self.x[t]
            xs_2 = self.x[t+1]
            xs_to_add = [
                get_intermediate_value(xs_1[i], xs_2[i], z_1, z_2, z_to_add)
                for i in range(NUM_ODE_XS_OVER_Z)
            ]

            M_to_add = get_intermediate_value(self.M[t], self.M[t+1], 
                                              z_1, z_2, z_to_add)

            z_piece = numpy.concatenate(
                (self.z[t_0:t+1], numpy.array([z_to_add], dtype=complex))
            )
            x_piece = numpy.concatenate(
                (self.x[t_0:t+1], numpy.array([xs_to_add], dtype=complex))
            )
            M_piece = numpy.concatenate(
                (self.M[t_0:t+1], numpy.array([M_to_add], dtype=complex))
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
            bp, t_old, chi = int_data
            t_new = t_old + i + 1
            new_cuts_intersections.append([bp, t_new, chi])
        
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
    slope = (v_2 - v_1) / (z_2 - z_1)
    v_med = v_1 + slope * (z_med - z_1)
    return v_med


def get_s_wall_seeds(sw, theta, branch_point, config,):
    ### S-walls are seeded from branch points.
    ### Each branch point has a number of ramification 
    ### points lying above it.
    ### Regardless of the representation, it is sufficient
    ### to consider one of these ramification points
    ### to extract the seed data.
    ### We thus stick to (any)one ramification point of the 
    ### fundamental representation to get the seeds.
    rp = branch_point.ffr_ramification_points[0]
    delta = config['accuracy']
    dt = config['size_of_small_step']

    ###
    # 1. Find the first-order approximations of the starting points
    # of S-walls around a given branch point, which is of the form
    # \Delta z_i = c_i / (\lambda_0)^(rp.i/(rp.i+1))
    # at a branch point from a ramification point, and
    # \Delta z_i = c_i / (\lambda_0)^(rp.i) exp(rp.i theta I)
    # at a branch point from a massless regular puncture
    ###

    # 1.1 find the coefficient and the exponent of the leading term
    # of the SW differential at the ramification point.
    lambda_0, diff_e = get_local_sw_diff(sw, rp, ffr=True)

    # 1.2 find c_i, a phase factor for each S-wall.
    omega_1 = exp(2*pi*1j/rp.i)
    omega = [omega_1**k for k in range(rp.i)]

    beta_1 = exp(rp.i*2*pi*1j/(rp.i+1))
    beta = [beta_1**k for k in range(rp.i+1)]

    cs = []
    if diff_e == sympy.Rational(1, rp.i):
        # the branch point is a ramification point
        # go over pairs of omegas that differ by \omega_1^i
        for i in range(1, rp.i):
            new_locs = []
            # and go over all the omegas
            for j in range(rp.i):
                if j+i < rp.i:
                    new_loc = 1/(omega[j]-omega[j+i])
                else:
                    new_loc = 1/(omega[j]-omega[j+i-rp.i])
                new_loc = cpow(new_loc, rp.i, rp.i+1)
                new_locs += [new_loc*beta_i for beta_i in beta]
            new_locs = remove_duplicate(new_locs,
                                        lambda l1, l2: abs(l1 - l2) < delta)
            cs += [(c*exp(theta*1j*sympy.Rational(rp.i, rp.i+1)) /
                    cpow(lambda_0, rp.i, rp.i+1)) for c in new_locs]

    elif diff_e == -1 + sympy.Rational(1, rp.i):
        # the branch point is a massless regular puncture
        # go over pairs of omegas that differ by \omega_1^i
        for i in range(1, rp.i):
            cs.append(exp(rp.i*theta*1j) /
                      cpow(((omega[0]-omega[i])*lambda_0), rp.i))

    else:
        logging.error('unknown form of sw.diff at rp ({}, {}): '
                      'diff_e  = {}'.format(rp.z, rp.x, diff_e))
        raise GetSWallSeedsError(diff_e)

    gathered_cs = gather(cs, lambda c1, c2: abs(c1 - c2) < delta)
    logging.debug('list of c = %s, # = %d', gathered_cs, len(gathered_cs))

    # 2. Now calculate \Delta z_i for each S-wall and
    # find the two points on the curve that are projected onto it.
    seeds = []
    for cv, cs in gathered_cs.iteritems():
        # cv is the value of c, cm is its multiplicity.
        cm = len(cs)
        # resize to the size of the small step
        Delta_z = cv/abs(cv)*delta
        z_0 = rp.z + Delta_z
        xs_at_z_0 = find_xs_at_z_0(sw, z_0, rp.x, rp.i, ffr=True)
        dev_phases = [pi for i in range(len(xs_at_z_0)**2)] 
        for i in range(len(xs_at_z_0)):
            diffx = sw.diff.num_v.subs(z, z_0)
            v_i = complex(diffx.subs(x, xs_at_z_0[i]))
            for j in range(len(xs_at_z_0)):
                if i == j:
                    continue
                else:
                    v_j = complex(diffx.subs(x, xs_at_z_0[j]))
                    delta_z = exp(1j*theta)/(v_i - v_j)*dt
                    # flattened index
                    fij = i*len(xs_at_z_0) + j
                    dev_phases[fij] = phase((delta_z/Delta_z))
        min_dev_indices = n_nearest_indices(dev_phases, 0.0, cm)
        for k in min_dev_indices:
            i, j = unravel(k, len(xs_at_z_0))
            M_0 = 0
            seeds.append(
                [z_0, [xs_at_z_0[i], xs_at_z_0[j]], M_0]
            )

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

    alpha_1 = s_wall_1.get_root_at_t(t_1)
    alpha_2 = s_wall_2.get_root_at_t(t_2)

    if is_root(alpha_1 + alpha_2, sw_data.g_data):
        return Joint(z, s_wall_1, s_wall_2, t_1, t_2, label, sw_data)
    else:
        return None

