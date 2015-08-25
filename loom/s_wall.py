import logging
import numpy
import sympy
import pdb

from cmath import exp, pi, phase
from math import floor

from geometry import get_local_sw_diff, find_xs_at_z_0
from misc import (gather, cpow, remove_duplicate, unravel, ctor2, r2toc,
                  GetSWallSeedsError, n_nearest_indices, is_root)

x, z = sympy.symbols('x z')

# TODO: Generalize this.
# Number of x's at a fixed z
num_x_over_z = 2

class Joint:
    def __init__(self, z=None, s_wall_1=None, s_wall_2=None,
                t_1=None, t_2=None, label=None, sw_data=None):
        alpha_1 = s_wall_1.get_root_at_t(t_1)
        alpha_2 = s_wall_2.get_root_at_t(t_2) 

        self.z = z
        self.root = alpha_1 + alpha_2
        self.M = s_wall_1.M[t_1] + s_wall_2.M[t_2]
        self.parents = [s_wall_1.label, s_wall_2.label]
        self.indices = [t_1, t_2]
        self.label = [s_wall_1.label, s_wall_2.label]

        xs_at_z = find_xs_at_z_0(sw_data, z)
        new_wall_weight_pairs = sw_data.g_data.ordered_weight_pairs(self.root)
        w_p_0 = new_wall_weight_pairs[0]
        x_i_s = [xs_at_z[w_p[0]] for w_p in new_wall_weight_pairs]
        x_j_s = [xs_at_z[w_p[1]] for w_p in new_wall_weight_pairs]
        self.x_s = [[x_i_s[k], x_j_s[k]] for k in range(len(x_i_s))]

        ffr_xs_at_z = find_xs_at_z_0(sw_data, z, ffr=True)
        ffr_new_wall_weight_pairs = sw_data.g_data.ordered_weight_pairs(
                                                        self.root, ffr=True)
        ffr_w_p_0 = ffr_new_wall_weight_pairs[0]
        ffr_x_i = ffr_xs_at_z[ffr_w_p_0[0]]
        ffr_x_j = ffr_xs_at_z[ffr_w_p_0[1]]
        self.ffr_x = [ffr_x_i, ffr_x_j]

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
                raise ValueError('Cannot pick ij, jk pairs ' +
                                'from colliding S-walls\n{}\n{}'.format(
                                    weight_pairs_1, weight_pairs_2)
                                )
        self.combos = combos


    def __eq__(self, other):
        return self.label == other.label

    def get_json_data(self):
        json_data = {
            'z': ctor2(self.z),
            'M': self.M,
            'parents': [parent for parent in self.parents],
            'label': self.label,
        }
        return json_data

    def set_json_data(self, json_data):
        self.z = r2toc(json_data['z'])
        self.M = json_data['M']
        self.parents = [parent for parent in json_data['parents']]
        self.label = json_data['label']

    def is_equal_to(self, other, accuracy):
        if(abs(self.z - other.z) > accuracy):
            return False
        if(abs(self.M - other.M) > accuracy):
            return False
        if(len(self.x_s) != len(other.x_s)):
            return False
        for i in range(len(self.x_s)):
            if ((abs(self.x_s[i][0] - other.x_s[i][0]) > accuracy)
                or (abs(self.x_s[i][1] - other.x_s[i][1]) > accuracy)):
                return False
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
            # self.x = numpy.empty(n_steps+1, (complex, num_x_over_z))
            self.x = numpy.empty(n_steps+1, (complex, 2))
            #self.M = numpy.empty(n_steps+1, float)
            self.M = numpy.empty(n_steps+1, complex)
            self.z[0] = z_0
            self.x[0] = x_0
            self.M[0] = M_0
        self.parents = parents
        self.label = label

        self.cuts_intersections = []
        self.splittings = []
        self.local_roots = []
        self.local_weight_pairs = []        

    def __setitem__(self, t, data):
        """
        Set the data of the S-wall at t, where
            data = [z[t], x[t][0], x[t][1], ..., M[t]]
        """
        self.z[t] = data[0]
        self.x[t] = data[1:num_x_over_z+1]
        self.M[t] = data[num_x_over_z+1]


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
        self.x.resize((size, num_x_over_z))
        self.M.resize((size))


    def get_json_data(self):
        json_data = {
            'z': numpy.array([self.z.real, self.z.imag]).T.tolist(),
            'M': numpy.array(self.M).T.tolist(),
            'x': numpy.rollaxis(
                numpy.array([self.x.real, self.x.imag]), 0, 3
            ).tolist(),
            'parents': [parent for parent in self.parents],
            'label': self.label,
        }
        return json_data


    def set_json_data(self, json_data):
        self.z = numpy.array([r2toc(z_t) for z_t in json_data['z']])
        self.M = numpy.array(json_data['M'])
        self.x = numpy.array(
            [[r2toc(x_i) for x_i in x_t] for x_t in json_data['x']]
        )
        self.parents = [parent for parent in json_data['parents']]
        self.label = json_data['label']


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
        z_range_limits = config['z_range_limits']
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

        if z_range_limits is not None:
            [z_real_min, z_real_max], [z_imag_min, z_imag_max] = z_range_limits

        while ode.successful() and step < num_of_steps:
            step += 1
            # Stop if z is inside a cutoff of a puncture.
            if len(ppzs) > 0:
                min_d = min([abs(z_i - ppz) for ppz in ppzs])
                if min_d < size_of_puncture_cutoff:
                    self.resize(step)
                    break

            # Stop if z is ouside the range limit.
            if z_range_limits is not None:
                if (z_i.real < z_real_min or
                    z_i.real > z_real_max or
                    z_i.imag < z_imag_min or
                    z_i.imag > z_imag_max):
                    self.resize(step)
                    break

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
            M_i = y_i[num_x_over_z+1]
            self[step] = y_i

    def get_root_at_t(self, t):
        """
        Given an integer t which parametrizes a point 
        of the trajectory in proper time, return the local 
        root at that point.
        """
        if t < 0 or t > (len(self.z) - 1):
            raise ValueError
        else:
            closed_splittings = self.splittings + [len(self.z) - 1]
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
            closed_splittings = self.splittings + [len(self.z) - 1]
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

    def enhance_at_cuts(self):
        # Add the intersection points of Swalls and branch cuts
        # also update the intersection data accordingly
        # NOTE: we update the .z and .M attributes of the SWall 
        # class, but we don't update the .x attribute, as it seems 
        # unnecessary for now
        wall_pieces_z = []
        wall_pieces_M = []

        # split the wall into piececs, at the end of each 
        # piece add the corresponding intersection point
        old_cuts_intersections = [i for i in self.cuts_intersections]
        for i, int_data in enumerate(old_cuts_intersections):
            bp, t_0, chi = int_data
            point_to_add = get_intermediate_point(
                                        self.z[t_0], self.z[t_0+1], bp.z)
            mass_to_add = get_intermediate_mass(
                                                self.M[t_0], self.M[t_0+1], 
                                                self.z[t_0], self.z[t_0+1], 
                                                point_to_add
                                                )
            if i == 0:
                piece = list(self.z[:t_0+1])
                piece.append(point_to_add)
                mass_piece = list(self.M[:t_0+1])
                mass_piece.append(mass_to_add)
            else:
                t_0_prev = old_cuts_intersections[i-1][1]
                piece = list(self.z[t_0_prev+1:t_0+1])
                piece.append(point_to_add)
                mass_piece = list(self.M[t_0_prev+1:t_0+1])
                mass_piece.append(mass_to_add)

            wall_pieces_z.append(piece)
            wall_pieces_M.append(mass_piece)

        # get the last piece of the wall
        last_t_0 = old_cuts_intersections[-1][1]
        last_piece = list(self.z[last_t_0+1:])
        wall_pieces_z.append(last_piece)
        last_mass_piece = list(self.M[last_t_0+1:])
        wall_pieces_M.append(last_mass_piece)

        # now assemble the new pieces
        splittings = []
        new_s_wall_z = []
        new_s_wall_M = []
        
        for piece in wall_pieces_z:
            new_s_wall_z += piece
            #recall that the newly added point is at
            # the end of each piece
            splittings.append(len(new_s_wall_z) - 1)
        # now remove the last splitting, because at the end 
        # of the last piece there is no actual branch cut
        self.splittings = splittings[:len(splittings)-1]

        for mass_piece in wall_pieces_M:
            new_s_wall_M += mass_piece

        # update the z coordinates
        self.z = numpy.array(new_s_wall_z)

        # update the mass
        self.M = numpy.array(new_s_wall_M)

        # update the intersection data
        new_cuts_intersections = []
        for i, int_data in enumerate(old_cuts_intersections):
            bp, t_old, chi = int_data
            t_new = self.splittings[i]
            new_cuts_intersections.append([bp, t_new, chi])
        
        self.cuts_intersections = new_cuts_intersections

        pass






def get_intermediate_point(z_1, z_2, z_med):
    """
    get the intermediate point between z_1 and z_2
    in correspondence of the real part of z_med
    """
    x_1 = z_1.real
    y_1 = z_1.imag
    x_2 = z_2.real
    y_2 = z_2.imag
    x_med = z_med.real
    # if not x_2 < x_med < x_1 or x_1 < x_med < x_2:
    #     print 'ERROR: x_1={}, x_2={}, x_med={}'.format(x_1,x_2,x_med)
    slope = (y_2 - y_1) / (x_2 - x_1)
    y_med = y_1 + slope * (x_med - x_1)
    return x_med + 1j * y_med

def get_intermediate_mass(m_1, m_2, z_1, z_2, z_med):
    """
    get the intermediate mass between m_1 and m_2
    in correspondence of z_med
    """
    # if not x_2 < x_med < x_1 or x_1 < x_med < x_2:
    #     print 'ERROR: x_1={}, x_2={}, x_med={}'.format(x_1,x_2,x_med)
    slope = (m_2 - m_1) / (z_2 - z_1)
    m_med = m_1 + slope * (z_med - z_1)
    return m_med


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

    cs = gather(cs, lambda c1, c2: abs(c1 - c2) < delta)
    logging.debug('list of c = %s, # = %d', cs, len(cs))

    # 2. Now calculate \Delta z_i for each S-wall and
    # find the two points on the curve that are projected onto it.
    seeds = []
    for c in cs:
        cv = c[0]   # value of c
        cm = c[1]   # multiplicity of c
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


# def get_joint(z, x1_i, x2_i, x1_j, x2_j, M1, M2, parent_i, parent_j,
#               accuracy=None, xs_at_z=None, g_data=None, label=None):
#     """
#     Return a joint if formed, otherwise return None.
#     """

#     if (abs(x1_i - x2_j) < accuracy and abs(x1_j - x2_i) < accuracy):
#         return None
#     elif (abs(x2_i - x1_j) < accuracy):
#         if differ_by_root(
#             x1_i, x2_j, accuracy=accuracy, xs=xs_at_z, g_data=g_data,
#         ):
#             return Joint(z, [x1_i, x2_j], M1+M2, [parent_i, parent_j], label)
#         else:
#             return None
#     elif (abs(x2_j - x1_i) < accuracy):
#         if differ_by_root(
#             x1_j, x2_i, accuracy=accuracy, xs=xs_at_z, g_data=g_data,
#         ):
#             return Joint(z, [x1_j, x2_i], M1+M2, [parent_j, parent_i], label)
#         else:
#             return None

def get_joint(z, s_wall_1, s_wall_2, t_1, t_2, 
              sw_data=None, label=None):
    """
    Return a joint if formed, otherwise return None.
    """

    alpha_1 = s_wall_1.get_root_at_t(t_1)
    alpha_2 = s_wall_2.get_root_at_t(t_2)

    if is_root(alpha_1 + alpha_2, sw_data.g_data):
        return Joint(z, s_wall_1, s_wall_2, t_1, t_2, label, sw_data)
    else:
        return None


def differ_by_root(x1, x2, accuracy=None, xs=None, g_data=None):
    root_system = g_data.root_system
    k = g_data.fundamental_representation_index
    # NOTE: A shortcut. Use this after checking this function
    # works correctly.
    if root_system[0] == 'A':
        if k == 1:
            return True
    elif root_system[0] == 'D' and k == 1:
        if abs(x1-(-x2)) < accuracy:
            return False
        else:
            return True
    else:
        # TODO
        raise NotImplementedError
        
    

