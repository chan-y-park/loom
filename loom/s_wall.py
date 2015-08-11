import logging
import numpy
import sympy
import pdb

from cmath import exp, pi, phase

from geometry import get_local_sw_diff
from misc import (gather, cpow, remove_duplicate, unravel, ctor2, r2toc,
                  GetSWallSeedsError, n_nearest_indices, find_xs_at_z_0,)

x, z = sympy.symbols('x z')

# Number of x's at a fixed z
num_x_over_z = 2

class Joint:
    def __init__(self, z=None, x=None, M=None, parents=None, label=None,):
        self.z = z
        self.x = x
        self.M = M
        self.parents = parents
        self.label = label

    def __eq__(self, other):
        return self.label == other.label

    def get_json_data(self):
        json_data = {
            'z': ctor2(self.z),
            'x': [ctor2(x_i) for x_i in self.x],
            'M': self.M,
            'parents': [parent for parent in self.parents],
            'label': self.label,
        }
        return json_data

    def set_json_data(self, json_data):
        self.z = r2toc(json_data['z'])
        self.x = [r2toc(x_i) for x_i in json_data['x']]
        self.M = json_data['M']
        self.parents = [parent for parent in json_data['parents']]
        self.label = json_data['label']

    def is_equal_to(self, other, accuracy):
        if(abs(self.z - other.z) > accuracy):
            return False
        if(abs(self.M - other.M) > accuracy):
            return False
        if(len(self.x) != len(other.x)):
            return False
        for i in range(len(self.x)):
            if (abs(self.x[i] - other.x[i]) > accuracy):
                return False
        return True


class SWall(object):
    def __init__(self, z_0=None, x_0=None, M_0=None, parents=None,
                 label=None, n_steps=None,):
        """
        SWall.z is a NumPy array of length n_steps+1,
        where z[t] is the base coordinate.

        SWall.x is a Numpy array of the fiber coordinates at t, i.e.
            SWall.x[t] = [x[t][0], x[t][1], ...].
        """
        if n_steps is None:
            self.z = []
            self.x = []
            self.M = []
        else:
            self.z = numpy.empty(n_steps+1, complex)
            self.x = numpy.empty(n_steps+1, (complex, num_x_over_z))
            self.M = numpy.empty(n_steps+1, float)
            self.z[0] = z_0
            self.x[0] = x_0
            self.M[0] = M_0
        self.parents = parents
        self.label = label
        # XXX: interface for marking branch-cut crossings.
        self.splitting = []


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
            z_real_min, z_real_max, z_imag_min, z_imag_max = z_range_limits

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


def get_s_wall_seeds(sw, theta, ramification_point, config,):
    rp = ramification_point
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
    lambda_0, diff_e = get_local_sw_diff(sw, rp)

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
        xs_at_z_0 = find_xs_at_z_0(sw.curve.num_eq, z_0, rp.x, rp.i)
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


def get_joint(z, x1_i, x2_i, x1_j, x2_j, M1, M2, parent_i, parent_j, accuracy=None,
              xs_at_z=None, g_data=None, label=None):
    """
    Return a joint if formed, otherwise return None.
    """

    if (abs(x1_i - x2_j) < accuracy and abs(x1_j - x2_i) < accuracy):
        return None
    elif (abs(x2_i - x1_j) < accuracy):
        if differ_by_root(
            x1_i, x2_j, accuracy=accuracy, xs=xs_at_z, g_data=g_data,
        ):
            return Joint(z, [x1_i, x2_j], M1+M2, [parent_i, parent_j], label)
        else:
            return None
    elif (abs(x2_j - x1_i) < accuracy):
        if differ_by_root(
            x1_j, x2_i, accuracy=accuracy, xs=xs_at_z, g_data=g_data,
        ):
            return Joint(z, [x1_j, x2_i], M1+M2, [parent_j, parent_i], label)
        else:
            return None


def differ_by_root(x1, x2, accuracy=None, xs=None, g_data=None):
    root_system = g_data["root_system"]
    k = g_data["representation"]
    # NOTE: A shortcut. Use this after checking this function
    # works correctly.
    if root_system[0] == 'A':
        if k == 1:
            return True
        else:
            # XXX: This part is not implemented yet.
            return False
    elif root_system[0] == 'D' and k == 1:
        if abs(x1-(-x2)) < accuracy:
            return False
        else:
            return True


