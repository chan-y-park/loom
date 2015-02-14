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
    def __init__(self, z=None, x1=None, x2=None, parents=None, label=None,):
        self.z = z
        self.x1 = x1
        self.x2 = x2
        self.parents = parents 
        self.label = label

    def __eq__(self, other):
        return self.label == other.label

    def get_json_data(self):
        json_data = {
            'z': ctor2(self.z),
            'x1': ctor2(self.x1),
            'x2': ctor2(self.x2),
            'parents': [parent for parent in self.parents],
            'label': self.label,
        }
        return json_data

    def set_json_data(self, json_data):
        self.z = r2toc(json_data['z'])
        self.x1 = r2toc(json_data['x1'])
        self.x2 = r2toc(json_data['x2'])
        self.parents = [parent for parent in json_data['parents']]
        self.label = json_data['label']

    def is_equal_to(self, other, accuracy):
        if(abs(self.z - other.z) < accuracy and
           abs(self.x1 - other.x1) < accuracy and
           abs(self.x2 - other.x2) < accuracy):
            return True
        else:
            return False


class SWall(object):
    def __init__(self, z_0=None, x_0=None, parents=None,
                 label=None, n_steps=None,):
        """
        SWall.data is a NumPy array of [z_i, x_i] with length n_steps+1,
        where z_i is the base coordinate and x_i is a Numpy array of 
        the fiber coordinates at t = t_i, i.e. 
            SWall.data[i] = [z_i, [x_i[0], ...]].
        """
        if n_steps is None:
            self.z = []
            self.x = []
        else:
            self.z = numpy.empty(n_steps+1, complex)
            self.x = numpy.empty(n_steps+1, (complex, num_x_over_z))
            self.z[0] = z_0
            self.x[0] = x_0
        self.parents = parents
        self.label = label


    def __setitem__(self, t, data):
        """
        Set the data of the S-wall at t, where data = [z, x[0], ...]
        """
        self.z[t] = data[0]
        self.x[t] = data[1:]


    def __getitem__(self, t):
        """
        Get the data of the S-wall at t, where data = [z, x[0], ...]
        """
        return numpy.concatenate([[self.z[t]], self.x[t]])


    def resize(self, size):
        """
        Resize z & x arrays to discard garbage data.
        """
        self.z.resize((size))
        self.x.resize((size, num_x_over_z))


    def get_json_data(self):
        json_data = {
            'z': [ctor2(z_t) for z_t in self.z],
            'x': [[ctor2(x_i) for x_i in x_t] for x_t in self.x],
            'parents': [parent for parent in self.parents],
            'label': self.label,
        }
        return json_data


    def set_json_data(self, json_data):
        self.z = numpy.array([r2toc(z_t) for z_t in json_data['z']]) 
        self.x = numpy.array(
            [[r2toc(x_i) for x_i in x_t] for x_t in json_data['x']]
        )
        self.parents = [parent for parent in json_data['parents']]
        self.label = json_data['label']
    

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

        step = 0
        z_i = self.z[0]
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

            # Adjust the step size if z is near a branch point.
            if (len(rpzs) > 0 and 
                min([abs(z_i - rpz) for rpz in rpzs]) < size_of_neighborhood):
                dt = size_of_small_step
            else:
                dt = size_of_large_step

            y_i = ode.integrate(ode.t + dt)
            z_i = y_i[0]
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
            seeds.append(
                [z_0, [xs_at_z_0[i], xs_at_z_0[j]]]
            )
        
    return seeds


def get_joint(z, x1_i, x2_i, x1_j, x2_j, parent_i, parent_j, accuracy, 
              root_system=None, label=None):
    """
    Return a joint if formed, otherwise return None.
    """
    if (abs(x1_i - x2_j) < accuracy and abs(x1_j - x2_i) < accuracy):
        return None
    elif (abs(x2_i - x1_j) < accuracy):
        if (root_system[0] == 'D' and abs(x1_i-(-x2_j)) < accuracy):
            return None
        else:
            return Joint(z, x1_i, x2_j, [parent_i, parent_j], label)
    elif (abs(x2_j - x1_i) < accuracy):
        if (root_system[0] == 'D' and abs(x1_j-(-x2_i)) < accuracy):
            return None
        else:
            return Joint(z, x1_j, x2_i, [parent_j, parent_i], label)


