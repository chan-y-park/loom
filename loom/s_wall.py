import logging
import sympy
import pdb

from cmath import exp, pi, phase

from curve import get_local_sw_diff
from misc import (gather, cpow, remove_duplicate, unravel, ctor2, r2toc,
                  GetSWallSeedsError, n_nearest_indices, find_xs_at_z_0,)

x, z = sympy.symbols('x z')

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
    def __init__(self, z_0=None, x1_0=None, x2_0=None, parents=None,
                 label=None,):
        self.data = [[z_0, x1_0, x2_0]]
        self.parents = parents
        self.label = label

    def get_json_data(self):
        json_data = {
            'data': [[ctor2(z), ctor2(x1), ctor2(x2)] 
                     for z, x1, x2 in self.data],
            'parents': [parent for parent in self.parents],
            'label': self.label,
        }
        return json_data

    def set_json_data(self, json_data):
        self.data = [[r2toc(z), r2toc(x1), r2toc(x2)] 
                      for z, x1, x2 in json_data['data']]
        self.parents = [parent for parent in json_data['parents']]
        self.label = json_data['label']

    def get_zs(self, ti=0, tf=None):
        """
        return a list of (z.real, z.imag)
        """
        if tf is None:
            tf = len(self.data)
        zs = []
        for t in range(ti, tf):
            z = self.data[t][0]
            zs.append(z)
        return zs

    def get_zxzys(self, ti=0, tf=None):
        """
        return a list of (z.real, z.imag)
        """
        if tf is None:
            tf = len(self.data)
        zxzys = []
        for t in range(ti, tf):
            z = self.data[t][0]
            zxzys.append((z.real, z.imag))
        return zxzys

#    def out_of_range(self, z_range_limits):
#        z_f = self.data[-1][0]
#        z_real_min, z_real_max, z_imag_min, z_imag_max = z_range_limits
#        if (z_f.real < z_real_min or
#            z_f.real > z_real_max or
#            z_f.imag < z_imag_min or
#            z_f.imag > z_imag_max):
#            return True
#        else:
#            return False

    def grow(
        self,
        ode,
        ramification_point_zs,
        puncture_point_zs,
        z_range_limits=None,
        num_of_steps=None,
        size_of_small_step=None,
        size_of_large_step=None,
        size_of_neighborhood=None,
    ):
        steps = 0
        rpzs = ramification_point_zs
        ppzs = puncture_point_zs
        z_i, x1_i, x2_i = self.data[-1]
        ode.set_initial_value([z_i, x1_i, x2_i])

        if z_range_limits is not None:
            z_real_min, z_real_max, z_imag_min, z_imag_max = z_range_limits

        while ode.successful() and steps < num_of_steps:
            if (len(rpzs) > 0 and 
                min([abs(z_i - rpz) for rpz in rpzs]) < size_of_neighborhood):
                dt = size_of_small_step
            else:
                dt = size_of_large_step
            z_i, x1_i, x2_i = ode.integrate(ode.t + dt)
            self.data.append([z_i, x1_i, x2_i])

            if (z_i.real < z_real_min or 
                z_i.real > z_real_max or
                z_i.imag < z_imag_min or
                z_i.imag > z_imag_max):
                break

            steps += 1


def get_s_wall_seeds(sw_curve, sw_diff, theta, ramification_point,
                     config_data,):
    rp = ramification_point
    delta = config_data.accuracy
    dt = config_data.size_of_small_step

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
    lambda_0, diff_e = get_local_sw_diff(sw_curve, sw_diff, rp)

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
        logging.error('unknown form of sw_diff at rp ({}, {}): '
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
        xs_at_z_0 = find_xs_at_z_0(sw_curve.num_eq, z_0, rp.x, rp.i)
        dev_phases = [pi for i in range(len(xs_at_z_0)**2)] 
        for i in range(len(xs_at_z_0)):
            diffx = sw_diff.num_v.subs(z, z_0) 
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
            seeds.append((z_0, xs_at_z_0[i], xs_at_z_0[j]))
        
    return seeds


def get_joint(z, x1_i, x2_i, x1_j, x2_j, parent_i, parent_j, accuracy, 
              spectral_network_type='A', label=None):
    """
    return a joint if formed, otherwise return None.
    """
    # TODO: implementation of joint rule changes according to the spectral
    # network type.
    if(abs(x1_i - x2_j) < accuracy and abs(x1_j - x2_i) < accuracy):
        return None
    elif(abs(x2_i - x1_j) < accuracy):
        return Joint(z, x1_i, x2_j, [parent_i.label, parent_j.label], label)
    elif(abs(x2_j - x1_i) < accuracy):
        return Joint(z, x1_j, x2_i, [parent_j.label, parent_i.label], label)


