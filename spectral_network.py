import numpy
import sympy
import logging
import pdb

from itertools import combinations
from cmath import exp, pi, phase
from fractions import Fraction

from curve import RamificationPoint, SWCurve, SWDiff
from s_wall import SWall
from misc import (cpow, gather, remove_duplicate, n_nearest, LocalDiffError,)

x, z = sympy.symbols('x z')

class SpectralNetwork:
    def __init__(self, sw_curve, sw_diff, theta, config_data):
        self.sw_curve = sw_curve
        self.sw_diff = sw_diff
        self.theta = theta
        self.config_data = config_data
        self.set_ode_f()
        self.generate_from_ramification_points()

    def set_ode_f(self):
        logging.debug('need to implement')

# NOTE: tried SymPy exact numbers but didn't work.
#    def get_s_wall_seeds(self, ramification_point):
#        rp = ramification_point
#        theta = self.theta
#
#        omega_1 = sympy.exp(sympy.sympify('2*pi*I/{}'.format(rp.i)))
#        omega = [omega_1**k for k in range(rp.i)]
#        logging.debug('omega = %s', omega)
#
#        beta_1 = sympy.sympify('exp({}*2*pi*I/({}+1))'.format(rp.i, rp.i))
#        beta = [beta_1**k for k in range(rp.i+1)]
#        logging.debug('beta = %s', beta)
#
#        cs = []
#
#        for i in range(1, rp.i):
#            new_locs = []
#            for j in range(rp.i): 
#                if j+i < rp.i:
#                    new_loc = 1/(omega[j]-omega[j+i])
#                else:
#                    new_loc = 1/(omega[j]-omega[j+i-rp.i])
#                new_loc = (new_loc**sympy.Rational(rp.i, rp.i+1)*
#                           sympy.exp(theta*sympy.I*sympy.Rational(rp.i, 
#                                                                  rp.i+1))
#                          )
#                new_locs += [new_loc*beta_i for beta_i in beta]
#            # End of j loop
#            logging.debug('len(new_locs) = %d', len(new_locs))
#            pdb.set_trace()
#            new_locs = remove_duplicate(new_locs, 
#                             lambda l1, l2: sympy.simplify(l1 - l2) == 0)
#            pdb.set_trace()
#            logging.debug('len(new_locs) = %d', len(new_locs))
#            cs += new_locs
#
#        logging.debug('cs = %s, # = %d', cs, len(cs))
#
#        for c in cs:
#            logging.debug('c = %s', complex(c))


    def get_s_wall_seeds(self, ramification_point):
        rp = ramification_point
        theta = self.theta
        delta = self.config_data.accuracy
        dz = self.config_data.size_of_small_step

        ###
        # 1. Find the first-order approximations of the starting points 
        # of S-walls around a given ramification point, which is of the form
        # \Delta z_i = c_i / (\lambda_0)^(rp.i/(rp.i+1))
        ###

        # 1.1 find c_i, a phase factor for each S-wall.
        omega_1 = exp(2*pi*1j/rp.i)
        omega = [omega_1**k for k in range(rp.i)]
        #logging.debug('omega = %s', omega)

        beta_1 = exp(rp.i*2*pi*1j/(rp.i+1))
        beta = [beta_1**k for k in range(rp.i+1)]
        #logging.debug('beta = %s', beta)

        cs = []

        # go over pairs of omegas that differ by \omega_1^i
        for i in range(1, rp.i):
            new_locs = []
            # and go over all the omegas
            for j in range(rp.i): 
                if j+i < rp.i:
                    new_loc = 1/(omega[j]-omega[j+i])
                else:
                    new_loc = 1/(omega[j]-omega[j+i-rp.i])
                new_loc = (cpow(new_loc, rp.i, rp.i+1)*
                           exp(theta*1j*Fraction(rp.i, rp.i+1)))
                new_locs += [new_loc*beta_i for beta_i in beta]
            # End of j loop
            #logging.debug('len(new_locs) = %d', len(new_locs))
            new_locs = remove_duplicate(new_locs, 
                             lambda l1, l2: abs(l1 - l2) < delta)
            #logging.debug('len(new_locs) = %d', len(new_locs))
            cs += new_locs

        cs = gather(cs, lambda c1, c2: abs(c1 - c2) < delta)
        logging.debug('list of c = %s, # = %d', cs, len(cs))

        # 1.2 find the coefficient and the exponent of the leading term
        # of the SW differential at the ramification point.
        lambda_0, diff_e = get_local_sw_diff(self.sw_curve,
                                                  self.sw_diff, rp)
        if(diff_e != sympy.Rational(1, rp.i)):
            logging.error('unknown form of sw_diff at rp ({}, {}): '
                          '{} != 1/rp.i = 1/{}'.format(rp.z, rp.x, diff_e, 
                                                       rp,i)
                         )
            raise LocalDiffError(diff_e)

        # 2. Now calculate \Delta z_i for each S-wall and 
        # find the two points on the curve that are projected onto it.  
        for c in cs:
            Delta_z = c[0] / cpow(lambda_0, rp.i, rp.i+1)
            # resize to the size of the small step 
            Delta_z = Delta_z/abs(Delta_z)*dz
            z_0 = rp.z + Delta_z
            fx_at_z_0 = self.sw_curve.num_eq.subs(z, z_0)
            fx_at_z_0_coeffs = map(complex, 
                                  sympy.Poly(fx_at_z_0, x).all_coeffs())
            xs_at_z_0 = sorted(numpy.roots(fx_at_z_0_coeffs),
                               lambda x1, x2: cmp(abs(x1 - rp.x),
                                                  abs(x2 - rp.x)))[:rp.i]
            dev_phases = [pi for i in range(len(xs_at_z_0)**2)] 
            for i in range(len(xs_at_z_0)):
                diffx = self.sw_diff.num_v.subs(z, z_0) 
                v_i = diffx.subs(x, xs_at_z_0[i])
                for j in range(len(xs_at_z_0)):
                    if i == j:
                        continue
                    else:
                        v_j = diffx.subs(x, xs_at_z_0[j]) 
                        delta_z = exp(self.theta)/(v_i - v_j)
                        # flattened index
                        fij = i*len(xs_at_z_0) + j
                        dev_phases[fij] = phase((delta_z/Delta_z))
            min_dev_indices = n_nearest(dev_phases, 0, len(c))
            pdb.set_trace()

    def generate_from_ramification_points(self):
        for rp in self.sw_curve.ramification_points:
            self.get_s_wall_seeds(rp)

    def generate_from_punctures(self, **params):
        logging.debug('need to implement')

    def find_joints(self):
        logging.debug('need to implement')

    def generate_from_joints(self, **params):
        logging.debug('need to implement')

    def grow(self, **params):
        logging.debug('need to implement')

    def plot(self):
        logging.debug('need to implement')

def get_local_sw_diff(sw_curve, sw_diff, ramification_point):
    rp = ramification_point
    local_curve = sw_curve.num_eq.series(x, rp.x, rp.i+1)
    local_curve = local_curve.series(z, rp.z, 2).removeO()
    #logging.debug('local curve at rp = ({}, {}): '
    #              '{}'.format(rp.z, rp.x, local_curve))
    # curve_at_rp = a(z - rp.z) + b(x - rp.x)^(rp.i)
    a = local_curve.coeff(z).coeff(x, 0)
    b = local_curve.coeff(x**rp.i).coeff(z, 0)
    #logging.debug('a = {}, b = {}'.format(a, b))
    local_x = rp.x + (-(a/b)*(z - rp.z))**sympy.Rational(1, rp.i)
    # substitute x with x(z)
    local_diff = sw_diff.num_v.subs(x, local_x)
    # series expansion in z at rp.z
    local_diff = local_diff.series(z, rp.z, 1)
    if local_diff.getO() is None:
        # series expansion didn't occur.
        # translate z such that rp.z = 0
        local_diff = local_diff.subs(z, z+rp.z)
    # get the coefficient and the exponent of the leading term
    (diff_c, diff_e) = local_diff.leadterm(z)
    if diff_e == 0:
        # remove the constant term from the local_diff
        local_diff -= local_diff.subs(z, 0)
        (diff_c, diff_e) = local_diff.leadterm(z)

    return (complex(diff_c), diff_e)

#def get_min_devs(dev_phases, num):
#    """
#    From the list of deviation phases, pick num phases that are closest to 0
#    and return there unflattened indices.
#    """
#    
#    sorted_indices = sorted(range(len(dev_phases)),
#                            lambda k1, k2: cmp(abs(dev_phases[k1]),
#                                               abs(dev_phases[k2])),
#                            key=lambda k: s[k])
#
#    min_indices = sorted_indices[:num]
#
#    size = sympy.sqrt(len(dev_phases))
#    if not size.is_Integer:
#        logging.error('dev_phases is not a square matrix.')
#        raise GetMinDevsError(size)
#    result = []
#    for k in min_indices:
#        i = k % size
#        j = k - i*size
#        result.append((i, j))
#
#    return result
#

def generate_spectral_network(config_data):
    sw_curve = SWCurve(config_data)
    sw_curve.find_ramification_points()
    sw_diff = SWDiff(config_data)

    logging.info('\nList of ramification points')
    for rp in sw_curve.ramification_points:
        logging.info('%s', rp)
    logging.info('\n')

    if(config_data.single_network == True):
        spectral_network = SpectralNetwork(sw_curve, sw_diff, 
                                           config_data.phase, config_data) 


