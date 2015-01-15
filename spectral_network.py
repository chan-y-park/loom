import sympy
import logging
import pdb

from itertools import combinations
from cmath import exp, pi
from fractions import Fraction

from curve import RamificationPoint, SWCurve, SWDiff
from s_wall import SWall
from misc import cpow, remove_duplicate

class SpectralNetwork:
    def __init__(self, sw_curve, sw_diff, theta, config_data):
        self.sw_curve = sw_curve
        self.sw_diff = sw_diff
        self.theta = theta
        self.set_ode_f()
        self.accuracy = config_data.accuracy
        self.generate_from_ramification_points()

    def set_ode_f(self):
        logging.debug('need to implement')

    def generate_from_ramification_points(self):
        #for rp in self.sw_curve.ramification_points:
        #    self.get_s_wall_seeds(rp)
        # for debugging
        rp = RamificationPoint(0, 0, 3)
        self.get_s_wall_seeds(rp)

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

        omega_1 = exp(2*pi*1j/rp.i)
        omega = [omega_1**k for k in range(rp.i)]
        logging.debug('omega = %s', omega)

        beta_1 = exp(rp.i*2*pi*1j/(rp.i+1))
        beta = [beta_1**k for k in range(rp.i+1)]
        logging.debug('beta = %s', beta)

        cs = []

        for i in range(1, rp.i):
            new_locs = []
            for j in range(rp.i): 
                if j+i < rp.i:
                    new_loc = 1/(omega[j]-omega[j+i])
                else:
                    new_loc = 1/(omega[j]-omega[j+i-rp.i])
                new_loc = (cpow(new_loc, rp.i, rp.i+1)*
                           exp(theta*1j*Fraction(rp.i, rp.i+1))
                          )
                new_locs += [new_loc*beta_i for beta_i in beta]
            # End of j loop
            logging.debug('len(new_locs) = %d', len(new_locs))
            #pdb.set_trace()
            new_locs = remove_duplicate(new_locs, 
                             lambda l1, l2: abs(l1 - l2) < self.accuracy)
            #pdb.set_trace()
            logging.debug('len(new_locs) = %d', len(new_locs))
            cs += new_locs

        logging.debug('cs = %s, # = %d', cs, len(cs))

        for c in cs:
            print('{{{}, {}}},'.format(c.real, c.imag))

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


