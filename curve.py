import pdb
import sympy
import numpy
import logging

from misc import get_root_multiplicity

x, z = sympy.symbols('x z')

class RamificationPoint:
    def __init__(self, z, x, i):
        self.z = z
        self.x = x
        self.i = i

    def __str__(self):
        return 'z = {}, x = {}, i = {}'.format(self.z, self.x, self.i)

class PuncturePoint:
    def __init__(self, z, cutoff):
        self.z = z
        self.cutoff = cutoff

class SWDiff:
    """
        Define a Seiberg-Witten differential of the form 
            \lambda = v(x, z) dz
    """
    def __init__(self, config_data):
        self.sym_v = sympy.sympify(config_data.sw_diff_v_string)
        self.parameters = config_data.sw_parameters 
        self.num_v = self.sym_v.subs(self.parameters)
        logging.info('\nSeiberg-Witten differential: %s dz\n',
                     sympy.latex(self.num_v))

class SWCurve:
    """
    a class containing a Seiberg-Witten curve and relevant information.

    Attributes:
        eq_string: equation for the SW curve in string
        equation: equation as a SymPy function
        ramification_points: list of RamificationPoint instances
        puncture_points: list of PuncturePoint instances
    Methods:
        set_curve_parameters
        find_ramification_points
    """
    def __init__(self, config_data):
        self.sym_eq = sympy.sympify(config_data.sw_curve_eq_string)
        self.parameters = config_data.sw_parameters 
        self.num_eq = self.sym_eq.subs(self.parameters)
        logging.info('\nSeiberg-Witten curve: %s = 0\n',
                     sympy.latex(self.num_eq))
        self.accuracy = config_data.accuracy
        self.ramification_points = []
        self.puncture_points = []

    def find_ramification_points(self):
        f = self.num_eq
        # NOTE: solve_poly_system vs. solve
        #sols = sympy.solve_poly_system([f, f.diff(x)], z, x)
        #logging.debug('sympy.solve_poly_system: %s\n', sols)
        sols = sympy.solve([f, f.diff(x)], z, x)
        logging.debug('sympy.solve: %s\n', sols)
        for z0, x0 in sols:
            logging.debug('sympy: z0 = %s, x0 = %s', z0, x0)
            fx_at_z0 = f.subs(z, z0)
            fx_at_z0_coeffs = map(complex, 
                                  sympy.Poly(fx_at_z0, x).all_coeffs())
            m = get_root_multiplicity(fx_at_z0_coeffs, complex(x0), 
                                      self.accuracy) 
            if m > 1:
                rp = RamificationPoint(complex(z0), complex(x0), m)
                logging.debug('rp = %s', rp)
                self.ramification_points.append(rp)

