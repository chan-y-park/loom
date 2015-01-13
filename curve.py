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
        self.sym_equation = sympy.sympify(config_data.curve_eq_string)
        self.parameters = config_data.curve_parameters 
        self.ramification_points = []
        self.puncture_points = []

    def find_ramification_points(self, accuracy):
        f = self.sym_equation.subs(self.parameters)
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
            m = get_root_multiplicity(fx_at_z0_coeffs, complex(x0), accuracy) 
            if m > 1:
                rp = RamificationPoint(complex(z0), complex(x0), m)
                logging.debug('rp = %s', rp)
                self.ramification_points.append(rp)

