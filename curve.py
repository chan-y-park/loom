import pdb
import sympy
import numpy
import logging

from misc import get_root_multiplicity

x, z = sympy.symbols('x z')

class RamificationPoint:
    def __init__(self, z, x, i, label=None, is_puncture=False):
        self.z = z
        self.x = x
        self.i = i
        self.label = label
        self.is_puncture = is_puncture

    def __str__(self):
        return 'z = {}, x = {}, i = {}'.format(self.z, self.x, self.i)

    def __eq__(self, other):
        return self.label == other.label


class PuncturePoint:
    def __init__(self, z, cutoff, label):
        self.z = z
        self.cutoff = cutoff
        self.label = label

    def __eq__(self, other):
        return self.label == other.label

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
        #logging.debug('sympy.solve: %s\n', sols)
        for z_0, x_0 in sols:
            #logging.debug('sympy: z_0 = %s, x_0 = %s', z_0, x_0)
            fx_at_z_0 = f.subs(z, z_0)
            fx_at_z_0_coeffs = map(complex, 
                                  sympy.Poly(fx_at_z_0, x).all_coeffs())
            m = get_root_multiplicity(fx_at_z_0_coeffs, complex(x_0), 
                                      self.accuracy) 
            if m > 1:
                label = 'ramification point #{}'.format(
                    len(self.ramification_points)
                )
                rp = RamificationPoint(complex(z_0), complex(x_0), m, label)
                self.ramification_points.append(rp)

