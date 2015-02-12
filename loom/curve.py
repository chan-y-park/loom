import sympy
import numpy
import logging
import pdb

from math import log10

from misc import ctor2, r2toc, get_root_multiplicity, PSL2C

x, z = sympy.symbols('x z')

class RamificationPoint:
    def __init__(self, z=None, x=None, i=None, label=None, is_puncture=False):
        self.z = z
        self.x = x
        self.i = i
        self.label = label
        self.is_puncture = is_puncture

    def __str__(self):
        return 'z = {}, x = {}, i = {}'.format(self.z, self.x, self.i)

    def __eq__(self, other):
        return self.label == other.label

    def get_json_data(self):
        json_data = {
            'z': ctor2(self.z),
            'x': ctor2(self.x),
            'i': ctor2(self.i),
            'label': self.label,
            'is_puncture': self.is_puncture
        }
        return json_data

    def set_json_data(self, json_data):
        self.z = r2toc(json_data['z'])
        self.x = r2toc(json_data['x'])
        self.i = r2toc(json_data['i'])
        self.label = json_data['label']
        self.is_puncture = json_data['is_puncture']


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
    def __init__(self, config):
        # PSL2C-transformed z & dz
        Cz = PSL2C(config['mt_params'], z) 
        dCz = Cz.diff(z)

        sw_diff_v = sympy.sympify(config['sw_diff_v'])
        v = sw_diff_v.subs(z, Cz) * dCz
        self.sym_v = sympy.simplify(v)
        self.parameters = config['sw_parameters']
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
    def __init__(self, config):
        #PSL2C transformed z
        Cz = PSL2C(config['mt_params'], z) 

        sw_curve = sympy.sympify(config['sw_curve'])
        self.sym_eq = sympy.simplify(sw_curve.subs(z, Cz))
        self.parameters = config['sw_parameters']
        self.num_eq = self.sym_eq.subs(self.parameters)
        logging.info('\nSeiberg-Witten curve: %s = 0\n',
                     sympy.latex(self.num_eq))
        self.accuracy = config['accuracy']

    def get_ramification_points(self, config):
        ramification_points = []
        punctures = [PSL2C(config['mt_params'], p, inverse=True)
                     for p in config['punctures']] 
        f = self.num_eq
        # NOTE: solve_poly_system vs. solve
        #sols = sympy.solve_poly_system([f, f.diff(x)], z, x)
        sols = sympy.solve([f, f.diff(x)], z, x)
        for z_0, x_0 in sols:
            if z_0 in punctures:
                continue
            fx_at_z_0 = f.subs(z, z_0)
            fx_at_z_0_coeffs = map(complex, 
                                  sympy.Poly(fx_at_z_0, x).all_coeffs())
            m = get_root_multiplicity(fx_at_z_0_coeffs, complex(x_0), 
                                      self.accuracy) 
            if m > 1:
                label = 'ramification point #{}'.format(
                    len(ramification_points)
                )
                rp = RamificationPoint(complex(z_0), complex(x_0), m, label)
                ramification_points.append(rp)

        return ramification_points


def get_local_sw_diff(sw_curve, sw_diff, ramification_point):
    rp = ramification_point
    # use Dz = z - rp.z & Dx = x - rp.x
    Dz, Dx = sympy.symbols('Dz, Dx')
    local_curve = (
        sw_curve.num_eq
        .subs(x, rp.x+Dx)
        .subs(z, rp.z+Dz)
        .series(Dx, 0, rp.i+1).removeO()
        .series(Dz, 0, 2).removeO()
    )
    # curve_at_rp = a(z - rp.z) + b(x - rp.x)^(rp.i)
    a = local_curve.n().coeff(Dz).coeff(Dx, 0)
    b = local_curve.n().coeff(Dx**rp.i).coeff(Dz, 0)
    # Dx = Dx(Dz)
    Dx_Dz = (-(a/b)*Dz)**sympy.Rational(1, rp.i)
    local_diff = (
        sw_diff.num_v
        .subs(x, rp.x+Dx_Dz)
        .subs(z, rp.z+Dz)
        .series(Dz, 0, 1).removeO()
    )
    # get the coefficient and the exponent of the leading term
    (diff_c, diff_e) = local_diff.leadterm(Dz)
    if diff_e == 0:
        # remove the constant term from the local_diff
        local_diff -= local_diff.subs(Dz, 0)
        (diff_c, diff_e) = local_diff.leadterm(Dz)

    return (complex(diff_c.n()), diff_e)



