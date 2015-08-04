import sympy
import numpy
import logging
import pdb
import subprocess

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

class SWCurve:
    def __init__(self, differentials=None, g_data=None):
        N = len(g_data['weights'])
        self.eq_str = 'x^{} '.format(N)
        for k, u_k in differentials.iteritems():
            self.eq_str += '+ ({}) '.format(u_k)
            if k != N:
                self.eq_str += '* x^{}'.format(N-k)
        self.sym_eq = None
        self.num_eq = None


class SWDiff:
    def __init__(self, v_str):
        self.v_str = v_str
        self.sym_v = None
        self.num_v = None


class SWData:
    """
    A class containing a Seiberg-Witten curve 
        \lambda^N + \sum_{k=2}^N \phi_k \lambda^{N-k}, 
    where \lambda is the Seiberg-Witten differential of the form 
        \lambda = x dz
    and \phi_k = u_k(z) dz^k.
    """
    def __init__(self, config):
        self.parameters = config['sw_parameters']
        self.differentials = eval(config['differentials'])

        self.g_data = sage_get_g_data(config)

        self.curve = SWCurve(
            differentials=self.differentials, 
            g_data=self.g_data,
        )
        self.diff = SWDiff('x')
        self.punctures = None

        # PSL2C-transformed z & dz
        Cz = PSL2C(config['mt_params'], z, inverse=True) 
        dCz = Cz.diff(z)

        # Seiberg-Witten curve
        self.curve.sym_eq = sympy.simplify(
            sympy.sympify(self.curve.eq_str).subs(z, Cz)
        )
        # num_eq is from sym_eq with its parameters 
        # substituted with numerical values.
        self.curve.num_eq = self.curve.sym_eq.subs(self.parameters)
        logging.info('\nSeiberg-Witten curve: %s = 0\n',
                     sympy.latex(self.curve.num_eq))

        # Seiberg-Witten differential
        # sym_v is a SymPy expression. 
        self.diff.sym_v = sympy.simplify(
            sympy.sympify(self.diff.v_str).subs(z, Cz) * dCz
        )
        # num_v is from sym_v with its parameters 
        # substituted with numerical values.
        self.diff.num_v = self.diff.sym_v.subs(self.parameters)
        logging.info('\nSeiberg-Witten differential: %s dz\n',
                     sympy.latex(self.diff.num_v))

        if config['punctures'] is None:
            self.punctures = []
        else:
            self.punctures = [
                PSL2C(config['mt_params'], p, numerical=True)
                for p in config['punctures']
            ]


def get_ramification_points(sw, accuracy):
    f = sw.curve.num_eq

    ramification_points = []

    # NOTE: solve_poly_system vs. solve
    #sols = sympy.solve_poly_system([f, f.diff(x)], z, x)
    #sols = sympy.solve([f, f.diff(x)], z, x)
    #if sols is None:
    #    # Use Sage instead
    #    sols = sage_solve_poly_system([f, f.diff(x)])
    sols = sage_solve_poly_system([f, f.diff(x)])
    for z_0, x_0 in sols:
        if (len(sw.punctures) > 0 and
            (min([abs(z_0 - p) for p in sw.punctures]) < accuracy)
        ):
            continue
        fx_at_z_0 = f.subs(z, z_0)
        fx_at_z_0_coeffs = map(complex, sympy.Poly(fx_at_z_0, x).all_coeffs())
        mx = get_root_multiplicity(fx_at_z_0_coeffs, complex(x_0), accuracy)
        if mx > 1:
            fz_at_x_0 = f.subs(x, x_0)
            fz_at_x_0_coeffs = map(complex, 
                                   sympy.Poly(fz_at_x_0, z).all_coeffs())
            
            mz = get_root_multiplicity(fz_at_x_0_coeffs, complex(z_0),
                                       accuracy)
            if mz > 1:
                continue
            label = 'ramification point #{}'.format(len(ramification_points))
            rp = RamificationPoint(complex(z_0), complex(x_0), mx, label)
            logging.info("{}: z = {}, x = {}, i = {}."
                         .format(label, rp.z, rp.x, rp.i))
            ramification_points.append(rp)
    return ramification_points


def get_local_sw_diff(sw, ramification_point):
    rp = ramification_point
    # use Dz = z - rp.z & Dx = x - rp.x
    Dz, Dx = sympy.symbols('Dz, Dx')
    local_curve = (
        sw.curve.num_eq
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
        sw.diff.num_v
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


def get_fibers(config, z_0):
    """
    Return a list of x-coordinates of the fibers over z.
    """
    sw = SWData(config)
    fx = sw.curve.num_eq.subs(z, z_0)
    xs = sympy.solve(fx, x)
    return map(complex, xs) 


def sage_solve_poly_system(poly_system):
    """
    Use sage to solve the given system of polynomial equations of x and z.
    """
    sols_str = subprocess.check_output(
        ["sage", "./loom/sage_scripts/solve_poly_system.sage"] +
        [str(poly) for poly in poly_system]
    )
    sols = eval(sols_str)
    return sols

def sage_get_g_data(config):
    root_system = config['root_system']
    g_data_str = subprocess.check_output(
        ["sage", "./loom/sage_scripts/get_g_data.sage", 
         root_system, config["representation"]]
    )
    g_data = eval(g_data_str)

    g_data['root_system'] = root_system
    g_data['representation'] = eval(config['representation'])
    g_data['type'] = root_system[0]
    g_data['rank'] = eval(root_system[1:])

    return g_data
