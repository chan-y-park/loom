import numpy
import scipy
import sympy
import logging
import pdb

from fractions import Fraction
from sympy import limit, oo
from cmath import exp, log


class LocalDiffError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class GetSWallSeedsError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class NNearestError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class UnravelError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def ctor2(complex_number):
    return (complex_number.real, complex_number.imag)


def r2toc(real_tuple):
    return (real_tuple[0]+1j*real_tuple[1])


def get_root_multiplicity(coefficients, root_0, accuracy):
    """
    Returns the multiplicity of a root of a polynomial
    defined by the list 'coefficients' according to numpy.roots
    """
    multiplicity = 0
    roots = numpy.roots(coefficients)
    for root_i in roots:
        if abs(root_0 - root_i) < accuracy:
            multiplicity += 1

    return multiplicity


def cpow(base, exponent_numerator, exponent_denominator=1):
    return complex(base)**complex(Fraction(exponent_numerator,
                                           exponent_denominator))


def gather(a_list, compare, result=None):
    if result is None:
        result = []
    e_0 = a_list[0]
    result.append([e_0, 0])
    next_list = []
    pop_list = []
    for i in range(len(a_list)):
        if compare(e_0, a_list[i]) is True:
            pop_list.append(i)
    for i, e_i in enumerate(a_list):
        if i in pop_list:
            result[-1][1] += 1
        else:
            next_list.append(e_i)
    if len(next_list) == 0:
        return result
    else:
        return gather(next_list, compare, result)


def remove_duplicate(a_list, compare):
    return [e[0] for e in gather(a_list, compare)]


def n_nearest(a_list, value, n):
    """
    Find n elements of a_list nearest to value and return them,
    by comparing the euclidean norms.
    """
    compare = lambda v1, v2: cmp(abs(v1 - value), abs(v2 - value))
    return sorted(a_list, compare)[:n]


def n_nearest_indices(a_list, value, n):
    """
    Find n elements of a_list nearest to value and return their indices,
    by comparing the euclidean norms.
    """
    compare = lambda v1, v2: cmp(abs(v1 - value), abs(v2 - value))
    key = lambda k: a_list[k]
    sorted_indices = sorted(range(len(a_list)), compare, key)

    return sorted_indices[:n]


def unravel(k, row_size, column_size=None):
    """
    return unflattened indices of a matrix.
    """
    if column_size is None:
        column_size = row_size

    if k < 0:
        logging.error('k = {} should not be negative.'.format(k))
        raise UnravelError(k)

    i = k / row_size
    j = k % row_size

    if i > column_size:
        logging.error('i = {} should not be greater than '
                      'the size of the column = {}.'.format(i, column_size))
        raise UnravelError(k)

    return (i, j)

#def find_xs_at_z_0(f_z_x, z_0, x_0=None, num_x=1):
#    """
#    solve f(x, z_0) = 0 and return num_x x's nearest to x_0.
#    """
#    x, z = sympy.symbols('x z')
#
#    f_x_at_z_0 = f_z_x.subs(z, z_0)
#    f_x_at_z_0_coeffs = map(complex, sympy.Poly(f_x_at_z_0, x).all_coeffs())
#    xs_at_z_0 = numpy.roots(f_x_at_z_0_coeffs)
#    if x_0 is None:
#        return xs_at_z_0
#    else:
#        return sorted(xs_at_z_0,
#                      lambda x1, x2: cmp(abs(x1 - x_0), abs(x2 - x_0)))[:num_x]

def PSL2C(C, z, inverse=False, numerical=False):
    """
    Apply linear fractional transformation defined by C onto z.
    """
    if C is None:
        C = [[1, 0], [0, 1]]

    if inverse is True:
        a = C[1][1]
        b = -C[0][1]
        c = -C[1][0]
        d = C[0][0]
    else:
        a = C[0][0]
        b = C[0][1]
        c = C[1][0]
        d = C[1][1]

    if numerical is True:
        return (complex(a) * z + complex(b))/(complex(c) * z + complex(d))
    else:
        u = sympy.symbols('u')
        Cu = (a*u+b)/(c*u+d)

        if z == oo:
            Cz = limit(Cu, u, z)
        else:
            Cz = Cu.subs(u, z)

        return Cz


def put_on_cylinder(z, mt_params=None):
    """
    Put PSL2C-transformed z-coords
    onto the original cylinder.
    """
    return log(PSL2C(mt_params, z, inverse=True, numerical=True))/1.0j


def get_ode(sw, phase, accuracy):
    x, z = sympy.symbols('x z')
    ode_absolute_tolerance = accuracy

    f = sw.curve.num_eq
    df_dz = f.diff(z)
    df_dx = f.diff(x)
    # F = -(\partial f/\partial z)/(\partial f/\partial x)
    F = sympy.lambdify((z, x), -df_dz/df_dx)
    v = sympy.lambdify((z, x), sw.diff.num_v)

    def ode_f(t, zx1x2M):
        z_i = zx1x2M[0]
        x1_i = zx1x2M[1]
        x2_i = zx1x2M[2]
        dz_i_dt = exp(phase*1j)/(v(z_i, x1_i) - v(z_i, x2_i))
        dx1_i_dt = F(z_i, x1_i) * dz_i_dt
        dx2_i_dt = F(z_i, x2_i) * dz_i_dt
        dM_dt = 1
        return [dz_i_dt, dx1_i_dt, dx2_i_dt, dM_dt]

    ode = scipy.integrate.ode(ode_f)
    ode.set_integrator(
        'zvode',
        #method='adams',
        atol=ode_absolute_tolerance,
    )

    return ode


### chan: TODO: use numba?
def delete_duplicates(l):
    seen = set()
    uniq = []
    for x in l:
        if x not in seen:
            uniq.append(x)
            seen.add(x)
    return uniq


def clock(direction):
    if direction == 'left':
        return 'ccw'
    elif direction == 'right':
        return 'cw'
    else:
        logging.info('\nCannot read direction!\n')


def left_right(l, point):
    """
    given the list 
    l = [..., z, ...]
    and a point in the list (specified by the corresponding integer),
    determines whether x increases or decreases at that point, 
    returning repsectively 'left' or 'right'
    """
    if point > len(l)-1:
        logging.info('Cant determine direction, point doesnt belong to list!')
    elif point > 0:
        if l[point-1].real < l[point].real:
            return 'right'
        else:
            return 'left'
    elif point == 0:
        if l[point].real < l[point+1].real:
            return 'right'
        else:
            return 'left'


def is_root(np_array, g_data):
    ans = False
    for rt in list(g_data.roots):
        if (np_array == rt).all():
            ans = True
            break
        else:
            pass
    return ans