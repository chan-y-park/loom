import numpy
import sympy
import logging
import pdb

from fractions import Fraction
from sympy import limit, oo


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

#def remove_duplicate(a_list, compare, i0=0):
#    for i in range(i0, len(a_list)):
#        for j in range(i+1, len(a_list)):
#            if compare(a_list[i], a_list[j]) is True:
#                a_list.pop(j)
#                return remove_duplicate(a_list, compare, i)
#
#    return a_list

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

def find_xs_at_z_0(f_z_x, z_0, x_0=None, num_x=1):
    """
    solve f(x, z_0) = 0 and return num_x x's nearest to x_0.
    """
    x, z = sympy.symbols('x z')

    f_x_at_z_0 = f_z_x.subs(z, z_0)
    f_x_at_z_0_coeffs = map(complex, sympy.Poly(f_x_at_z_0, x).all_coeffs())
    xs_at_z_0 = numpy.roots(f_x_at_z_0_coeffs)
    if x_0 is None:
        return xs_at_z_0
    else:
        return sorted(xs_at_z_0,
                      lambda x1, x2: cmp(abs(x1 - x_0), abs(x2 - x_0)))[:num_x]

def PSL2C(C, z, inverse=False):
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

    
    u = sympy.symbols('u')
    Cu = (a*u+b)/(c*u+d)

    if z == oo:
        Cz = limit(Cu, u, z) 
    else:
        Cz = Cu.subs(u, z)
    return Cz 

#def n_nearest(a_list, value, n):
#    """
#    Find n elements of a_list nearest to value and return their indices.
#    """
#    
#    compare = lambda v1, v2: cmp(abs(v1 - value), abs(v2 - value))
#    key = lambda k: a_list[k]
#    sorted_indices = sorted(range(len(a_list)), compare, key)
#
#    min_indices = sorted_indices[:n]
#
#    size = sympy.sqrt(len(a_list))
#    if not size.is_Integer:
#        logging.error('a_list is not a square matrix.')
#        raise NNearestError(size)
#    else:
#        size = int(size)
#    result = []
#    for k in min_indices:
#        i = k / size
#        j = k % size
#        result.append((i, j))
#
#    return result

#from random import randint

#l = [randint(1, 9) for i in range(9)]
#print l
#print n_nearest(l, 2, 3)
#print gather(l, lambda x, y: x == y)
#print remove_duplicate(l, lambda x, y: x == y)

