import numpy
#import scipy
import sympy
import logging
#import pdb

from fractions import Fraction
from sympy import limit, oo
#from cmath import exp, log
from cmath import log


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
    return sorted(a_list, cmp=compare)[:n]


def n_nearest_indices(a_list, value, n):
    """
    Find n elements of a_list nearest to value and return their indices,
    by comparing the euclidean norms.
    """
    compare = lambda v1, v2: cmp(abs(v1 - value), abs(v2 - value))
    key = lambda k: a_list[k]
    sorted_indices = sorted(range(len(a_list)), cmp=compare, key=key)
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


def split_with_overlap(np_1d_array, splittings):
    """
    Return a list of splitted numpy 1d array with an overlapping element.
    """
    a = np_1d_array
    ss = splittings
    n_segs = len(ss) + 1

    if n_segs == 1:
        return [a]

    # start with initial piece
    segs = [a[:ss[0]+1]]
    for i in range(len(ss)-1):
        s_0 = ss[i]
        s_1 = ss[i+1] + 1
        if s_0 < 0 or s_1 < 0:
            raise ValueError("split_with_overlap(): overlap is too large.")
        segs.append(a[s_0:s_1])
    # add the last piece
    s_f = ss[-1]
    if s_f < 0:
        raise ValueError("split_with_overlap(): overlap is too large.")
    segs.append(a[s_f:])

    return segs
 
