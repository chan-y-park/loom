import numpy
import sympy
import logging

from fractions import Fraction

def get_root_multiplicity(coefficients, root_0, accuracy):
    """
    Returns the multiplicity of a root of a polynomial
    defined by the list 'coefficients' according to numpy.roots
    """
    multiplicity = 0
    #logging.debug('accuracy = %f', accuracy)
    roots = numpy.roots(coefficients)
    #logging.debug('numpy.roots -> %s', roots)
    for root_i in roots:
        if abs(root_0 - root_i) < accuracy:
            multiplicity += 1

    #logging.debug('multiplicity = %d', multiplicity)
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
    Find n elements of a_list nearest to value and return their indices.
    """
    
    compare = lambda v1, v2: cmp(abs(v1 - value), abs(v2 - value))
    key = lambda k: a_list[k]
    sorted_indices = sorted(range(len(a_list)), compare, key)

    min_indices = sorted_indices[:n]

    size = sympy.sqrt(len(a_list))
    if not size.is_Integer:
        logging.error('a_list is not a square matrix.')
        raise NNearestError(size)
    else:
        size = int(size)
    result = []
    for k in min_indices:
        i = k / size
        j = k % size
        result.append((i, j))

    return result

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

#from random import randint

#l = [randint(1, 9) for i in range(9)]
#print l
#print n_nearest(l, 2, 3)
#print gather(l, lambda x, y: x == y)
#print remove_duplicate(l, lambda x, y: x == y)

