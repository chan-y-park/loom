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
    logging.debug('accuracy = %f', accuracy)
    roots = numpy.roots(coefficients)
    logging.debug('numpy.roots -> %s', roots)
    for root_i in roots:
        if abs(root_0 - root_i) < accuracy:
            multiplicity += 1

    logging.debug('multiplicity = %d', multiplicity)
    return multiplicity

def cpow(base, exponent_numerator, exponent_denominator=1):
    return complex(base)**complex(Fraction(exponent_numerator,
                                           exponent_denominator))

def remove_duplicate(a_list, compare, i0=0):
    for i in range(i0, len(a_list)):
        for j in range(i+1, len(a_list)):
            if compare(a_list[i], a_list[j]) is True:
                a_list.pop(j)
                return remove_duplicate(a_list, compare, i)

    return a_list

#from random import randint

#l = [randint(1, 20) for i in range(20)]
#print remove_duplicate(l, lambda x, y: x == y)
