import numpy
import logging

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

