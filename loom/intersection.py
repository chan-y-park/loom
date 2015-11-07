"""
A general intersection module.

Objects and functions to find intersections of real 1-dim curves
on a real 2-dim plane.
"""
#import logging
import numpy
from scipy.interpolate import interp1d
from scipy.optimize import brentq, newton
from sympy import Interval
from itertools import combinations

import pdb

class NoIntersection(Exception):
    """
    An exception class to be raised when failed
    to find an intersection.
    """
    def __init__(self, value=''):
        self.value = value

    def __str__(self):
        return repr(self.value)


def remove_duplicate_intersection(new_ilist, old_ilist):
    """
    Remove any duplicate in new_ilist1, then remove intersections
    of new_ilist that already exist in old_ilist
    """
    temp_ilist = new_ilist

    for intersection1, intersection2 in combinations(temp_ilist, 2):
        if intersection1 == intersection2:
            new_ilist.remove(intersection2)
        else:
            continue

    temp_ilist = new_ilist

    for new_intersection in temp_ilist:
        for intersection in old_ilist:
            if new_intersection == intersection:
                new_ilist.remove(new_intersection)


def find_curve_range_intersection(curve_1, curve_2, cut_at_inflection=False):
    """
    Return intersections of x- and y-ranges of two real curves,
    which are parametric curves on the xy-plane given as
    (x_array, y_array), a tuple of NumPy arrays.
    """
    x1, y1 = curve_1
    x2, y2 = curve_2

    if cut_at_inflection is True:
        x1_min, x1_max = sorted([x1[0], x1[-1]])
        x2_min, x2_max = sorted([x2[0], x2[-1]])

        y1_min, y1_may = sorted([y1[0], y1[-1]])
        y2_min, y2_may = sorted([y2[0], y2[-1]])
    else:
        x1_min, x1_max = numpy.sort(x1)[[0, -1]]
        x2_min, x2_max = numpy.sort(x2)[[0, -1]]

        y1_min, y1_may = numpy.sort(y1)[[0, -1]]
        y2_min, y2_may = numpy.sort(y2)[[0, -1]]

    x1_interval = Interval(x1_min, x1_max)
    x2_interval = Interval(x2_min, x2_max)

    y1_interval = Interval(y1_min, y1_may)
    y2_interval = Interval(y2_min, y2_may)

    x_range = x1_interval.intersect(x2_interval)
    y_range = y1_interval.intersect(y2_interval)

    return (x_range, y_range)


def find_intersection_of_segments(segment_1, segment_2, accuracy=1e-1,
                                  newton_maxiter=5):
    """
    Find an intersection of two segments of curves.

    First find interpolations of segments using scipy.interp1d and
    use SciPy's Brent method to find an intersection. When this
    doesn't work, use SciPy's polynomial interpolation and then use
    secant method to find an intersection.

    Parameters
    ----------
    segment_1, segment_2: Segments to find their intersection. Each
        segment is (x_array, y_array), a tuple of NumPy arrays.
    newton_maxiter: Maximum number of iterations for secant method.
        When increased, this gives a better accuracy of the
        intersection but it also greatly reduces the performance
        due to many cases when there is no intersection but
        the module keeps trying to find one.
    """
    # First check if the two segments share any x- and y-range.
    x_range, y_range = find_curve_range_intersection(
        segment_1, segment_2, cut_at_inflection=True
    )
    if (x_range.is_EmptySet or y_range.is_EmptySet or x_range.is_FiniteSet or
            y_range.is_FiniteSet):
        # The segments and the bin do not share a domain and therefore
        # there is no intersection.
        raise NoIntersection()

    f1 = interp1d(*segment_1)
    f2 = interp1d(*segment_2)
    delta_f12 = lambda x: f1(x) - f2(x)

    try:
        #logging.debug('try brentq.')
        intersection_x = brentq(delta_f12, x_range.start, x_range.end)
        intersection_y = f1(intersection_x)

    except ValueError:

        x0 = 0.5*(x_range.start + x_range.end)

        try:
            #logging.debug('try newton with x0 = %.8f.', x0)
            intersection_x = newton(delta_f12, x0, maxiter=newton_maxiter)
            #logging.debug('intersection_x = %.8f.', intersection_x)
        except ValueError:
            # newton() searches for x outside the interpolation domain.
            # Declare no intersection.
            raise NoIntersection()
        except RuntimeError:
            # Newton's method fails to converge; declare no intersection
            raise NoIntersection()

        # Verify the solution returned by newton().
        if abs(delta_f12(intersection_x)) > accuracy:
            raise NoIntersection()

        # Check if the intersection is within the curve range.
        # If not, the intersecion is not valid.
        if intersection_x not in x_range:
            raise NoIntersection()
        intersection_y = f1(intersection_x)
        if intersection_y not in y_range:
            raise NoIntersection()

    return [float(intersection_x), float(intersection_y)]
