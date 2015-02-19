"""
A general intersection module.

Objects and functions to find intersections of real 1-dim curves
on a real 2-dim plane.

This module is based on the following libraries:
NumPy 1.8.2
SciPy 0.15.1
SymPy 0.7.6
"""
import logging
import warnings
import numpy
import pdb

#from numpy import __version__ as numpy_version
#from scipy import __version__ as scipy_version
#from sympy import __version__ as sympy_version
from scipy.interpolate import interp1d
from scipy.interpolate import BarycentricInterpolator
from scipy.interpolate import UnivariateSpline
from scipy.optimize import brentq, newton
from sympy import Interval, Intersection

from math import floor
from itertools import combinations
from warnings import warn

class NoIntersection(Exception):
    """
    An exception class to be raised when failed
    to find an intersection.
    """
    def __init__(self, value=''):
        self.value = value

    def __str__(self):
        return repr(self.value)

# NOTE: It seems that using HitTable takes too much RAM resource.
# We will not use it when finding joints of S-walls, instead we will
# find joints by first cutting S-walls into segments of monotonic
# curves and then find intersections between the segments.
class HitTable:
    """
    Construct a hash table of a coarse-grained plane.

    The table is implemented as a Python dictionary, whose item
    has a key of the form (x_ind, y_ind), i.e. (1, 2), and the
    corresponding value is a bin of the plane, which in turn
    is a dictionary of
       key: value =  ind: [[t1_i, t1_f], [t2_i, t2_f], ...],
    where c = list_of_curves[ind] gives a curve that hits the
    bin, c[tn_i] is the starting coordinate of the nth segment
    of the binned curve, and similarly c[tn_f] gives its endpoint.

    Parameters
    ----------
    bin_size : size of one square bin
    bin_fill_offset: When set to 0, each segment goes over the
        boundaries of bins, whereas a nonzero value corresponds to
        trimming each end of segment by that amount. A nonzero offset
        can result in missing an intersection when it is near the
        boundary of a bin.

    Methods
    -------
    get_bin_key(location): returns the key of the bin at the location.
    get_bin_location(bin_key): returns the center location of the bin
        corresponding to bin_key.
    put(bin_key, curve_index, segment_range): put a segment of a curve
        into the bin corresponding to bin_key.
    fill(list_of_curves): fill HitTable with the curves
        in list_of_curves.
    """

    def __init__(self, bin_size=1.0, bin_fill_offset=0):
        self._bin_size = bin_size
        self._bin_fill_offset = bin_fill_offset
        self._hit_table = {}

    def __getitem__(self, key):
        return self._hit_table[key]

    def __iter__(self):
        return self._hit_table.iterkeys()

    def __len__(self):
        return len(self._hit_table)

    def get_json_data(self):
        data = {}
        for bin_key, value in self._hit_table.iteritems():
            str_key = str(bin_key)
            data[str_key] = value
        return data

    def load_json_data(self, data):
        for str_key, value in data.iteritems():
            self._hit_table[eval(str_key)] = value

    def get_bin_size(self):
        return self._bin_size

    def get_bin_key(self, location):
        x_coord, y_coord = location

        bin_key_x, bin_key_y = map(
            lambda coord: int(floor(float(coord)/float(self._bin_size))),
            location
        )

        if not isinstance(bin_key_x, int):
            raise TypeError('bin_key_x = {} is not an integer'
                            '.'.format(bin_key_x))

        if not isinstance(bin_key_y, int):
            raise TypeError('bin_key_y = {} is not an integer'
                            '.'.format(bin_key_y))

        return (bin_key_x, bin_key_y)

    def get_bin_location(self, bin_key):

        location = map(
            lambda k: (k * self._bin_size + 0.5*self._bin_size),
            bin_key
        )

        return location

    def put(self, bin_key, curve_index, segment_range):

        if bin_key not in self._hit_table:
            bin_key_x, bin_key_y = bin_key
            # Check if bin_key is not a tuple of two integers.
            if (not isinstance(bin_key_x, int) or
                    not isinstance(bin_key_y, int)):
                raise KeyError('bin_key ' + str(bin_key) +
                               ' not a tuple of two integers.')
            # bin_key is legitimate; create an empty dict of curve segments.
            self._hit_table[bin_key] = {}
        if curve_index not in self._hit_table[bin_key]:
            # Check if curve_index is an integer
            if not isinstance(curve_index, int):
                raise KeyError('curve_index ' + str(curve_index) + ' not an '
                               'integer.')
            # Create an empty list for segments
            self._hit_table[bin_key][curve_index] = []

        self._hit_table[bin_key][curve_index].append(segment_range)

    def fill(self, curve_index, curve):
        """
        return a list of bin_key's that are put into the table.
        """
        t_i = t_0 = 0
        x_0, y_0 = curve[t_0]
        prev_bin_key = self.get_bin_key([x_0, y_0])
        new_bin_keys = []
        for t_n, [x_n, y_n] in enumerate(curve):
            try:
                bin_key = self.get_bin_key([x_n, y_n])
            except TypeError:
                raise
            # Cut the curve into segments where it goes over the current
            # bin or where it has a turning point.
            if prev_bin_key != bin_key:
                # curve[t_n - 1] and curve[t_n] are in different bins.
                # NOTE: bin_fill_offset = 0 gives segments that go over the
                # boundaries of bins. A nonzero offset can result in
                # missing an intersection when it is near the boundary
                # of a bin.
                try:
                    t_f = t_n - self._bin_fill_offset
                    if t_i > 0:
                        # Again, this makes the segment to go over the
                        # current bin by one step.
                        t_i -= 1
                    if t_i < t_f:
                        self.put(prev_bin_key, curve_index, [t_i, t_f])
                        new_bin_keys.append(prev_bin_key)
                except KeyError as e:
                    print str(e)
                    pass
                t_i = t_n
                prev_bin_key = bin_key
            elif is_turning_point(curve, t_n):
                # NOTE: when a curve has a turning point inside a bin,
                # there will be two segments in the bin from the same curve.
                # This is OK when we don't check a self-intersection
                # of a curve, which is fine when there is no branch cut
                # inside a bin.
                try:
                    self.put(prev_bin_key, curve_index, [t_i, t_n])
                    new_bin_keys.append(prev_bin_key)
                except KeyError as e:
                    print str(e)
                    pass
                t_i = t_n
            else:
                # curve[t_n] is not the end of the segment.
                continue
        # Take care of the last segment, if any.
        if t_i < t_n:
            try:
                self.put(prev_bin_key, curve_index, [t_i, t_n])
                new_bin_keys.append(prev_bin_key)
            except KeyError as e:
                print str(e)
                pass

        return new_bin_keys


# Only used by HitTable
def is_turning_point(curve, t):
    """Check whether curve[t] is the turning point or not."""
    t_max = len(curve) - 1
    if t_max < 2 or t == 0 or t == t_max:
        return False

    curve_near_t = curve[t-1:t+2]
    if len(curve_near_t) != 3:
        return False

    x, y = [list(c) for c in zip(*curve_near_t)]

    # Check if dx/dy = 0
    if (x[1] - x[0]) * (x[2] - x[1]) < 0:
        return True
    # Check if dy/dx = 0
    elif (y[1] - y[0]) * (y[2] - y[1]) < 0:
        return True
    else:
        return False
 

#def get_turning_points(curve):
#    """
#    Return a list of indices of turning points,
#    i.e. dx/dy = 0 or dy/dx = 0.
#
#    curve is (x_array, y_array), a tuple of NumPy arrays.
#    """
#    tps = []
#
#    if len(curve) < 3:
#        return tps
#
#    x, y = curve
#
#    if len(x) != len(y):
#        logging.error('Incorrect form of curve: len(x) = {}, len(y) = {}'
#                      .format(len(x), len(y)))
#        return tps
#
#    for t in range(1, len(x)-1):
#        if ((x[t] - x[t-1]) * (x[t+1] - x[t]) < 0 or
#            (y[t] - y[t-1]) * (y[t+1] - y[t]) < 0):
#            tps.append(t)
#
#    return tps
        

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
        x1_min, x1_max = numpy.argsort(x1)[[0, -1]]
        x2_min, x2_max = numpy.argsort(x2)[[0, -1]]

        y1_min, y1_may = numpy.argsort(y1)[[0, -1]]
        y2_min, y2_may = numpy.argsort(y2)[[0, -1]]

    x1_interval = Interval(x1_min, x1_max)
    x2_interval = Interval(x2_min, x2_max)

    y1_interval = Interval(y1_min, y1_may)
    y2_interval = Interval(y2_min, y2_may)

    x_range = x1_interval.intersect(x2_interval)
    y_range = y1_interval.intersect(y2_interval)

    return (x_range, y_range)


def find_intersection_of_segments(segment_1, segment_2, accuracy,
                                  bin_center=None, bin_size=None,
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

    if bin_center is not None and bin_size is not None:
        bin_center_x, bin_center_y = bin_center
        bin_x_interval = Interval(bin_center_x - 0.5*bin_size,
                                  bin_center_x + 0.5*bin_size)
        bin_y_interval = Interval(bin_center_y - 0.5*bin_size,
                                  bin_center_y + 0.5*bin_size)
        x_range = x_range.intersect(bin_x_interval)
        y_range = y_range.intersect(bin_y_interval)

    if (x_range.is_EmptySet or y_range.is_EmptySet or x_range.is_FiniteSet or
            y_range.is_FiniteSet):
        # The segments and the bin do not share a domain and therefore
        # there is no intersection.
        raise NoIntersection()

    f1 = interp1d(*segment_1)
    f2 = interp1d(*segment_2)
    delta_f12 = lambda x: f1(x) - f2(x)

    try:
        logging.debug('try brentq.')
        intersection_x = brentq(delta_f12, x_range.start, x_range.end)
        intersection_y = f1(intersection_x)
    except ValueError:
        """
        (f1 - f2) has the same sign at x_range.start & x_range.end
        use Newton's method instead, and for that purpose interpolate
        curves using polynomial interpolation.
        """
        # TODO: maybe spline interpolation is better because it can
        # provide derivatives of interpolatioin functions, but couldn't make
        # it work yet.
        """
        # NOTE: cubic spline interpolation requires more than 3 points.
        if len(segment_1) <= 3 or len(segment_2) <= 3:
            # not enough data; stop finding an intersection
            raise NoIntersection
        # UnivariateSpline requires x to be an increasing array.
        segment_1.sort(0)
        segment_2.sort(0)
        f1 = UnivariateSpline(*zip(*segment_1))
        f1_prime = f1.derivative()
        f2 = UnivariateSpline(*zip(*segment_2))
        f2_prime = f2.derivative()
        x0 = 0.5*(x_range.start + x_range.end)

        delta_f12 = lambda x: f1(x) - f2(x)
        delta_f12_prime = lambda x: f1_prime(x) - f2_prime(x)

        intersection_x = newton(delta_f12, x0, delta_f12_prime(x0))
        """
        logging.debug('try BarycentricInterpolator.')
        f1 = BarycentricInterpolator(*segment_1)
        f2 = BarycentricInterpolator(*segment_2)
        delta_f12 = lambda x: f1(x) - f2(x)

        x0 = 0.5*(x_range.start + x_range.end)

        try:
            logging.debug('try newton with x0 = %.8f.', x0)
            intersection_x = newton(delta_f12, x0, maxiter=newton_maxiter)
            logging.debug('intersection_x = %.8f.', intersection_x)
        except RuntimeError:
            # Newton's method fails to converge; declare no intersection
            raise NoIntersection()
        #except RuntimeWarning:
        #    pdb.set_trace()
        #    pass

        # Check if the intersection is within the curve range.
        # If not, the intersecion is not valid.
        if intersection_x not in x_range:
            raise NoIntersection()
        intersection_y = f1(intersection_x)
        if intersection_y not in y_range:
            raise NoIntersection()

    return [float(intersection_x), float(intersection_y)]
