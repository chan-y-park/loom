"""
A general intersection module.

Objects and functions to find intersections of real 1-dim curves
on a real 2-dim plane.

This module is based on the following libraries:
NumPy 1.9.0
SciPy 0.14.0
SymPy 0.7.5
"""
import logging
from numpy import __version__ as numpy_version
from scipy import __version__ as scipy_version
from sympy import __version__ as sympy_version
from scipy.interpolate import interp1d
from scipy.interpolate import BarycentricInterpolator
from scipy.interpolate import UnivariateSpline
from scipy.optimize import brentq, newton
# NOTE:this module requires SymPy 0.7.5.
from sympy import Interval, Intersection

from math import floor
from itertools import combinations
from warnings import warn

### temporarily muting these warnings, 
### they show up when I'm loading f_recover
### which is a bit annoying when doing debugging
### TO FIX

# Library version checks.
# if numpy_version < '1.9.0':
#     message = ('Current NumPy version ' + str(numpy_version) +
#                ' is lower than 1.9.0; '
#                'this module may not work properly.')
#     warn(message, Warning)

# if scipy_version < '0.14.0':
#     message = ('Current SciPy version ' + str(scipy_version) +
#                ' is lower than 0.14.0; '
#                'this module may not work properly.')
#     warn(message, Warning)

# if sympy_version < '0.7.5':
#     message = ('Current SymPy version ' + str(sympy_version) +
#                ' is lower than 0.7.5; '
#                'this module may not work properly.')
#     warn(message, Warning)


class NoIntersection(Exception):
    """
    An exception class to be raised when failed
    to find an intersection.
    """
    def __init__(self, value=''):
        self.value = value

    def __str__(self):
        return repr(self.value)


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

    def get_bin_size(self):
        return self._bin_size

    def get_bin_key(self, location):
        x_coord, y_coord = location

        bin_key_x, bin_key_y = map(
            lambda coord: int(floor(float(coord)/float(self._bin_size))),
            location
        )

        if not isinstance(bin_key_x, int):
            logging.debug('x_coord = %.8f, bin_size = %s',
                          x_coord, self._bin_size)
            raise TypeError('bin_key_x = {} is not an integer'
                            '.'.format(bin_key_x))

        if not isinstance(bin_key_y, int):
            logging.debug('y_coord = %.8f, bin_size = %s',
                          y_coord, self._bin_size)
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
        t_i = t_0 = 0
        x_0, y_0 = curve[t_0]
        prev_bin_key = self.get_bin_key([x_0, y_0])
        for t_n, [x_n, y_n] in enumerate(curve):
            try:
                bin_key = self.get_bin_key([x_n, y_n])
            except TypeError:
                logging.debug('t_n, x_n, y_n = %d, %.8f, %.8f',
                              t_n, x_n, y_n)
                raise
            # Cut the curve into a segment where it goes over the current
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
                except KeyError as e:
                    logging.debug('prev_bin_key, bin_key = %s, %s',
                                  prev_bin_key, bin_key)
                    print str(e)
                    pass
                t_i = t_n
                prev_bin_key = bin_key
            elif is_turning_point(curve, t_n):
                try:
                    self.put(prev_bin_key, curve_index, [t_i, t_n])
                except KeyError as e:
                    logging.debug('prev_bin_key, bin_key = %s, %s',
                                  prev_bin_key, bin_key)
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
            except KeyError as e:
                print str(e)
                pass


def is_turning_point(curve, t):
    """Check whether curve[t] is the turning point or not."""
    t_max = len(curve) - 1
    if t_max < 2 or t == 0 or t == t_max:
        return False

    x, y = [list(c) for c in zip(*curve[t-1:t+2])]

    # Check if dx/dy = 0
    if (x[1] - x[0]) * (x[2] - x[1]) < 0:
        return True
    # Check if dy/dx = 0
    elif (y[1] - y[0]) * (y[2] - y[1]) < 0:
        return True
    else:
        return False


def find_intersection_of_segments(segment_1, segment_2, bin_center, bin_size,
                                  newton_maxiter=5):
    """
    Find an intersection of two segments of curves in the same bin.

    First find interpolations of segments using scipy.interp1d and
    use SciPy's Brent method to find an intersection. When this
    doesn't work, use SciPy's polynomial interpolation and then use
    secant method to find an intersection.

    Parameters
    ----------
    segment_1, segment_2: Segments to find their intersection. Each
        segment is a list of [x, y].
    bin_center, bin_size: center location and the size of the bin
        that contains the segments.
    newton_maxiter: Maximum number of iterations for secant method.
        When increased, this gives a better accuracy of the
        intersection but it also greatly reduces the performance
        due to many cases when there is no intersection but
        the module keeps trying to find one.
    """
    # First check if the two segments share any x- and y-range.
    x1_i, y1_i = segment_1[0]
    x1_f, y1_f = segment_1[-1]

    x2_i, y2_i = segment_2[0]
    x2_f, y2_f = segment_2[-1]

    logging.debug('x1_i, x1_f = %.8f, %.8f', x1_i, x1_f)
    logging.debug('y1_i, y1_f = %.8f, %.8f', y1_i, y1_f)
    logging.debug('x2_i, x2_f = %.8f, %.8f', x2_i, x2_f)
    logging.debug('y2_i, y2_f = %.8f, %.8f', y2_i, y2_f)
    bin_center_x, bin_center_y = bin_center

    x1_interval = Interval(*sorted([x1_i, x1_f]))
    x2_interval = Interval(*sorted([x2_i, x2_f]))
    bin_x_interval = Interval(bin_center_x - 0.5*bin_size,
                              bin_center_x + 0.5*bin_size)

    y1_interval = Interval(*sorted([y1_i, y1_f]))
    y2_interval = Interval(*sorted([y2_i, y2_f]))
    bin_y_interval = Interval(bin_center_y - 0.5*bin_size,
                              bin_center_y + 0.5*bin_size)

    x_range = Intersection(x1_interval, x2_interval, bin_x_interval)
    y_range = Intersection(y1_interval, y2_interval, bin_y_interval)

    if (x_range.is_EmptySet or y_range.is_EmptySet or x_range.is_FiniteSet or
            y_range.is_FiniteSet):
        # The segments and the bin do not share a domain and therefore
        # there is no intersection.
        logging.debug('x_range = %s, y_range = %s', x_range, y_range)
        raise NoIntersection()

    logging.debug('x_range = [%.8f, %.8f]', x_range.start, x_range.end)
    logging.debug('y_range = [%.8f, %.8f]', y_range.start, y_range.end)

    f1 = interp1d(*zip(*segment_1))
    f2 = interp1d(*zip(*segment_2))
    delta_f12 = lambda x: f1(x) - f2(x)

    try:
        logging.debug('try brentq.')
        intersection_x = brentq(delta_f12, x_range.start, x_range.end)
        intersection_y = f1(intersection_x)
    except ValueError:
        logging.debug('f1(x_range.start), f1(x_range.end) = %.8f, %.8f',
                      f1(x_range.start), f1(x_range.end))
        logging.debug('f2(x_range.start), f2(x_range.end) = %.8f, %.8f',
                      f2(x_range.start), f2(x_range.end))
        logging.debug('delta_f12(x_range.start), delta_f12(x_range.end) '
                      '= %.8f, %.8f', delta_f12(x_range.start),
                      delta_f12(x_range.end))
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
        f1 = BarycentricInterpolator(*zip(*segment_1))
        f2 = BarycentricInterpolator(*zip(*segment_2))
        delta_f12 = lambda x: f1(x) - f2(x)

        x0 = 0.5*(x_range.start + x_range.end)

        try:
            logging.debug('try newton with x0 = %.8f.', x0)
            intersection_x = newton(delta_f12, x0, maxiter=newton_maxiter)
            logging.debug('intersection_x = %.8f.', intersection_x)
        except RuntimeError:
            # Newton's method fails to converge; declare no intersection
            raise NoIntersection()

        # Check if the intersection is within the curve range.
        # If not, the intersecion is not valid.
        if intersection_x not in x_range:
            raise NoIntersection()
        intersection_y = f1(intersection_x)
        if intersection_y not in y_range:
            raise NoIntersection()

    return [intersection_x, float(intersection_y)]
