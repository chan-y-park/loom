import logging
from math import floor


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
