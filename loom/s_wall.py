import logging
import numpy
import scipy
import sympy
# import cmath

import constants

from math import ceil
from cmath import exp, pi
from math import floor
from scipy import interpolate
from sympy import oo

from geometry import (
    find_xs_at_z_0, align_sheets_for_e_6_ffr, SHEET_NULL_TOLERANCE
)
from geometry import BranchPoint
from geometry import get_x

from ctypes_api import CTypesSWall
from ctypes_api import Message

from misc import (
    cpow, remove_duplicate, ctor2, r2toc, delete_duplicates, is_root,
    get_descendant_roots, sort_roots, n_nearest
)
from misc import nearest_index

use_numba = True
try:
    import numba
except ImportError, e:
    use_numba = False

x, z = sympy.symbols('x z')

# Number of x's at a fixed z
NUM_ODE_XS_OVER_Z = 2

# Desired precision on the phase of seeds
# XXX Warning: setting it too small will bring the seeding point
# too close to a branch point.
SEED_PHASE_PRECISION = 0.001
SEED_PRECISION_MAX_DEPTH = 5

# Minimum number of data points for each S-wall.
MIN_NUM_OF_DATA_PTS = 3


class Joint:
    def __init__(self, z=None, M=None, ode_xs=None, parents=None, roots=None,
                 logger_name='loom',):

        self.data_attributes = ['z', 'M', 'parents', 'label', 'roots',
                                'ode_xs']
        self.z = None
        self.M = None
        self.parents = None
        self.label = None
        self.roots = None
        self.ode_xs = None

        if z is None:
            # Return without setting attributes.
            return

        self.z = z
        self.M = M
        self.parents = parents
        self.roots = roots
        self.ode_xs = ode_xs

    def set_z_rotation(self, z_rotation):
        self.z *= complex(z_rotation)

    def __eq__(self, other):
        return self.label == other.label

    def get_json_data(self):
        json_data = {
            'z': ctor2(self.z),
            'M': self.M,
            'parents': [parent.label for parent in self.parents],
            'label': self.label,
            'ode_xs': [ctor2(x) for x in self.ode_xs],
        }
        if self.roots is not None:
            json_data['roots'] = [root.tolist() for root in self.roots]
        elif self.root is not None:
            json_data['roots'] = [self.root.tolist()]
        return json_data

    def set_from_json_data(self, json_data):
        self.z = r2toc(json_data['z'])
        # XXX: The following fixes the bug in the old data.
        if isinstance(json_data, list):
            self.M = json_data['M'][0]
        else:
            self.M = json_data['M']
        # NOTE: Joint.parents is loaded from JSON with labels of parents,
        # but is replaced with parent objects
        # in SpectralNetwork.set_from_json_data()
        self.parents = json_data['parents']
        self.label = json_data['label']
        # TODO: remove Joint.root & the following exception handling.
        try:
            self.roots = [numpy.array(root) for root in json_data['roots']]
        except KeyError:
            self.root = numpy.array(json_data['root'])
        self.ode_xs = [r2toc(x) for x in json_data['ode_xs']]

    def is_equal_to(self, other, accuracy, x_accuracy=None):
        if(abs(self.z - other.z) > accuracy):
            return False

        if(abs(self.M - other.M) > accuracy):
            return False

        if len(self.roots) != len(other.roots):
            return False

        for i in range(len(self.roots)):
            if numpy.array_equal(self.roots[i], other.roots[i]) is not True:
                return False

        # FIXME: Probably the following comparison is not needed
        # after comparing the roots of two S-walls.
        # Decide whether to do the comparison or not.
        if x_accuracy is None:
            x_accuracy = accuracy
        for i in range(NUM_ODE_XS_OVER_Z):
            if abs(self.ode_xs[i] - other.ode_xs[i]) > x_accuracy:
                return False

        return True


class SWall(object):
    def __init__(self, z_0=None, x_0=None, M_0=None, parents=None,
                 parent_roots=None, label=None, n_steps=None,
                 logger_name='loom'):
        """
        SWall.z is a NumPy array of length n_steps+1,
        where z[t] is the base coordinate.

        SWall.x is a Numpy array of the fiber coordinates at t, i.e.
        SWall.x[t] = [x[t][0], x[t][1], ...].
        (NOTE: for now, we just use two sheets, they serve only for
        numerical ode evolution. To get the sheets use the method
        SWall.get_sheets_at_t() instead.)
        """
        self.logger_name = logger_name
        # FIXME: self.zs & self.xs instead of z & x?
        if n_steps is None:
            self.z = []
            self.x = []
            self.M = []
        else:
            self.z = numpy.empty(n_steps + 1, numpy.complex128)
            self.x = numpy.empty(
                n_steps + 1, (numpy.complex128, NUM_ODE_XS_OVER_Z)
            )
            self.M = numpy.empty(n_steps + 1, numpy.float64)
            self.z[0] = z_0
            self.x[0] = x_0
            self.M[0] = M_0
        self.parents = parents
        self.parent_roots = parent_roots
        self.label = label

        # cuts_intersections = [[branch_point, i, '(c)cw'], ...]
        self.cuts_intersections = []
        # roots_basepoint = [t_0, z_0, x_0, root]
        self.roots_basepoint = []
        self.local_roots = []
        self.multiple_local_roots = None
        # local_weight_pairs is a list of pair of intgers.
        self.local_weight_pairs = []

        self.data_attributes = [
            'z', 'M', 'x', 'parents', 'parent_roots', 'label',
            'cuts_intersections',
            'local_roots', 'multiple_local_roots', 'local_weight_pairs',
        ]

    def is_trivialized(self):
        return (len(self.local_roots) > 0)

    def set_z_rotation(self, z_rotation):
        self.z *= complex(z_rotation)

    def __setitem__(self, t, data):
        """
        Set the data of the S-wall at t, where
            data = [z[t], x[t][0], x[t][1], ..., M[t]]
        """
        self.z[t] = data[0]
        self.x[t] = data[1:NUM_ODE_XS_OVER_Z + 1]
        self.M[t] = data[NUM_ODE_XS_OVER_Z + 1]

    def __getitem__(self, t):
        """
        Get the data of the S-wall at t, where
            data = [z[t], x[t][0], x[t][1], ..., M[t]]
        """
        return [self.z[t], self.x[t][0], self.x[t][1], self.M[t]]

    def resize(self, size):
        """
        Resize z & x arrays to discard garbage data.
        """
        self.z.resize((size))
        self.x.resize((size, NUM_ODE_XS_OVER_Z))
        self.M.resize((size))

    def get_json_data(self):
        json_data = {
            'z': numpy.array([self.z.real, self.z.imag]).T.tolist(),
            'M': self.M.tolist(),
            'x': numpy.rollaxis(
                numpy.array([self.x.real, self.x.imag]), 0, 3
            ).tolist(),
            'parents': [parent.label for parent in self.parents],
            'label': self.label,
            'cuts_intersections': [
                [br_loc.label, t, d]
                for br_loc, t, d in self.cuts_intersections
            ],
            'local_roots': [root.tolist() for root in self.local_roots],
            'local_weight_pairs': self.local_weight_pairs,
        }
        if self.parent_roots is not None:
            json_data['parent_roots'] = [
                root.tolist() for root in self.parent_roots
            ]

        if len(self.roots_basepoint) > 0:
            json_data['roots_basepoint'] = [
                self.roots_basepoint[0],
                ctor2(self.roots_basepoint[1]),
                [ctor2(x) for x in list(self.roots_basepoint[2])],
                list(self.roots_basepoint[3])
            ]
        if self.multiple_local_roots is not None:
            json_data['multiple_local_roots'] = [
                [root.tolist() for root in multiple_roots]
                for multiple_roots in self.multiple_local_roots
            ]
        return json_data

    def set_from_json_data(self, json_data):
        self.z = numpy.array([r2toc(z_t) for z_t in json_data['z']])
        # XXX: The following fixes the bug in the old data.
        if isinstance(json_data, list):
            self.M = numpy.array([M_t[0] for M_t in json_data['M']])
        else:
            self.M = numpy.array(json_data['M'])
        self.x = numpy.array(
            [[r2toc(x_i) for x_i in x_t] for x_t in json_data['x']]
        )
        # NOTE: SWall.parents is loaded from JSON with labels of parents,
        # but is replaced with parent objects
        # in SpectralNetwork.set_from_json_data()
        self.parents = json_data['parents']
        try:
            self.parent_roots = [
                numpy.array(root)
                for root in json_data['parent_roots']
            ]
        except KeyError:
            pass
        self.label = json_data['label']
        # NOTE: SWall.cuts_intersections is loaded from JSON
        # with labels of branch points,
        # but is replaced with branch point objects
        # in SpectralNetwork.set_from_json_data()
        self.cuts_intersections = json_data['cuts_intersections']
        self.local_roots = [
            numpy.array(root)
            for root in json_data['local_roots']
        ]
        try:
            self.multiple_local_roots = [
                [numpy.array(root) for root in multiple_roots]
                for multiple_roots in json_data['multiple_local_roots']
            ]
        except KeyError:
            pass
        self.local_weight_pairs = json_data['local_weight_pairs']

        try:
            roots_basepoint_data = json_data['roots_basepoint']
            self.roots_basepoint = [
                roots_basepoint_data[0],
                r2toc(roots_basepoint_data[1]),
                numpy.array(
                    [r2toc(x) for x in list(roots_basepoint_data[2])]
                ),
                numpy.array(roots_basepoint_data[3])
            ]
        except KeyError:
            pass

    def set_refs(self, obj_dict):
        self.parents = [
            obj_dict[parent_label]
            for parent_label in self.parents
        ]
        self.cuts_intersections = [
            [obj_dict[br_loc_label], t, d]
            for br_loc_label, t, d
            in self.cuts_intersections
        ]

    def Z(self):
        Z = 0
        for i in range(len(self.z) - 1):
            dz = self.z[i + 1] - self.z[i]
            Dx = self.x[i][0] - self.x[i][1]
            Z += Dx * dz
        return Z

    def is_primary(self):
        return (
            len(self.parents) == 1 and
            isinstance(self.parents[0], BranchPoint)
        )

    def get_splits(self, endpoints=False):
        splits = [t for bp, t, d in self.cuts_intersections]
        if endpoints is True:
            # Add values of t's of endpoints to the return value.
            return [0] + splits + [len(self.z) - 1]
        else:
            return splits

    def get_segments(self):
        splits = self.get_splits(endpoints=True)
        segments = []
        for i in range(len(splits) - 1):
            segments.append([splits[i], splits[i + 1]])
        return segments

    def downsample(self, ratio=None):
        # Resample data.
        n_data_pts = len(self.z)

        min_n_data_pts = constants.S_WALL_MIN_N_DATA_PTS
        if n_data_pts < min_n_data_pts:
            return None

        max_n_data_pts = constants.S_WALL_MAX_N_DATA_PTS
        if ratio is None:
            ratio = int(ceil(1.0 * n_data_pts /  max_n_data_pts))
        self.z = self.z[::ratio]
        self.x = self.x[::ratio]
        self.M = self.M[::ratio]

        # Adjust attributes accordingly.
        if len(self.cuts_intersections) > 0:
            cuts_intersections = []
            for bp, i, chi in self.cuts_intersections:
                cuts_intersections = [bp, i // ratio, chi]
            self.cuts_intersections = cuts_intersections
            self.enhance_at_cuts()

    def grow(
        self,
        branch_point_zs=[],
        puncture_point_zs=[],
        config=[],
        libs=None,
        use_scipy_ode=True,
        method=None,
        twist_lines=None,
    ):
        logger = logging.getLogger(self.logger_name)
        logger.info('Growing {}...'.format(self.label))

        bpzs = branch_point_zs
        ppzs = puncture_point_zs
        size_of_small_step = config['size_of_small_step']
        size_of_large_step = config['size_of_large_step']
        size_of_bp_neighborhood = config['size_of_bp_neighborhood']
        size_of_pp_neighborhood = config['size_of_pp_neighborhood']
        size_of_puncture_cutoff = config['size_of_puncture_cutoff']
        mass_limit = config['mass_limit']
        accuracy = config['accuracy']

        array_size = len(self.z)
        
        if method is not None:
            method_in_p_nbhd = constants.LIB_SCIPY_ODE
            method_out_p_nbhd = method
        else:
            method_in_p_nbhd = constants.LIB_SCIPY_ODE
            if use_scipy_ode:
                method_out_p_nbhd = constants.LIB_SCIPY_ODE
            elif libs.ctypes_s_wall.is_available():
                method_out_p_nbhd = constants.LIB_C
            elif libs.numba_grow is not None:
                method_out_p_nbhd = constants.LIB_NUMBA
            else:
                method_out_p_nbhd = constants.LIB_PYTHON

        current_method = None
        step = 0
        finished = False
        # Failure counter. Stop when this becomes less than zero.
        count = 3
        while(step < array_size and not finished and count > 0):
            logger.debug('step = {}'.format(step))

            z_i, _, _, _ = self[step]
            min_d_from_bps = get_min_d(z_i, bpzs)
            min_d_from_pps = get_min_d(z_i, ppzs)
            if (
                min_d_from_bps < size_of_bp_neighborhood or
                min_d_from_pps < size_of_pp_neighborhood
            ):
                if current_method is None:
                    current_method = method_in_p_nbhd
                stop_condition = constants.OUT_P_NBHD
            else:
                if current_method is None:
                    current_method = method_out_p_nbhd
                stop_condition = constants.IN_P_NBHD

            if current_method == constants.LIB_C:
                logger.debug(
                    'Growing {} using C libraries...'
                    .format(self.label)
                )
                msg = libs.ctypes_s_wall.message
                msg.s_wall_size = array_size
                msg.step = step
                msg.stop_condition = stop_condition 
                msg.rv = constants.NEWTON_MAX_STEPS

                libs.ctypes_s_wall.grow(
                    msg,
                    self.z,
                    self.x,
                    self.M,
                )

            elif current_method == constants.LIB_NUMBA:
                logger.debug(
                    'Growing {} using Numba...'
                    .format(self.label)
                )
                msg = libs.numba_grow(
                    self.z, self.x, self.M,
                    step=step,
                    stop_condition=stop_condition,
                    phi_k_czes=libs.phi_k_czes,
                    max_steps=constants.NEWTON_MAX_STEPS,
                    c_dz_dt=libs.c_dz_dt,
                    bpzs=bpzs,
                    ppzs=ppzs,
                    size_of_small_step=size_of_small_step,
                    size_of_large_step=size_of_large_step,
                    size_of_bp_neighborhood=size_of_bp_neighborhood,
                    size_of_pp_neighborhood=size_of_pp_neighborhood,
                    size_of_puncture_cutoff=size_of_puncture_cutoff,
                    mass_limit=mass_limit,
                    accuracy=accuracy,
                    twist_lines=libs.ctypes_s_wall.tl,
                )

            elif(
                current_method == constants.LIB_SCIPY_ODE or
                current_method == constants.LIB_PYTHON
            ):
                msg = Message()
                msg.s_wall_size = len(self.z)
                msg.step = step
                msg.stop_condition = stop_condition
                msg.rv = 0

                get_xs = libs.get_xs

                if current_method == constants.LIB_SCIPY_ODE:
                    ode = libs.ode
                    ode.set_initial_value(self[step])
                    name = 'SciPy ODE'
                else:
                    dz_dt = libs.dz_dt
                    name = 'Python libraries'

                logger.debug(
                    'Growing {} using {}...'
                    .format(self.label, name)
                )

                while step < (array_size - 1):
                    z_i, x_i_1, x_i_2, M_i = self[step]

                    min_d_from_bps = get_min_d(z_i, bpzs)
                    min_d_from_pps = get_min_d(z_i, ppzs)

                    if step >= MIN_NUM_OF_DATA_PTS:
                        if (min_d_from_pps < size_of_puncture_cutoff):
                            # Stop if z is inside a cutoff of a puncture.
                            msg.step = step
                            msg.rv = constants.NEAR_PUNCTURE
                            break
                        elif (
                            (mass_limit is not None) and
                            (M_i > mass_limit)
                        ):
                            # Stop if M exceeds mass limit.
                            msg.step = step
                            msg.rv = constants.MASS_LIMIT
                            break

                    # Adjust the step size if z is near a branch point.
                    step_size_factor = min([1.0, abs(x_i_1 - x_i_2)])
                    if min_d_from_bps < size_of_bp_neighborhood:
                        if stop_condition == constants.IN_P_NBHD:
                            msg.step = step
                            msg.rv = constants.IN_P_NBHD
                            break
                        else:
                            dt = size_of_small_step * step_size_factor
                    elif min_d_from_pps < size_of_pp_neighborhood:
                        if stop_condition == constants.IN_P_NBHD:
                            msg.step = step
                            msg.rv = constants.IN_P_NBHD
                            break
                        else:
                            dt = size_of_large_step * step_size_factor
                    else:
                        if stop_condition == constants.OUT_P_NBHD:
                            msg.step = step
                            msg.rv = constants.OUT_P_NBHD
                            break
                        else:
                            dt = size_of_large_step * step_size_factor

                    if (current_method == constants.LIB_SCIPY_ODE):
                        y_n = ode.integrate(ode.t + dt)

                        if not ode.successful():
                            msg.step = step
                            msg.rv = constants.ERROR_SCIPY_ODE
                            logger.warning(
                                '{} grow(): ode.integrate() failed at '
                                't = {}; z_i = {}, x_i = ({}, {}). '
                                .format(
                                    self.label, step,
                                    z_i, x_i_1, x_i_2,
                                )
                            )
                            break

                        z_n, x_n_1, x_n_2, M_n = y_n
                        # SciPy ODE returns M as a complex number.
                        M_n = M_n.real
                        y_n = [z_n, x_n_1, x_n_2, M_n]

                        # Check if this S-wall is near a twist line.
                        if twist_lines is not None:
                            if (
                                (z_i.imag * z_n.imag) < 0 and
                                twist_lines.contains(
                                    (z_i.real + z_n.real) * 0.5
                                )
                            ):
                                xs_at_z_n = get_xs(z_n)
                                i_1 = nearest_index(xs_at_z_n, (-1 * x_i_2))
                                i_2 = nearest_index(xs_at_z_n, (-1 * x_i_1))
                                x_n_1 = xs_at_z_n[i_1]
                                x_n_2 = xs_at_z_n[i_2]
                                y_n = [z_n, x_n_1, x_n_2, M_n]

                                if i_1 == i_2:
                                    msg.rv = constants.ERROR_SAME_XS
                                else:
                                    ode.set_initial_value(y_n)

                    elif(current_method == constants.LIB_PYTHON):
                        # Dx_i = x_i_1 - x_i_2
                        # z_n = z_i + dt * exp(
                        #     1j*(cmath.phase(libs.c_dz_dt) - cmath.phase(Dx_i))
                        # )
                        # M_n = M_i + abs(Dx_i) * dt
                        z_n = z_i + dt * dz_dt(z_i, x_i_1, x_i_2)
                        M_n = M_i + dt

                        if (
                            twist_lines is not None and
                            (z_i.imag * z_n.imag) < 0 and
                            twist_lines.contains(
                                (z_i.real + z_n.real) * 0.5
                            )
                        ):
                            x_i_1, x_i_2 = (-1 * x_i_2), (-1 * x_i_1)

                        xs_at_z_n = get_xs(z_n)
                        i_1 = nearest_index(xs_at_z_n, x_i_1)
                        i_2 = nearest_index(xs_at_z_n, x_i_2)

                        x_n_1 = xs_at_z_n[i_1]
                        x_n_2 = xs_at_z_n[i_2]
                        y_n = [z_n, x_n_1, x_n_2, M_n]

                        if i_1 == i_2:
                            msg.rv = constants.ERROR_SAME_XS

                    step += 1
                    self[step] = y_n

                    if msg.rv == constants.ERROR_SAME_XS:
                        msg.step = step
                        logger.warning(
                            '{} grow(): failed to get x\'s at '
                            't = {} near a twist line; '
                            'z_n = {}, x_i = ({}, {}), '
                            'xs_at_z_n = {}, i_1 == i_2 == {}.'
                            .format(
                                self.label, step,
                                z_n, x_i_1, x_i_2, xs_at_z_n, i_1,
                            )
                        )
                        break

                # End of inner while()

                if step == (array_size - 1):
                    msg.step = step
                    finished = True
            else:
                logger.warning('SWall.grow(): no grow method specified.')
                finished = True
                break

            logger.debug(
                'Growing {} stopped at t = {}: {}.'
                .format(self.label, msg.step, msg)
            )

            if msg.rv == 0:
                # Sucessfully finished.
                finished = True
                step = msg.step

            elif msg.out_p_nbhd():
                finished = False
                current_method = method_out_p_nbhd
                step = msg.step

            elif msg.in_p_nbhd():
                finished = False
                current_method = method_in_p_nbhd
                step = msg.step

            elif msg.near_puncture() or msg.mass_limit():
                finished = True
                new_size = msg.step + 1
                if new_size < array_size:
                    self.resize(new_size)
                break

            elif msg.error_same_xs():
                finished = False
                count -= 1
                logger.warning(
                    'Failed in growing {}: {} at t = {}, z = {}.'
                   .format(self.label, msg, msg.step, self.z[msg.step])
                )

                logger.warning('Try using SciPy ODE.')
                current_method = constants.LIB_SCIPY_ODE
                step = msg.step - 1

            elif msg.rv == constants.ERROR_SCIPY_ODE:
                finished = False
                count -= 1
                logger.warning('Try using a different method...')
                use_scipy_ode = False

            else:
                raise RuntimeError

        # End of outer while()

        if count == 0:
            logger.warning(
                '{}: SWall.grow() failed at step = {}'
                .format(self.label, step)
            )
            self.resize(step + 1)


    def determine_root_types(self, sw_data, cutoff_radius=0,):
        """
        1- Determine at which points the wall crosses a cut,
        for instance [55, 107, 231] would mean that
        it changes root-type 3 times.
        2- Then, pick a suitable point along the swall, away
        from branch points or singularities, and determine
        the root there.
        3- Finally, extend the determination
        of the root type to other segments by following
        the wall across the various splits induced by cuts,
        both forward and backwards, using the Weyl monodromy.
        """
        logger = logging.getLogger(self.logger_name)
        logger.info('Determining the root type of {}...'
                    .format(self.label))

        g_data = sw_data.g_data
        # branching will occur at branch points or irregular singularities
        branch_loci = sw_data.branch_points + sw_data.irregular_singularities
        # Adding a minimal radius of 1.0 is necessary
        # in case there is only a single branch point at z=0,
        # otherwise max_radius would be 0.
        max_radius = 2 * max(
            [abs(c_l.z) for c_l in branch_loci] + [1.0]
        )
        # Parametrize the z-coordinate of the k-wall's coordinates
        # as a function of the index.
        num_of_zs = len(self.z)
        traj_t = numpy.arange(num_of_zs)
        traj_z_r = self.z.real

        # Scan over branching loci cuts, see if path ever crosses one
        # based on x-coordinates only, and create a list of intersections,
        # each element being [br_loc, t, 'cw'|'ccw'].
        cuts_intersections = []
        for br_loc in branch_loci:
            br_loc_x = br_loc.z.real
            br_loc_y = br_loc.z.imag

            if num_of_zs > MIN_NUM_OF_DATA_PTS:
                # If the length of the S-wall's coordinates
                # is greater than 3, use the B-spline of SciPy
                # to find the intersections between cuts and the S-wall.
                g = interpolate.splrep(traj_t, traj_z_r - br_loc_x, s=0)

                # Produce a list of integers corresponding to indices of
                # the S-wall's coordinate list that seem to cross branch-cuts
                # based on the z-coordinate's real part.
                t_zeros = interpolate.sproot(g).tolist()

            elif 1 < num_of_zs <= MIN_NUM_OF_DATA_PTS:
                # There are two or three data points on the S-wall.
                # Use a linear interpolation for every pair of data points.
                t_zeros = []
                for i in range(num_of_zs - 1):
                    y_a, y_b = (traj_z_r - br_loc_x)[i:i + 2]
                    if y_a * y_b < 0:
                        t_zeros.append(y_a / (y_a - y_b))
            else:
                raise RuntimeError('{} has only one data point.'
                                   .format(self.label))

            intersection_ts = []
            for t_zero in t_zeros:
                t = int(floor(t_zero))
                # Remove a duplicate, if any.
                if (t in intersection_ts):
                    continue
                # Enforce imaginary-part of z-coordinate intersection
                # criterion: branch cuts extend vertically.
                if self.z[t].imag > br_loc_y:
                    intersection_ts.append(t)

            for t in intersection_ts:
                # Drop intersections of a primary S-wall with the
                # branch cut emanating from its parent branch-point
                # if such intersections happens within a short
                # distance from the starting point.
                if (
                    isinstance(br_loc, BranchPoint) and
                    br_loc == self.parents[0] and
                    (abs(br_loc.z - self.z[t]) < cutoff_radius)
                ):
                    continue

                # Check that the intersection actually happens
                # and is not an artifact of the interpolation used above
                # which could become an extrapolation
                x_p = self.z[t].real
                x_n = self.z[t + 1].real
                if not ((x_p < br_loc_x < x_n) or (x_n < br_loc_x < x_p)):
                    logger.warning(
                        '*** warning *** Drop a fake cut intersection.'
                    )
                    continue

                # Add
                # [the branch-point, t,
                #  the direction (either 'cw' or 'ccw')]
                # to each intersection.
                cuts_intersections.append(
                    [br_loc, t, clock(left_right(self.z, t))]
                )

        # Now sort intersections according to where they happen in proper
        # time; recall that the elements of cuts_intersections are
        # organized  as      [..., [branch_point_idx, t, 'ccw'] ,...]
        # where 't' is the integer of proper time at the intersection.
        self.cuts_intersections = sorted(
            cuts_intersections, cmp=lambda k1, k2: cmp(k1[1], k2[1])
        )
        logger.debug(
            'S-wall {} intersects the following cuts at the points\n{}.'
            .format(self.label, self.cuts_intersections)
        )

        # Add the actual intersection point to the S-wall
        # then update the attribute SWall.cuts_intersections accordingly
        #self.enhance_at_cuts(sw_data)
        self.enhance_at_cuts()

        # Choose a suitable point along the wall
        # we pick the one whose z coordinate's real part is
        # farthest from critical loci of the fibration, including
        # branch points and all types of punctures.
        # Ideally, we would like this basepoint to be not too far away
        # on the C-plane, because getting close to infinity
        # means colliding sheets usually.
        # But some walls are born beyond the max_radius in general
        # in that case, we just choose the t=0 coordinate
        # Also we keep t < 500, to avoid numerical errors induced by
        # numerics of ode_int, which tend to spoil the values of
        # s_wall.x[t] far along the trajectory.
        # FIXME: this fact of choosing t < 500 can cause trouble
        # if we have a primary S-wall that is almost vertical
        # and (running upwards), because then its
        # points are very close to a branch cut, and
        # it is better to take t as large as possible to take
        # advantage of the distance from the cut, to compute the root.
        # This can be handled by checking explicitly for such cases.
        t_0 = sorted(
            ([
                [t_i, min(
                    z_r_distance_from_critical_loci(z_i, sw_data)
                )]
                for t_i, z_i in enumerate(self.z) if (
                    (abs(z_i) < max_radius or t_i == 0) and t_i < 500
                )
            ]), cmp=lambda a, b: cmp(a[1], b[1])
        )[-1][0]

        z_0 = self.z[t_0]
        xs_0 = self.x[t_0]

        # Determine the initial root-type
        initial_root = get_s_wall_root(z_0, xs_0, sw_data,)

        if is_root(initial_root, sw_data.g_data) is False:
            logging.warning(
                'Could not assign a root to {}'
                .format(self.label)
            )

        self.roots_basepoint = [t_0, z_0, xs_0, initial_root]

        # A list of ordered pairs [...[i, j]...]
        # such that weights[j] - weights[i] = root
        initial_weight_pairs = (
            sw_data.g_data.ordered_weight_pairs(initial_root,)
        )

        self.local_roots = [initial_root]
        self.local_weight_pairs = [initial_weight_pairs]

        if len(self.cuts_intersections) > 0:
            # Note that we reverse the time-ordering!
            intersections_before_t_0 = [
                i_p for i_p in self.cuts_intersections if i_p[1] < t_0
            ][::-1]

            intersections_after_t_0 = [
                i_p for i_p in self.cuts_intersections if i_p[1] > t_0
            ]

            # Fill in the root types that occur after the basepoint
            for k in range(len(intersections_after_t_0)):
                br_loc, t, direction = intersections_after_t_0[k]

                current_root = self.local_roots[-1]
                new_root = g_data.weyl_monodromy(
                    current_root, br_loc, direction
                )
                new_weight_pairs = g_data.ordered_weight_pairs(new_root)

                self.local_roots.append(new_root)
                self.local_weight_pairs.append(new_weight_pairs)

            # Fill in the root types that occur before the basepoint
            # recall that their time-ordering has already been reversed
            # so the first one in the list is the closest to t_0, and so on
            for k in range(len(intersections_before_t_0)):
                br_loc, t, direction = intersections_before_t_0[k]

                current_root = self.local_roots[0]
                new_root = g_data.weyl_monodromy(
                    current_root, br_loc, direction, reverse=True
                )
                new_weight_pairs = g_data.ordered_weight_pairs(new_root)

                self.local_roots.insert(0, new_root)
                self.local_weight_pairs.insert(0, new_weight_pairs)

        # check that the roots obtained through comparison with the
        # trivialization coincides with the roots of the joint
        # or branch point sourcing the S-wall.
        if self.parent_roots is not None:
            if not (
                any(
                    (
                        (self.local_roots[0] == x).all() or
                        (self.local_roots[0] == -x).all()
                    )
                    for x in self.parent_roots
                )
            ):
                logger.warning(
                    'The root type of the S-wall is incompatible '
                    'with the sum of roots of its parents!'
                )
                return 'Rebuild S-wall'

        # Finished handling the single 'base' root of the S-wall.
        # Now handle the multiple local roots instead.
        root_0 = self.local_roots[0]
        if len(self.cuts_intersections) > 0:
            t_0 = int(floor(self.cuts_intersections[0][1] / 2))
        z_0 = self.z[t_0]
        ode_xs_0 = self.x[t_0]
        Dx_0 = ode_xs_0[0] - ode_xs_0[1]

        if (t_0 + 1) < len(self.z):
            ode_xs_1 = self.x[t_0 + 1]
        elif t_0 > 0:
            ode_xs_1 = self.x[t_0 - 1]
        else:
            RuntimeError('Length of z = {}: too short.'.format(len(self.z)))
        dx = max(abs(ode_xs_0[0] - ode_xs_1[0]),
                 abs(ode_xs_0[1] - ode_xs_1[1]),)

        ffr_xs_at_z_0 = sw_data.get_sheets_at_z(z_0, ffr=True).values()

        if self.parent_roots is None or len(self.parent_roots) == 1:
            self.multiple_local_roots = [[root] for root in self.local_roots]
        elif len(self.parent_roots) > 1:
            multiple_local_roots_0 = [root_0]
            all_roots_from_parents = (self.parent_roots +
                                      [-root for root in self.parent_roots])
            descendant_roots = get_descendant_roots(all_roots_from_parents,
                                                    sw_data.g_data,)
            all_roots_from_parents += descendant_roots
            for root in all_roots_from_parents:
                if numpy.array_equal(root_0, root):
                    continue
                else:
                    wall_weight_pairs = (
                        sw_data.g_data.ordered_weight_pairs(
                            root, ffr=True
                        )
                    )
                    ffr_w_p_0 = wall_weight_pairs[0]
                    ode_x1 = ffr_xs_at_z_0[ffr_w_p_0[0]]
                    ode_x2 = ffr_xs_at_z_0[ffr_w_p_0[1]]
                    Dx = ode_x1 - ode_x2
                    if abs(Dx - Dx_0) < dx:
                        multiple_local_roots_0.append(root)
            self.multiple_local_roots = [
                sort_roots(multiple_local_roots_0, sw_data.g_data)
            ]
            # We prepared all the base roots,
            # now we find how they change
            # as the S-wall crosses cuts.
            for k in range(len(self.cuts_intersections)):
                br_loc, t, direction = self.cuts_intersections[k]
                new_roots = []
                for current_root in self.multiple_local_roots[k]:
                    new_root = g_data.weyl_monodromy(
                        current_root, br_loc, direction
                    )
                    new_roots.append(new_root)
                self.multiple_local_roots.append(new_roots)
        else:
            raise RuntimeError('Parent of {} has no root.'
                               .format(self.label))

    def get_roots_at_t(self, t):
        """
        Given an integer t which parametrizes a point
        of the trajectory in proper time, return the local
        root at that point.
        """
        t_max = len(self.z) - 1
        if t < 0 or t > t_max:
            raise RuntimeError(
                'get_roots_at_t(): '
                't = {}, should be between 0 and {}.'
                .format(t, t_max)
            )
        else:
            closed_splits = self.get_splits() + [len(self.z) - 1]
            for i, sp in enumerate(closed_splits):
                if t <= sp:
                    if self.multiple_local_roots is not None:
                        return self.multiple_local_roots[i]
                    else:
                        return [self.local_roots[i]]
                    break
                else:
                    pass

    def get_weight_pairs_at_t(self, t):
        """
        Given an integer t which parametrizes a point
        of the trajectory in proper time, return the local
        pair of weights at that point.
        """
        t_max = len(self.z) - 1
        if t < 0 or t > t_max:
            raise RuntimeError(
                'get_weight_pairs_at_t(): '
                't = {}, should be between 0 and {}.'
                .format(t, t_max)
            )
        else:
            closed_splits = self.get_splits() + [len(self.z) - 1]
            for i, sp in enumerate(closed_splits):
                if t <= sp:
                    return self.local_weight_pairs[i]
                    break
                else:
                    pass

    def get_sheets_at_t(self, t, sw_data):
        """
        Given an integer t which parametrizes a point
        of the trajectory in proper time, return
        a list of pairs of values.
        For sheet labels [... [i, j] ...] at given t,
        return [... [x_i, x_j] ...]
        """
        z = self.z[t]
        # the following is a dictionary
        xs_at_z = sw_data.get_sheets_at_z(z)
        weight_pairs = self.get_weight_pairs_at_t(t)
        return [[xs_at_z[w_p[0]], xs_at_z[w_p[1]]] for w_p in weight_pairs]

    #def enhance_at_cuts(self, sw_data):
    def enhance_at_cuts(self):
        """
        Add the intersection points of Swalls and branch cuts
        also update the intersection data accordingly
        """
        if len(self.cuts_intersections) == 0:
            return None

        wall_pieces_z = []
        wall_pieces_x = []
        wall_pieces_M = []

        # split the wall into pieces, at the end of each
        # piece add the corresponding intersection point
        t_0 = 0
        for int_data in self.cuts_intersections:
            br_loc, t, chi = int_data

            z_1 = self.z[t]
            z_2 = self.z[t + 1]
            z_to_add = get_intermediate_z_point(z_1, z_2, br_loc.z)

            xs_1 = self.x[t]
            xs_2 = self.x[t + 1]
            xs_to_add = [
                get_intermediate_value(xs_1[i], xs_2[i], z_1, z_2, z_to_add)
                for i in range(NUM_ODE_XS_OVER_Z)
            ]

            M_to_add = get_intermediate_value(
                self.M[t], self.M[t + 1], z_1.real, z_2.real, z_to_add.real
            )

            z_piece = numpy.concatenate(
                (self.z[t_0:t + 1], numpy.array([z_to_add], dtype=complex))
            )
            x_piece = numpy.concatenate(
                (self.x[t_0:t + 1], numpy.array([xs_to_add], dtype=complex))
            )
            M_piece = numpy.concatenate(
                (self.M[t_0:t + 1], numpy.array([M_to_add], dtype=complex))
            )
            wall_pieces_z.append(z_piece)
            wall_pieces_x.append(x_piece)
            wall_pieces_M.append(M_piece)
            t_0 = t + 1

        # Get the last piece of the wall.
        wall_pieces_z.append(self.z[t_0:])
        wall_pieces_x.append(self.x[t_0:])
        wall_pieces_M.append(self.M[t_0:])

        # Make updates.
        self.z = numpy.concatenate(wall_pieces_z)
        self.x = numpy.concatenate(wall_pieces_x)
        self.M = numpy.concatenate(wall_pieces_M)

        # Update the intersection data.
        new_cuts_intersections = []
        for i, int_data in enumerate(self.cuts_intersections):
            br_loc, t_old, chi = int_data
            t_new = t_old + i + 1
            new_cuts_intersections.append([br_loc, t_new, chi])
        self.cuts_intersections = new_cuts_intersections

    def get_generation(self, generation=None):
        """
        Return the generation of this S-wall,
        starting with 1 when this is a primary S-wall.
        """
        generations = []
        if generation is None:
            generation = 0

        if len(self.parents) == 0:
            raise RuntimeError

        generation += 1
        for parent in self.parents:
            if isinstance(parent, BranchPoint):
                generations.append(generation)
            else:
                generations.append(parent.get_generation(generation))

        return max(generations)


# End of the definition of class SWall


class GrowLibs:
    def __init__(
        self,
        config=None,
        sw_data=None,
        phase=None,
        logger_name='loom',
    ):
        # logger = logging.getLogger(logger_name)
        # self.logger_name = logger_name
        self.f = sw_data.ffr_curve.num_eq
        self.v = sw_data.diff.num_v
        self.phase = phase
        self.accuracy = config['accuracy']
        self.c_dz_dt = None
        self.phi_k_czes = None
        # Method using SciPy ODE solver.
        self.ode = None
        # Method using Python routines.
        self.dz_dt = None
        self.get_xs = sw_data.ffr_curve.get_xs
        # Method using C libraries.
        self.ctypes_s_wall = None
        # Method using numba
        self.numba_grow = None

        self.default_lib = None

        x, z = sympy.symbols('x z')

        # NOTE: The following assumes that \lambda = x dz.
        c_v = self.v.coeff(x, n=1)
        self.c_dz_dt = complex(exp(phase * 1j) / c_v)

        self.phi_k_czes = sw_data.ffr_curve.get_phi_k_czes()

        # v = sympy.lambdify((z, x), self.v)
        def dz_dt(z, x1, x2):
            # Dv = (v(z, x1) - v(z, x2))
            Dv = x1 - x2
            if abs(Dv) < self.accuracy:
                raise RuntimeError(
                    'dz_dt(): Dv is too small, Dv={}, accuracy={}.'
                    .format(Dv, self.accuracy)
                )
            return self.c_dz_dt / Dv
        self.dz_dt = dz_dt

        df_dz = self.f.diff(z)
        df_dx = self.f.diff(x)
        # NOTE: F = -(\partial f / \partial z) / (\partial f / \partial x).
        F = sympy.lambdify((z, x), sympy.simplify(-df_dz / df_dx))

        def ode_f(t, z_x1_x2_M):
            z_i = z_x1_x2_M[0]
            x1_i = z_x1_x2_M[1]
            x2_i = z_x1_x2_M[2]
            dz_i_dt = dz_dt(z_i, x1_i, x2_i)
            dx1_i_dt = F(z_i, x1_i) * dz_i_dt
            dx2_i_dt = F(z_i, x2_i) * dz_i_dt
            dM_dt = 1
            return [dz_i_dt, dx1_i_dt, dx2_i_dt, dM_dt]

        self.ode = scipy.integrate.ode(ode_f)
        self.ode.set_integrator('zvode')

        self.ctypes_s_wall = CTypesSWall(
            config=config,
            sw_data=sw_data,
            phase=phase,
            c_dz_dt=self.c_dz_dt,
            phi_k_czes=self.phi_k_czes,
            logger_name=logger_name,
        )

        if use_numba:
            self.numba_grow = _grow
            # XXX Failed to compile a Numba function from _grow, try later.
            # self.numba_grow = numba.jit(nopython=True)(_grow)
            #self.numba_f_df_at_zx = numba.jit(nopython=True)(_f_df_at_zx)
            #self.numba_get_x = numba.jit(nopython=True)(_get_x)
        else:
            self.numba_grow = None

        if config['use_scipy_ode']:
            self.default_lib = constants.LIB_SCIPY_ODE
        elif self.ctypes_s_wall.is_available():
            self.default_lib = constants.LIB_C
        elif self.numba_grow is not None:
            self.default_lib = constants.LIB_NUMBA
        else:
            # Default S-wall grow library uses SciPy ODE.
            self.default_lib = constants.LIB_SCIPY_ODE


# XXX: Numba JIT complier fails to compile the following,
def _grow(
    zs, xs, Ms,
    step=0,
    stop_condition=None,
    phi_k_czes=None,
    max_steps=constants.NEWTON_MAX_STEPS,
    c_dz_dt=None,
    bpzs=None,
    ppzs=None,
    size_of_small_step=None,
    size_of_large_step=None,
    size_of_bp_neighborhood=None,
    size_of_pp_neighborhood=None,
    size_of_puncture_cutoff=None,
    mass_limit=None,
    accuracy=None,
    twist_lines=None,
):
    numba_rv = 0
    N, phi_k_n_czes, phi_k_d_czes = phi_k_czes
    array_size = len(zs)
    while step < (array_size - 1):
        z_i = zs[step]
        x_i_1, x_i_2 = xs[step]
        M_i = Ms[step]
        
        min_d_from_pps = constants.P_INF
        for ppz in ppzs:
            d = abs(z_i - ppz)
            if d < min_d_from_pps:
                min_d_from_pps = d

        min_d_from_bps = constants.P_INF
        for bpz in bpzs:
            d = abs(z_i - ppz)
            if d < min_d_from_bps:
                min_d_from_bps = d

        if step >= MIN_NUM_OF_DATA_PTS:
            # Stop if z is inside a cutoff of a puncture.
            if min_d_from_pps < size_of_puncture_cutoff:
                numba_rv = constants.NEAR_PUNCTURE
                break

            # Stop if M exceeds mass limit.
            if (mass_limit is not None) and (M_i > mass_limit):
                numba_rv = constants.MASS_LIMIT
                break

        # Adjust the step size if z is near a branch point.
        Dx_i = x_i_1 - x_i_2
        f_dt = abs(Dx_i)
        if f_dt > 1.0 : f_dt = 1.0
        if (
            min_d_from_bps < size_of_bp_neighborhood or
            min_d_from_pps < size_of_pp_neighborhood
        ):
            if stop_condition == constants.IN_P_NBHD:
                numba_rv = constants.IN_P_NBHD
                break
            else:
                dt = size_of_small_step * f_dt
        else:
            if stop_condition == constants.OUT_P_NBHD:
                numba_rv = constants.OUT_P_NBHD
                break
            else:
                dt = size_of_large_step * f_dt

        count = 0
        same_xs = False
        while count < constants.SAME_XS_MAX_STEPS:
            # z_n = z_i + dt * exp(1j * (cmath.phase(c_dz_dt) - 
            #                      cmath.phase(Dx_i)))
            # M_n = M_i + abs(Dx_i) * dt
            z_n = z_i + dt * c_dz_dt / Dx_i / (2.0 ** count)
            M_n = M_i + dt / (2.0 ** count)
#            import warnings
#            warnings.filterwarnings('error')
#            try:
#                z_n = z_i + dt * c_dz_dt / Dx_i / (2.0 ** count)
#                M_n = M_i + dt / (2.0 ** count)
#            except RuntimeWarning:
#                import pdb
#                pdb.set_trace()

            if (twist_lines is not None and (z_i.imag * z_n.imag) < 0):
                avg_z_r = (z_i.real + z_n.real) * 0.5
                for s, e in twist_lines:
                    if (s <= avg_z_r and avg_z_r <= e):
                        x_i_1, x_i_2 = (-1 * x_i_2), (-1 * x_i_1)

            x_n_1 = get_x(
                N, phi_k_n_czes, phi_k_d_czes, z_n, x_i_1, accuracy,
                max_steps=max_steps,
            )
            x_n_2 = get_x(
                N, phi_k_n_czes, phi_k_d_czes, z_n, x_i_2, accuracy,
                max_steps=max_steps,
            )
            if abs(x_n_1 - x_n_2) < accuracy:
                same_xs = True
            else:
                break
            count += 1

        step += 1
        zs[step] = z_n
        Ms[step] = M_n
        xs[step] = (x_n_1, x_n_2)

        if same_xs and count >= constants.SAME_XS_MAX_STEPS:
#            def __get_x(z, x, max_steps):
#                return get_x(
#                    N, phi_k_n_czes, phi_k_d_czes, complex(z), 
#                    complex(x), accuracy, max_steps=int(max_steps),
#                )
#            import pdb
#            pdb.set_trace()
            numba_rv = constants.ERROR_SAME_XS
            break

    msg = Message()
    msg.s_wall_size = array_size
    msg.step = step
    msg.stop_condition = stop_condition
    msg.rv = numba_rv

    return msg 


# XXX: Numba JIT complier fails to compile the following,
# Left for future use.
# if use_numba:
#     numba_grow = numba.jit(nopython=True)(_grow)

def _get_min_d(z, pzs):
    min_d = constants.P_INF
    for pz in pzs:
        d = abs(z - pz)
        if d < min_d:
            min_d = d
    return min_d


if use_numba:
    get_min_d = numba.jit(nopython=True)(_get_min_d)
else:
    get_min_d = _get_min_d


def get_s_wall_root(z, ffr_xs, sw_data):
    x_i, x_j = ffr_xs

    # Recall that S-wall numerical evolution
    # is based on the first fundamental representation.
    # In particular, the xs above are values of sheets
    # from the 1st fundamental rep cover, and should
    # be compared with the corresponding trivialization.
    # This is taken care of by setting the key argument
    # ffr=True when calling get_sheets_at_z.
    # The following is a dictionary
    sheets_at_z = sw_data.get_sheets_at_z(z, ffr=True)
    xs_at_z = sheets_at_z.values()

    if (
        (sw_data.g_data.type == 'D' or sw_data.g_data.type == 'E') and
        (abs(x_i) < SHEET_NULL_TOLERANCE or abs(x_j) < SHEET_NULL_TOLERANCE)
    ):
        if sw_data.g_data.type == 'D':
            n_w = 2
        elif sw_data.g_data.type == 'E':
            n_w = 3
        if abs(x_i) < SHEET_NULL_TOLERANCE:
            # Several sheets matching x_i
            closest_to_x_i = sorted(
                [[k, x_k] for k, x_k in enumerate(xs_at_z)],
                key=lambda y: abs(y[1] - x_i)
            )[0:n_w]
            i_s = [y[0] for y in closest_to_x_i]
            # Sheet matching x_j
            closest_to_x_j = sorted(xs_at_z, key=lambda x: abs(x - x_j))[0]
            j = (
                [k for k, v in sheets_at_z.iteritems()
                    if v == closest_to_x_j][0]
            )
            for k in i_s:
                # NOTE: Assuming that there is precisely ONE pair
                # whose difference is a root. In D-type it's really two pairs.
                # But either is fine, and they will be automatically
                # orthogonal to each other, so they are part of the same
                # "multi-root" set.
                # For E-type should introduce a little check to see if
                # there is actually more than one combo that works.
                alpha = (
                    sw_data.g_data.ffr_weights[j] -
                    sw_data.g_data.ffr_weights[k]
                )
                if is_root(alpha, sw_data.g_data):
                    i = k
                    break

        elif abs(x_j) < SHEET_NULL_TOLERANCE:
            # Sheet matching x_i
            closest_to_x_i = sorted(xs_at_z, key=lambda x: abs(x - x_i))[0]
            i = (
                [k for k, v in sheets_at_z.iteritems()
                    if v == closest_to_x_i][0]
            )
            # Several sheets matching x_j
            closest_to_x_j = sorted(
                [[k, x_k] for k, x_k in enumerate(xs_at_z)],
                key=lambda y: abs(y[1] - x_j)
            )[0:n_w]
            j_s = [y[0] for y in closest_to_x_j]
            for k in j_s:
                # NOTE: Assuming that there is precisely ONE pair
                # whose difference is a root. In D-type it's really two pairs.
                # But either is fine, and they will be automatically
                # orthogonal to each other, so they are part of the same
                # "multi-root" set.
                # For E-type should introduce a little check to see if
                # there is actually more than one combo that works.
                alpha = (
                    sw_data.g_data.ffr_weights[k] -
                    sw_data.g_data.ffr_weights[i]
                )
                if is_root(alpha, sw_data.g_data):
                    j = k
                    break

    else:
        # Sheet matching x_i
        closest_to_x_i = sorted(xs_at_z, key=lambda x: abs(x - x_i))[0]
        i = [k for k, v in sheets_at_z.iteritems() if v == closest_to_x_i][0]
        # Sheet matching x_j
        closest_to_x_j = sorted(xs_at_z, key=lambda x: abs(x - x_j))[0]
        j = [k for k, v in sheets_at_z.iteritems() if v == closest_to_x_j][0]

    return sw_data.g_data.ffr_weights[j] - sw_data.g_data.ffr_weights[i]


def get_s_wall_seeds(sw, theta, branch_point, config, logger_name):
    """
    S-walls are seeded from branch points.
    Each branch point has a number of ramification
    points lying above it.
    Regardless of the representation, it is sufficient
    to consider one of these ramification points
    to extract the seed data.
    We thus stick to (any)one ramification point of the
    fundamental representation to get the seeds.
    It is assumed that we are working in the first
    fundamental representation of the Lie algebra.
    """
    logger = logging.getLogger(logger_name)

    # FIXME: reintroduce the handling of massless punctures
    # see previous versions of this function, left above in comment.

    accuracy = config['accuracy']
    # initial_seed_size = config['size_of_small_step']
    initial_seed_size = accuracy 
    seeds = []
    min_dt = 1.0

    for rp in branch_point.ffr_ramification_points:
        z_0 = rp.z
        x_0 = rp.x
        r_i = rp.i
        rp_type = rp.ramification_type
        sw_diff_coeff = rp.sw_diff_coeff

        logger.debug('Finding seeds at ramification point (z,x)={}:'
                     .format([z_0, x_0]))
        logger.debug('\tramification index = {}'.format(r_i))
        logger.debug('\tramification type = {}'.format(rp_type))
        logger.debug(
            '\tleading coefficient of SW diff: {}'
            .format(sw_diff_coeff)
        )

        # Construct the seeding points for the branch point
        # by studying the type of ramification structure of the r.p.
        if rp_type == 'type_AD':
            continue
        elif rp_type == 'type_I':
            phases = [exp(2 * pi * 1j * float(i) / r_i) for i in range(r_i)]
            phi = [[p1 - p2 for p1 in phases] for p2 in phases]
            omega = exp(2.0 * pi * 1j * float(r_i) / float(r_i + 1))

            dz_phases = ([
                (1.0 / cpow((sw_diff_coeff), r_i, r_i + 1)) *
                exp(1j * theta * float(r_i) / (r_i + 1)) *
                ((-1.0 / phi[i][j]) ** (float(r_i) / (r_i + 1))) * (omega ** s)
                for i in range(r_i) for j in range(r_i)
                for s in range(r_i + 1) if i != j
            ])

            norm_dz_phases = [d / abs(d) for d in dz_phases]
            # these are the normalized phases of the seeds
            # with respect to the branch point:
            zetas = remove_duplicate(
                norm_dz_phases, lambda p1, p2: abs(p1 - p2) < (accuracy)
            )
        elif rp_type == 'type_II':
            if r_i % 2 == 1:
                raise RuntimeError(
                    'get_s_wall_seeds(): '
                    'Cannot have a type II ramification point '
                    'with odd ramification index.'
                )
            # defining this object just for enhanced readability of code
            # in comparing with notes on classification of ramifications
            r_k = r_i / 2
            phases = [
                exp(2 * pi * 1j * float(i) / (2.0 * r_k)) for i in range(r_k)
            ]
            phi = [[p1 - p2 for p1 in phases] for p2 in phases]
            psi = [[
                (phases[i] + phases[j]) * numpy.sign(i - j) for i in range(r_k)
            ] for j in range(r_k)]
            omega = exp(2.0 * pi * 1j * float(2 * r_k) / float(2 * r_k + 1))

            dz_phases = ([
                (1.0 / cpow(sw_diff_coeff, 2 * r_k, 2 * r_k + 1)) *
                exp(1j * theta * float(2 * r_k) / (2 * r_k + 1)) *
                ((-1.0 / phi[i][j]) ** (float(2 * r_k) / (2 * r_k + 1))) *
                (omega ** s)
                for i in range(r_k) for j in range(r_k)
                for s in range(2 * r_k + 1) if i != j
            ] + [
                (1.0 / cpow(sw_diff_coeff, 2 * r_k, 2 * r_k + 1)) *
                exp(1j * theta * float(2 * r_k) / (2 * r_k + 1)) *
                ((-1.0 / psi[i][j]) ** (float(2 * r_k) / (2 * r_k + 1))) *
                (omega ** s)
                for i in range(r_k) for j in range(r_k)
                for s in range(2 * r_k + 1) if i != j
            ])

            norm_dz_phases = [d / abs(d) for d in dz_phases]
            # these are the normalized phases of the seeds
            # with respect to the branch point:
            zetas = remove_duplicate(
                norm_dz_phases, lambda p1, p2: abs(p1 - p2) < accuracy
            )

        elif rp_type == 'type_III':
            if r_i % 2 == 1:
                raise RuntimeError(
                    'get_s_wall_seeds(): '
                    'Cannot have a type III ramification point '
                    'with odd ramification index.'
                )
            # defining this object just for enhanced readability of code
            # in comparing with notes on classification of ramifications
            r_k = r_i / 2

            phases = [
                exp(2 * pi * 1j * float(i) / (2.0 * (r_k - 1)))
                for i in range(r_k - 1)
            ] + [0.0]

            phi = [[p1 - p2 for p1 in phases] for p2 in phases]
            psi = [[
                (phases[i] + phases[j]) * numpy.sign(i - j)
                for i in range(r_k)
            ] for j in range(r_k)]

            omega = exp(
                2.0 * pi * 1j * float(2 * r_k - 2) / float(2 * r_k - 1)
            )

            dz_phases = ([
                (1.0 / cpow(sw_diff_coeff, 2 * r_k - 2, 2 * r_k - 1)) *
                exp(1j * theta * float(2 * r_k - 2) / (2 * r_k - 1)) *
                ((-1.0 / phi[i][j]) ** (float(2 * r_k - 2) / (2 * r_k - 1))) *
                (omega ** s)
                for i in range(r_k) for j in range(r_k)
                for s in range(2 * r_k - 1) if i != j
            ] + [
                (1.0 / cpow(sw_diff_coeff, 2 * r_k - 2, 2 * r_k - 1)) *
                exp(1j * theta * float(2 * r_k - 2) / (2 * r_k - 1)) *
                ((-1.0 / psi[i][j]) ** (float(2 * r_k - 2) / (2 * r_k - 1))) *
                (omega ** s)
                for i in range(r_k) for j in range(r_k)
                for s in range(2 * r_k - 1) if i != j
            ])

            norm_dz_phases = [d / abs(d) for d in dz_phases]
            # these are the normalized phases of the seeds
            # with respect to the branch point:
            zetas = remove_duplicate(
                norm_dz_phases, lambda p1, p2: abs(p1 - p2) < accuracy
            )

        elif rp_type == 'type_IV':
            # note: the following nomenclature of variables is
            # chose to match to the notes on seeding.
            # It would be better not to change it to avoid
            # potential confusion.

            # Note: the following are NOT the same as c_+/-
            # as in the notes, rather they are a 12-th root of those.
            sw_c_plus, sw_c_minus = sw_diff_coeff

            # computing the seed reference sheets -- see notes
            # for an explanation of what they are, and how they
            # will be employed
            seed_ref_xs = (
                [0.0, 0.0, 0.0] +
                [sw_c_plus * exp(2 * pi * 1j * float(i) / (12.0))
                 for i in range(12)] +
                [sw_c_minus * exp(2 * pi * 1j * float(i) / (12.0))
                 for i in range(12)]
            )

            aligned_seed_ref_xs = align_sheets_for_e_6_ffr(
                seed_ref_xs,
                sw.g_data.weights,
                near_degenerate_branch_locus=False
            )

            # now for each pair of reference sheets
            # which actually differs by a root, we consider
            # their difference to define a seed location
            delta_seed_xs = []
            for i, x_i in enumerate(aligned_seed_ref_xs):
                w_i = sw.g_data.weights[i]
                for j, x_j in enumerate(aligned_seed_ref_xs):
                    w_j = sw.g_data.weights[j]
                    if is_root(w_j - w_i, sw.g_data):
                        delta_seed_xs.append((x_j - x_i) / abs(x_j - x_i))

            omega = exp(2.0 * pi * 1j / 13.0)

            dz_phases = ([
                exp(1j * theta * 12.0 / 13.0) *
                ((-1.0 / delta_xs) ** (12.0 / 13.0)) *
                (omega ** s)
                for s in range(13) for delta_xs in delta_seed_xs
            ])

            norm_dz_phases = [d / abs(d) for d in dz_phases]
            # these are the normalized phases of the seeds
            # with respect to the branch point:
            zetas = remove_duplicate(
                norm_dz_phases, lambda p1, p2: abs(p1 - p2) < accuracy
            )

        # Now for each seeding point z_1 we identify two sheets
        # of the cover which match the phase of the displacement z_1-z_0
        for zeta in zetas:
            raise_precision = True
            precision_level = 0
            while (
                raise_precision is True and
                precision_level <= SEED_PRECISION_MAX_DEPTH
            ):
                logger.debug(
                    'Seeding precision level = {}'.format(precision_level)
                )
                dt = initial_seed_size * 0.1 ** (precision_level)
                if dt < min_dt:
                    min_dt = dt
                z_1 = z_0 + dt * zeta

                if rp_type == 'type_I':
                    all_x_s = find_xs_at_z_0(sw, z_1, x_0, r_i, ffr=True)
                    # just pick those sheets that are close enough
                    # to the ramification point
                    x_s = n_nearest(all_x_s, x_0, r_i)
                    # a list of the type
                    # [... [phase, [x_i, x_j]] ...]
                    x_i_x_j_phases = []
                    for i, x_i in enumerate(x_s):
                        for j, x_j in enumerate(x_s):
                            if i != j:
                                lambda_i = complex(
                                    sw.diff.num_v.subs(x, x_i).subs(z, z_1)
                                )
                                lambda_j = complex(
                                    sw.diff.num_v.subs(x, x_j).subs(z, z_1)
                                )
                                ij_factor = (
                                    -1.0 * exp(1j * theta) /
                                    (lambda_j - lambda_i)
                                )
                                x_i_x_j_phases.append(
                                    [(ij_factor) / abs(ij_factor), [x_i, x_j]]
                                )

                elif rp_type == 'type_II' or rp_type == 'type_III':
                    # we assume that the ramification index is maximal
                    # therefore we ask for all the sheets at z_1.
                    x_s = find_xs_at_z_0(sw, z_1, ffr=True)

                    # order of magnitude of expected separation
                    # of sheets at z_1
                    if rp_type == 'type_II':
                        dx = (
                            abs(cpow(sw_diff_coeff, 2 * r_k)) *
                            (dt ** (1.0 / float(r_i)))
                        )
                    elif rp_type == 'type_III':
                        dx = (
                            abs(cpow(sw_diff_coeff, 2 * r_k - 2)) *
                            (dt ** (1.0 / float(r_i)))
                        )
                    x_accuracy = min([accuracy, dx])

                    # a list of the type
                    # [... [phase, [x_i, x_j]] ...]
                    # where we eclude x_i=x_j and x_i=-x_j
                    # since in D-type there are no roots
                    # between a weight v and -v.
                    x_i_x_j_phases = []
                    for i, x_i in enumerate(x_s):
                        for j, x_j in enumerate(x_s):
                            if (
                                abs(x_j - x_i) > x_accuracy and
                                abs(x_j + x_i) > x_accuracy
                            ):
                                lambda_i = complex(
                                    sw.diff.num_v.subs(x, x_i).subs(z, z_1)
                                )
                                lambda_j = complex(
                                    sw.diff.num_v.subs(x, x_j).subs(z, z_1)
                                )
                                ij_factor = (
                                    -1.0 * exp(1j * theta) /
                                    (lambda_j - lambda_i)
                                )
                                x_i_x_j_phases.append(
                                    [(ij_factor) / abs(ij_factor), [x_i, x_j]]
                                )

                elif rp_type == 'type_IV':
                    x_s = find_xs_at_z_0(
                        sw, z_1, x_0, r_i,
                        ffr=True, use_sage=True
                    )
                    # Note: although we are working near a degenerate locus,
                    # it is correct to enforce the option
                    # near_degenerate_branch_locus=False
                    aligned_xs = align_sheets_for_e_6_ffr(
                        x_s,
                        sw.g_data.weights,
                        near_degenerate_branch_locus=False
                    )
                    # a list of the type
                    # [... [phase, [x_i, x_j]] ...]
                    x_i_x_j_phases = []
                    for i, x_i in enumerate(aligned_xs):
                        w_i = sw.g_data.weights[i]
                        for j, x_j in enumerate(aligned_xs):
                            w_j = sw.g_data.weights[j]
                            if is_root(w_j - w_i, sw.g_data):
                                lambda_i = complex(
                                    sw.diff.num_v.subs(x, x_i).subs(z, z_1)
                                )
                                lambda_j = complex(
                                    sw.diff.num_v.subs(x, x_j).subs(z, z_1)
                                )
                                ij_factor = (
                                    -1.0 * exp(1j * theta) /
                                    (lambda_j - lambda_i)
                                )
                                x_i_x_j_phases.append(
                                    [(ij_factor) / abs(ij_factor), [x_i, x_j]]
                                )

                closest_pair = sorted(
                    x_i_x_j_phases, key=lambda p: abs(p[0] - zeta)
                )[0][1]
                phase_mismatch = abs(sorted(
                    x_i_x_j_phases, key=lambda p: abs(p[0] - zeta)
                )[0][0] - zeta)
                if phase_mismatch < SEED_PHASE_PRECISION:
                    logger.debug(
                        'Reached desired precision for seed'
                        '\nThe Mismatch between its phase '
                        'and that of the displacement '
                        'is : {}'.format(phase_mismatch)
                    )
                    raise_precision = False
                else:
                    if precision_level < SEED_PRECISION_MAX_DEPTH:
                        precision_level += 1
                    else:
                        logger.warning(
                            'At {}, could not get the desired precision '
                            'on the seed of an S-wall.\nThe Mismatch between '
                            'the phase of a seed and that of the displacement '
                            'is : {}'.format(rp.label, phase_mismatch)
                        )
                        raise_precision = False
                        break
            M_0 = 0
            seeds.append([z_1, closest_pair, M_0])

    # for higher-index ramification points we need greater accuracy to
    # keep all the correct seeds, since dt is also their displacement
    # |z_1-z_0| we cannot just use dt, but must choose a small
    # fraction of it
    seeds = delete_duplicates(seeds, lambda s: s[0], accuracy=(min_dt / 100))
    logger.debug('Number of S-walls emanating = {}'.format(len(seeds)))
    logger.debug('these are the seeds {}\n'.format(seeds))
    branch_point.seeds = seeds
    return seeds


def z_r_distance_from_critical_loci(z, sw_data):
    critical_loci = (
        sw_data.branch_points + sw_data.irregular_singularities +
        sw_data.regular_punctures
    )
    return [abs(z.real - c_l.z.real)
            for c_l in critical_loci if c_l.z != oo]


def branch_locus_from_label(sw_data, br_loc_label):
    branch_loci = sw_data.branch_points + sw_data.irregular_singularities
    for br_loc in branch_loci:
        if br_loc.label == br_loc_label:
            return br_loc
    raise RuntimeError(
        'branch_locus_from_label(): '
        'Could not find any branching locus labeled {}'.format(br_loc_label)
    )


def clock(direction):
    if direction == 'left':
        return 'ccw'
    elif direction == 'right':
        return 'cw'
    else:
        raise RuntimeError('clock(): Cannot read direction.')


def left_right(l, point):
    """
    given the list
    l = [..., z, ...]
    and a point in the list (specified by the corresponding integer),
    determines whether x increases or decreases at that point,
    returning repsectively 'left' or 'right'
    """
    if point > len(l) - 1:
        raise RuntimeError(
            'left_right(): '
            'Cannot determine direction, '
            'point does not belong to the given list.'
        )
    elif point > 0:
        if l[point - 1].real < l[point].real:
            return 'right'
        else:
            return 'left'
    elif point == 0:
        if l[point].real < l[point + 1].real:
            return 'right'
        else:
            return 'left'


# TODO: Merge the following two functions?
def get_intermediate_z_point(z_1, z_2, bp_z_med):
    """
    get the intermediate point between z_1 and z_2
    in correspondence of the real part of z_med
    """
    # FIXME: division by zero may happend when x_1 and x_2 are too close.
    x_1 = z_1.real
    y_1 = z_1.imag
    x_2 = z_2.real
    y_2 = z_2.imag
    x_med = bp_z_med.real
    slope = (y_2 - y_1) / (x_2 - x_1)
    y_med = y_1 + slope * (x_med - x_1)
    return x_med + 1j * y_med


def get_intermediate_value(v_1, v_2, z_1, z_2, z_med):
    """
    Along an S-wall, get the intermediate value between v_1 and v_2
    in correspondence of z_med
    """
    # FIXME: division by zero may happend when z_1 and z_2 are too close.
    slope = (v_2 - v_1) / (z_2 - z_1)
    v_med = v_1 + slope * (z_med - z_1)
    return v_med
