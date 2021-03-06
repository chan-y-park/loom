import platform
import os
import logging
import ctypes
import numpy

from sympy import oo

import constants

numpy_ctypeslib_flags = ['C_CONTIGUOUS', 'ALIGNED']

array_1d_int = numpy.ctypeslib.ndpointer(
    dtype=numpy.int32,
    ndim=1,
    flags=numpy_ctypeslib_flags,
)

array_1d_float = numpy.ctypeslib.ndpointer(
    dtype=numpy.float64,
    ndim=1,
    flags=numpy_ctypeslib_flags,
)

array_2d_float = numpy.ctypeslib.ndpointer(
    dtype=numpy.float64,
    ndim=2,
    flags=numpy_ctypeslib_flags,
)

array_1d_complex = numpy.ctypeslib.ndpointer(
    dtype=numpy.complex128,
    ndim=1,
    flags=numpy_ctypeslib_flags,
)

array_2d_complex = numpy.ctypeslib.ndpointer(
    dtype=numpy.complex128,
    ndim=2,
    flags=numpy_ctypeslib_flags,
)


class Message(ctypes.Structure):
    _fields_ = [
        ('s_wall_size', ctypes.c_int),
        ('step', ctypes.c_int),
        ('stop_condition', ctypes.c_int),
        ('rv', ctypes.c_int),
    ]

    def error_same_xs(self):
        return self.rv == constants.ERROR_SAME_XS

    def near_puncture(self):
        return self.rv == constants.NEAR_PUNCTURE

    def mass_limit(self):
        return self.rv == constants.MASS_LIMIT

    def in_p_nbhd(self):
        return self.rv == constants.IN_P_NBHD

    def out_p_nbhd(self):
        return self.rv == constants.OUT_P_NBHD

    def __str__(self):
        if self.error_same_xs():
            msg = 'x1 == x2'
        elif self.near_puncture():
            msg = 'near a puncture'
        elif self.mass_limit():
            msg = 'reached the mass limit'
        elif self.in_p_nbhd():
            msg = 'inside the neighborhood of a point'
        elif self.out_p_nbhd():
            msg = 'outside the neighborhood of a point'
        elif self.rv == 0:
            msg = 'successfully finished'
        else:
            raise RuntimeError

        return msg


# class DiffParams(ctypes.Structure):
#     _fields_ = [
#         ('k', array_1d_int),
#         ('c', array_1d_complex),
#         ('e', array_1d_float),
#         ('n', ctypes.c_int),
#     ]


# class Points(ctypes.Structure):
#     _fields_ = [
#         ('z', array_1d_complex),
#         ('n', ctypes.c_int),
#     ]


class NumericalParameters(ctypes.Structure):
    _fields_ = [
        ('size_of_small_step', ctypes.c_double),
        ('size_of_large_step', ctypes.c_double),
        ('size_of_bp_neighborhood', ctypes.c_double),
        ('size_of_pp_neighborhood', ctypes.c_double),
        ('size_of_puncture_cutoff', ctypes.c_double),
        ('mass_limit', ctypes.c_double),
        ('accuracy', ctypes.c_double)
    ]


# class TwistLines(ctypes.Structure):
#     _fields_ = [
#         ('start', array_1d_float),
#         ('end', array_1d_float),
#         ('n', ctypes.c_int),
#     ]


class CTypesSWall:
    def __init__(
        self,
        config=None,
        sw_data=None,
        phase=None,
        c_dz_dt=None,
        phi_k_czes=None,
        logger_name='loom',
    ):
        self.clibs_s_wall_grow = None
        self.message = Message()

        self.N = None

        self.diff_n_k = None
        self.diff_n_c = None
        self.diff_n_e = None
        self.n_diff_n = None

        self.diff_d_k = None
        self.diff_d_c = None
        self.diff_d_e = None
        self.n_diff_d = None

        self.c_dz_dt = None

        self.bpz = None
        self.n_bpz = None

        self.ppz = None
        self.n_ppz = None

        self.np = NumericalParameters()

        self.tl = None
        self.n_tl = None

        logger = logging.getLogger(logger_name)
        try:
            clibs_s_wall = numpy.ctypeslib.load_library(
                's_wall',
                (os.path.dirname(os.path.realpath(__file__)) +
                 '/clibs/'),
            )
            self.clibs_s_wall_grow = clibs_s_wall.grow
        except OSError:
            logger.warning('Unable to load clib_s_wall_grow().')
            return None

        self.clibs_s_wall_grow.restype = ctypes.c_int

        self.clibs_s_wall_grow.argtypes = [
            ctypes.POINTER(Message),
            ctypes.c_int,       # int N
            array_1d_int,       # int* diff_n_k
            array_1d_complex,   # double complex* diff_n_c
            array_1d_float,     # double* diff_n_e
            ctypes.c_int,       # int n_diff_n
            array_1d_int,       # int* diff_d_k
            array_1d_complex,   # double complex* diff_d_c
            array_1d_float,     # double* diff_d_e
            ctypes.c_int,       # int n_diff_d
            array_1d_complex,   # double complex* c_dz_dt
            array_1d_complex,   # double complex* z
            array_2d_complex,   # double complex** x
            array_1d_float,     # double* M
            array_1d_complex,   # double complex* bpz
            ctypes.c_int,       # int n_bpz
            array_1d_complex,   # double complex* ppz
            ctypes.c_int,       # int n_ppz
            NumericalParameters,
            array_2d_float,     # twist_line* tl
            ctypes.c_int,       # int n_tl
        ]

        self.c_dz_dt = numpy.array([c_dz_dt], dtype=numpy.complex128)

        self.N, phi_k_n_czes, phi_k_d_czes = phi_k_czes

        self.diff_n_k = numpy.array(
            [k for k, c, e in phi_k_n_czes], dtype=numpy.int32
        )
        self.diff_n_c = numpy.array(
            [c for k, c, e in phi_k_n_czes], dtype=numpy.complex128
        )
        self.diff_n_e = numpy.array(
            [e for k, c, e in phi_k_n_czes], dtype=numpy.float64
        )
        self.n_diff_n = len(phi_k_n_czes)

        self.diff_d_k = numpy.array(
            [k for k, c, e in phi_k_d_czes], dtype=numpy.int32
        )
        self.diff_d_c = numpy.array(
            [c for k, c, e in phi_k_d_czes], dtype=numpy.complex128
        )
        self.diff_d_e = numpy.array(
            [e for k, c, e in phi_k_d_czes], dtype=numpy.float64
        )
        self.n_diff_d = len(phi_k_d_czes)

        bpzs = [
            p.z for p in
            #sw_data.ffr_ramification_points + sw_data.irregular_punctures
            sw_data.ffr_ramification_points
            if p.z != oo
        ]
        self.bpz = numpy.array(bpzs, dtype=numpy.complex128)
        self.n_bpz = int(len(bpzs))

        ppzs = [
            p.z for p in
            sw_data.regular_punctures + sw_data.irregular_punctures
            if p.z != oo
        ]
        self.ppz = numpy.array(ppzs, dtype=numpy.complex128)
        self.n_ppz = int(len(ppzs))

        self.np.size_of_small_step = config['size_of_small_step']
        self.np.size_of_large_step = config['size_of_large_step']
        self.np.size_of_bp_neighborhood = config['size_of_bp_neighborhood']
        self.np.size_of_pp_neighborhood = config['size_of_pp_neighborhood']
        self.np.size_of_puncture_cutoff = config['size_of_puncture_cutoff']
        self.np.mass_limit = config['mass_limit']
        self.np.accuracy = config['accuracy']

        config_twist_lines = config['twist_lines']
        if config_twist_lines is None:
            twist_lines = []
        else:
            twist_lines = eval(config_twist_lines)
        self.n_tl = len(twist_lines)
        self.tl = numpy.empty([len(twist_lines), 2], dtype=numpy.float64)
        for i, (s, e) in enumerate(twist_lines):
            if s == -oo:
                s = constants.N_INF
            if e == oo:
                e = constants.P_INF
            self.tl[i] = (s, e)

    def is_available(self):
        if self.clibs_s_wall_grow is None:
            return False
        else:
            return True

    def grow(self, msg, z, x, M):
        return self.clibs_s_wall_grow(
            msg,
            self.N,
            self.diff_n_k,
            self.diff_n_c,
            self.diff_n_e,
            self.n_diff_n,
            self.diff_d_k,
            self.diff_d_c,
            self.diff_d_e,
            self.n_diff_d,
            self.c_dz_dt,
            z,
            x,
            M,
            self.bpz,
            self.n_bpz,
            self.ppz,
            self.n_ppz,
            self.np,
            self.tl,
            self.n_tl,
        )


def libcgal_get_intersections():
    lib_name = 'libcgal_intersection'

    linux_distribution = platform.linux_distribution()[0]
    if linux_distribution == 'Ubuntu':
        lib_name += '_ubuntu'
    elif linux_distribution == 'debian':
        # NOTE: Anaconda Python returns 'debian' instead of 'Ubuntu'.
        lib_name += '_ubuntu'
    elif linux_distribution == 'Scientific Linux':
        lib_name += '_het-math2'
    else:
        raise OSError

    # Load CGAL shared library.
    libcgal_intersection = numpy.ctypeslib.load_library(
        lib_name,
        (os.path.dirname(os.path.realpath(__file__)) +
         '/cgal_intersection/'),
    )

    get_intersections = (libcgal_intersection.
                         find_intersections_of_curves)

    get_intersections.restype = ctypes.c_int
    get_intersections.argtypes = [
        array_1d_complex,
        ctypes.c_long,
        array_1d_complex,
        ctypes.c_long,
        array_2d_float, ctypes.c_int,
    ]

    return get_intersections
