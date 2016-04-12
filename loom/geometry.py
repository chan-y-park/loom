import sympy
import numpy
import logging
import copy
import pdb
import sympy.mpmath as mpmath

from sympy import oo, I
from sympy.mpmath import mp
# from sympy.mpmath.libmp.libhyper import NoConvergence
from itertools import combinations
from cmath import phase, pi
from matplotlib import cm as mpl_color_map

import sage_subprocess
from spectral_curve import get_ffr_curve_string
from misc import (ctor2, r2toc, PSL2C,
                  delete_duplicates, gather, parse_sym_dict_str,
                  n_remove_duplicate, spread_of_branch_points)

x, z = sympy.symbols('x z')

N_NULL_TRIPLES = 45
N_NULL_QUARTETS = 1008
NULL_TRIPLES_INDIVIDUAL = 5
NULL_QUARTETS_INDIVIDUAL = 72
SHEET_NULL_TOLERANCE = 0.001

ROOT_FINDING_MAX_STEPS = 50
ROOT_FINDING_PRECISION = 20

mp.dps = ROOT_FINDING_PRECISION

DEFAULT_LARGE_STEP_SIZE = 0.01

MIN_SPREAD = 0.1


class GData:
    """
    ffr_weights:
        List of weights of the FIRST fundamental representation,
            [v_0, ... , v_k , ...],
        where 'v_k' are numpy arrays corresponding to weights.
        - For g=A_n Lie algebras, the weights are given in IR^{n+1} as
            v_0 = (1,0,...,0),
            v_1 = (0,1,0,..), 
            ...,
            v_n = (0,...,0,1). 
          In this case, it does not matter how we identify weights
          with sheets, since the Weyl group acts by permuting all of 
          them freely.
        - For g=D_n, the weights are given in IR^{n} as
            v_0 = (1,0,...,0),
            v_1 = (0,1,...,0),
            v_{n-1} = (0,...,0,1),
            v_n = (-1,0,...,0),
            v_{n+1} = (0,-1,...,0),
            v_{2n-1} = (0,...,0,-1).
          In this case, we diivde the sheets into positive and negative ones,
          and assign the weights accordingly.
          The assignment of positive sheets is almost arbitrary: from each pair
          of positive/negative sheets one can pick either, as long as one makes
          an even number of "sign mistakes". We don't keep track of this,
          as a result there is an ambiguity in distinguishing one spinor 
          representation from the other
        - For other algebras, the order is imposed by SAGE.
    """
    def __init__(self, root_system=None, representation_str=None,
                 json_data=None, logger_name='loom'):
        self.logger_name = logger_name
        if json_data is not None:
            self.set_from_json_data(json_data)
        else:
            self.set_from_sage(root_system, representation_str)
        self.root_color_map = self.get_root_color_map()

        self.data_attributes = [
            'root_system', 'type', 'rank', 
            'fundamental_representation_index', 'highest_weight',
            'ffr_weights', 'roots', 'positive_roots', 'weights',
            'multiplicities', 'weight_basis', 'weight_coefficients',
            'root_color_map',
        ]

    def get_json_data(self):
        json_data = {
            'root_system': self.root_system,
            'type': self.type,
            'rank': self.rank,
            'fundamental_representation_index': (
                self.fundamental_representation_index
            ),
            'highest_weight': self.highest_weight,
            'ffr_weights': self.ffr_weights.tolist(),
            'roots': self.roots.tolist(),
            'positive_roots': self.positive_roots.tolist(),
            'weights': self.weights.tolist(),
            'multiplicities': self.multiplicities.tolist(),
            'weight_basis': self.weight_basis.tolist(),
            'weight_coefficients': self.weight_coefficients.tolist(),
        }

        return json_data

    def set_from_json_data(self, json_data):
        self.root_system = json_data['root_system']
        self.type = json_data['type']
        self.rank = json_data['rank']
        self.fundamental_representation_index = (
            json_data['fundamental_representation_index']
        )
        self.highest_weight = json_data['highest_weight']
        self.ffr_weights = numpy.array(json_data['ffr_weights'])
        self.roots = numpy.array(json_data['roots'])
        self.positive_roots = numpy.array(json_data['positive_roots'])
        self.weights = numpy.array(json_data['weights'])
        self.multiplicities = numpy.array(json_data['multiplicities'])
        self.weight_basis = numpy.array(json_data['weight_basis'])
        self.weight_coefficients = numpy.array(
            json_data['weight_coefficients']
        )

    def set_from_sage(self, root_system, representation_str):
        logger = logging.getLogger(self.logger_name)
        self.root_system = root_system
        # type is 'A', 'D', or 'E'.
        self.type = root_system[0]
        self.rank = eval(root_system[1:])

        representation = eval(representation_str)
        if isinstance(representation, int):
            # Representation is specified as an index
            # of a fundamental representation, i.e. n of \omega_n.
            self.fundamental_representation_index = representation
            self.highest_weight = [
                1 if i == (self.fundamental_representation_index - 1) else 0
                for i in range(self.rank)
            ]
        elif isinstance(representation, list):
            # Representation is specified in the coroot(Dynkin) basis.
            self.highest_weight = representation
            height = 0
            for i, n_i in enumerate(self.highest_weight):
                height += n_i
                if n_i == 1:
                    self.fundamental_representation_index = i + 1
            if height > 1:
                self.fundamental_representation_index = None 
                logger.warning('{} is not a fundamental representation.'
                               .format(representation_str))

        sage_data = sage_subprocess.get_g_data(
            root_system, 
            self.highest_weight,
        )

        self.ffr_weights = numpy.array(sage_data['ffr_weights'])
        self.roots = numpy.array(sage_data['roots'])
        self.positive_roots = numpy.array(sage_data['positive_roots'])
        self.weights = numpy.array(sage_data['weights'])
        self.multiplicities = numpy.array(sage_data['multiplicities'])
        self.weight_basis = numpy.array(sage_data['weight_basis'])
        # The i-th row of self.coefficients is the representation
        # of self.weights[i] in the self.basis.
        self.weight_coefficients = numpy.array(
            sage_data['weight_coefficients']
        )

    # XXX: rename this to get_ordered_weight_pairs()
    def ordered_weight_pairs(self, root, ffr=False):
        """
        Return list of pairs of weight indices.
        """
        pairs = []

        if ffr is False:
            weights = self.weights
        elif ffr is True:
            weights = self.ffr_weights

        for i, w_1 in enumerate(weights):
            for j, w_2 in enumerate(weights):
                # FIXME: Currently weights are arrays of floats,
                # the following comparison may have a numerical issue.
                if numpy.array_equal(w_2 - w_1, root):
                    pairs.append([i, j])

        return pairs

    # XXX: rename this to apply_weyl_monodromy()
    def weyl_monodromy(
        self, root, br_loc, direction, reverse=False, perm_matrix=None
    ):
        """
        Returns a new root of the segment of an S-wall
        when it crosses a branch cut from a brancing locus
        which could be a branch point or an irregular singularity.
        An option to use directly a sheet permutation matrix is also 
        given, by specifying the corresponding argument.
        """
        if perm_matrix is None:
            if reverse is False:
                if direction == 'ccw':
                    monodromy_matrix = br_loc.monodromy
                elif direction == 'cw':
                    monodromy_matrix = (
                        numpy.linalg.inv(br_loc.monodromy).astype(int)
                    )
            elif reverse is True:
                if direction == 'cw':
                    monodromy_matrix = br_loc.monodromy
                elif direction == 'ccw':
                    monodromy_matrix = (
                        numpy.linalg.inv(br_loc.monodromy).astype(int)
                    )
        else:
            if reverse is False:
                if direction == 'ccw':
                    monodromy_matrix = perm_matrix
                elif direction == 'cw':
                    monodromy_matrix = (
                        numpy.linalg.inv(perm_matrix).astype(int)
                    )
            elif reverse is True:
                if direction == 'cw':
                    monodromy_matrix = perm_matrix
                elif direction == 'ccw':
                    monodromy_matrix = (
                        numpy.linalg.inv(perm_matrix).astype(int)
                    )
                
        pair_0 = self.ordered_weight_pairs(root)[0]
        v_i_ind = pair_0[0]
        v_j_ind = pair_0[1]
        ordered_weights = self.weights
        new_v_i = sum([
            monodromy_matrix[k][v_i_ind] * v 
            for k, v in enumerate(ordered_weights)
        ])
        new_v_j = sum([
            monodromy_matrix[k][v_j_ind] * v 
            for k, v in enumerate(ordered_weights)
        ]) 

        new_root = new_v_j - new_v_i
        return new_root

    def get_root_color_map(self):
        """
        Create a mapping between a positive root and a color.
        A negative root will have the same color as the corresponding
        positive root.
        """
        p_roots = self.positive_roots
        num_p_roots = len(p_roots)
        p_root_colors = []
        for i in range(num_p_roots):
            try:
                color_map = mpl_color_map.viridis
            except AttributeError:
                color_map = mpl_color_map.jet
            r, g, b, alpha = color_map(
                (i / float(num_p_roots)), bytes=True
            )
            p_root_colors.append('#{:02x}{:02x}{:02x}'.format(r, g, b))

        return p_root_colors

    def get_root_color(self, root):
        logger = logging.getLogger(self.logger_name)

        for i, p_root in enumerate(self.positive_roots):
            if (
                numpy.array_equal(p_root, root) or
                numpy.array_equal(-p_root, root)
            ):
                return self.root_color_map[i]

        # If color is not assigned because this is not a true root
        # still give a color, so the plot can be displayed
        logger.warning('No color mapped for the root {}'
                       .format(root.tolist()))
        return self.root_color_map[0]


class RamificationPoint:
    def __init__(
        self, z=None, Ciz=None, x=None, i=None, label=None, json_data=None,
        is_puncture=False,
    ):
        if json_data is None:
            # z is the numerical value of the PSL2C-transformed z-coordinate.
            self.z = z
            # Ciz is the value of the z-coordinate 
            # before the PSL2C transformation.
            self.Ciz = Ciz
            self.x = x
            self.i = i
            self.label = label
            # For explanations on the following two attributes, 
            # see the function which analyzes ramification points
            # in the trivializatio module.
            self.ramification_type = None
            self.sw_diff_coeff = None
            self.sw_diff_coeffs_a_b = [None, None]
            self.is_puncture = is_puncture
        else:
            self.set_from_json_data(json_data)

        self.data_attributes = [
            'z', 'Ciz', 'x', 'i', 'label', 'ramification_type',
            'sw_diff_coeff', 'is_puncture',
        ]

    def set_z_rotation(self, z_rotation):
        if self.z != oo:
            self.z *= complex(z_rotation)

#    def __str__(self):
#        return 'z = {}, x = {}, i = {}'.format(self.z, self.x, self.i)

    def __eq__(self, other):
        return self.label == other.label

    def get_json_data(self):
        if self.sw_diff_coeff is None:
            sw_diff_coeff = None
        elif type(self.sw_diff_coeff) == list:
            # This handles the case of type IV ramification points
            sw_diff_coeff = [ctor2(x) for x in self.sw_diff_coeff]
        else:
            sw_diff_coeff = ctor2(self.sw_diff_coeff)
        json_data = {
            'z': ctor2(self.z),
            'Ciz': str(self.Ciz),
            'x': ctor2(self.x),
            'i': self.i,
            'label': self.label,
            'ramification_type': self.ramification_type,
            'sw_diff_coeff': sw_diff_coeff,
            'is_puncture': self.is_puncture
        }
        return json_data

    def set_from_json_data(self, json_data):
        self.z = r2toc(json_data['z'])
        self.Ciz = sympy.sympify(json_data['Ciz'])
        self.x = r2toc(json_data['x'])
        self.i = json_data['i']
        self.label = json_data['label']
        self.ramification_type = json_data['ramification_type']

        if json_data['sw_diff_coeff'] is not None:
            if self.ramification_type == 'type_IV':
                # In the case of type IV there are two coefficients to unpack
                sw_diff_coeff = [r2toc(x) for x in json_data['sw_diff_coeff']]
            else:
                sw_diff_coeff = r2toc(json_data['sw_diff_coeff'])
        else:
            sw_diff_coeff = None
        self.sw_diff_coeff = sw_diff_coeff 
        try:
            self.is_puncture = json_data['is_puncture']
        except KeyError:
            self.is_puncture = False


class Puncture:
    def __init__(self, z=None, Ciz=None, cutoff=None, label=None,
                 json_data=None):
        if json_data is None:
            # z is the numerical value of the PSL2C-transformed z-coordinate.
            self.z = z
            # Ciz is the value of the z-coordinate 
            # before the PSL2C transformation.
            self.Ciz = Ciz
            self.label = label
        else:
            self.set_from_json_data(json_data)

        self.data_attributes = ['z', 'Ciz', 'label']

    def set_z_rotation(self, z_rotation):
        if self.z != oo:
            self.z *= complex(z_rotation)

    def __eq__(self, other):
        return self.label == other.label

    def get_json_data(self):
        json_data = {
            'Ciz': str(self.Ciz),
            'label': self.label,
        }

        if self.z == oo:
            json_data['z'] = 'oo'
        else:
            json_data['z'] = ctor2(self.z)

        return json_data

    def set_from_json_data(self, json_data):
        if json_data['z'] == 'oo':
            self.z = oo
        else:
            self.z = r2toc(json_data['z'])

        self.Ciz = sympy.sympify(json_data['Ciz'])
        self.label = json_data['label']


class SWCurve:
    """
    SWCurve.sym_eq is the same curve as defined in the configuration,
    without the PSL2C transformation for the simplicity of its anlysis.

    SWCurve.num_eq is the curve with numerical coefficients &
    after the PSL2C transformation. This is the curve that is used
    in the numerical evolution of S-walls.
    """
    def __init__(self, casimir_differentials=None, g_data=None, 
                 diff_params=None, mt_params=None, 
                 z_rotation=sympy.sympify('1'), ffr=False,
                 expand=False):
        self.sym_eq = None
        self.num_eq = None

        if ffr is True:
            # Build a cover in the first fundamental representation.
            ffr_eq_str = get_ffr_curve_string(
                casimir_differentials, g_data.type, g_data.rank
            )

            try:
                if expand is True:
                    # TODO: SAGE is much faster in expanding an expression.
                    self.sym_eq = sympy.sympify(ffr_eq_str).expand()
                else:
                    self.sym_eq = sympy.sympify(ffr_eq_str)
            except:
                raise

            # NOTE: We apply PSL2C only to the numerical curve
            # for the simplicity of analysis.

            Ciz = PSL2C(mt_params, z_rotation * z, inverse=True) 
            self.num_eq = (
                self.sym_eq.subs(z, Ciz).subs(diff_params)
                .evalf(n=ROOT_FINDING_PRECISION, chop=True)
            )

        else:
            # TODO: Need to build a cover in a general representation
            # from the differentials, using symmetric polynomials.
            raise NotImplementedError(
                'class SWCurve with a general representation '
                'is not implemented yet.'
            )

        self.data_attributes = ['sym_eq', 'num_eq']

    def set_z_rotation(self, z_rotation):
        self.num_eq = (self.num_eq
                       .subs(z, z_rotation * z)
                       .evalf(n=ROOT_FINDING_PRECISION, chop=True))
        return self

    def get_xs(self, z_0, use_sage=False):
        """
        Return a numpy array of x-coordinates over z = z_0.
        """
        if self.num_eq is None:
            raise NotImplementedError

        if use_sage is False:
            fx = self.num_eq.subs(z, z_0)
            sym_poly = sympy.Poly(fx, x, domain='CC')
            coeff_list = map(complex, sym_poly.all_coeffs())
            return numpy.roots(coeff_list)
        else:
            f_x_eq = self.num_eq.subs(z, z_0).evalf(n=ROOT_FINDING_PRECISION)
            f_x_roots = sage_subprocess.solve_single_eq_single_var(
                f_x_eq,
                var='x',
                precision=ROOT_FINDING_PRECISION,
            )
            return map(complex, f_x_roots)


class SWDiff:
    def __init__(self, v_str, g_data=None, diff_params=None, mt_params=None,
                 z_rotation=sympy.sympify('1'),):
        # sym_v is a SymPy expression. 
        self.sym_v = sympy.sympify(v_str)

        # num_v is from sym_v with its parameters 
        # substituted with numerical values.
        # NOTE: We apply PSL2C only to the numerical curve
        # for the simplicity of analysis.
        Ciz = PSL2C(mt_params, z_rotation * z, inverse=True) 
        dCiz = Ciz.diff(z)
        # self.jac = dCiz
        self.num_v = (
            (self.sym_v.subs(z, Ciz) * dCiz).subs(diff_params)
            .evalf(n=ROOT_FINDING_PRECISION, chop=True)
        )

        self.data_attributes = ['sym_v', 'num_v']

    def set_z_rotation(self, z_rotation):
        self.num_v = (self.num_v
                      .subs(z, z_rotation * z)
                      .evalf(n=ROOT_FINDING_PRECISION, chop=True))
        return self


class SWDataBase(object):
    """
    A class containing the geometric data of a Seiberg-Witten curve
    in the first fundamental representation,
        \lambda^N + \sum_{k=2}^N a_k(z) dz^k \lambda^{N-k}, 
    where \lambda is the Seiberg-Witten differential of the form 
        \lambda = x dz.

    This class also includes ramification points of the curve, and 
    another curve in the representation given in the configuration.

    This is the base class of SWDataWithTrivialization, where 
    the trivialization data of the curve is contained.
    """
    def __init__(self, config, logger_name='loom', json_data=None,):
        self.logger_name = logger_name
        logger = logging.getLogger(self.logger_name)

        self.data_attributes = [
            'g_data', 'regular_punctures', 'irregular_punctures',
            'ffr_ramification_points', 'z_plane_rotation', 'accuracy',
            'ffr_curve', 'curve', 'diff',
        ]

        self.g_data = None
        self.regular_punctures = []
        self.irregular_punctures = []
        self.ffr_ramification_points = None
        # This rotation is applied when getting the trivialization.
        # TODO: When rotating the z-plane back to the original place,
        # set this to zero and introduce a new attribute to save
        # the angle of the branch cuts.
        self.z_plane_rotation = None
        self.accuracy = config['accuracy']

        self.ffr_curve = None
        self.curve = None
        self.diff = None

        if config['mt_params'] is not None:
            mt_params = sympy.sympify(config['mt_params'])
        else:
            mt_params = None

        casimir_differentials = {}
        for k, phi_k in parse_sym_dict_str(config['casimir_differentials']):
            casimir_differentials[eval(k)] = phi_k

        diff_params = {}
        for var, val in parse_sym_dict_str(config['differential_parameters']):
            diff_params[var] = sympy.sympify(val)

        if json_data is None:
            self.g_data = GData(root_system=config['root_system'],
                                representation_str=config['representation'],
                                logger_name=self.logger_name,)
            self.set_from_config(
                config,
                mt_params=mt_params,
                casimir_differentials=casimir_differentials,
                diff_params=diff_params,
            )
        else:
            self.set_from_json_data(json_data)
            self.ffr_curve = SWCurve(
                casimir_differentials=casimir_differentials, 
                g_data=self.g_data,
                diff_params=diff_params,
                mt_params=mt_params,
                z_rotation=self.z_plane_rotation,
                ffr=True,
            )
            self.diff = SWDiff(
                'x',
                g_data=self.g_data,
                diff_params=diff_params,
                mt_params=mt_params,
                z_rotation=self.z_plane_rotation,
            )

        # TODO: SWCurve in a general representation.
        if self.g_data.fundamental_representation_index == 1:
            self.curve = self.ffr_curve
        else:
            logger.warning(
                'Seiberg-Witten curve in a general representation '
                'is not implemented yet.'
            )
            self.curve = None

        # TODO: move these messages somewhere inside set_from_config
        # otherwise the messages about finding ramification points 
        # will show up in the output even before we prin the curve.
        # Display the Seiberg-Witten curve given by the configuration
        logger.info(
            'Seiberg-Witten curve in the 1st fundamental '
            'representation:\n(note: \lambda = x dz)'
            '\n{} = 0\n(numerically\n{}=0\n)'
            .format(sympy.latex(self.ffr_curve.sym_eq),
                    sympy.latex(self.ffr_curve.num_eq))
        )

        logger.info(
            'Seiberg-Witten differential:\n{} dz\n(numerically\n{} dz\n)'
            .format(sympy.latex(self.diff.sym_v),
                    sympy.latex(self.diff.num_v))
        )

        for rp in self.ffr_ramification_points:
            logger.info("{}: z = {}, x = {}, i = {}."
                        .format(rp.label, rp.z, rp.x, rp.i))

        for pct in self.regular_punctures + self.irregular_punctures:
            logger.info('{} at z={}'.format(pct.label, pct.z))

    def set_from_config(self, config, mt_params=None,
                        casimir_differentials=None, diff_params=None,):
        """
        Set attributes by calculating the corresponding 
        values using the configuration.
        """
        # Introduce a clockwise rotation of the z-plane,
        # after the PSL2C transformation, by the following phase.
        # Try rotating by different increments, up to pi/max_pi_div 

        logger = logging.getLogger(self.logger_name)

        min_abs_distance = None
        max_pi_div = 10
        rotate_z_plane = True
        pi_div = 0

        method = config['ramification_point_finding_method']
        logger.info('Ramification point finding method: {}'.format(method))
        if method == 'discriminant':
            get_ramification_points = (
                get_ramification_points_using_discriminant
            )
        elif method == 'system_of_eqs':
            get_ramification_points = (
                get_ramification_points_using_system_of_eqs
            )
        else:
            logger.warning(
                'Unknown or no method set to find ramification points.\n'
                'Use system_of_eqs by default.'
            )
            get_ramification_points = (
                get_ramification_points_using_system_of_eqs
            )

        regular_punctures = get_punctures_from_config(
            config['regular_punctures'], 'Regular puncture',
            diff_params, mt_params,
        )
        irregular_punctures = get_punctures_from_config(
            config['irregular_punctures'], 'Irregular puncture',
            diff_params, mt_params,
        )

        ffr_curve = SWCurve(
            casimir_differentials=casimir_differentials, 
            g_data=self.g_data,
            diff_params=diff_params,
            mt_params=mt_params,
            ffr=True,
        )

        logger.info(
            'Calculating ramification points of '
            'the Seiberg-Witten curve '
            'in the first fundamental rep.'
        )
    
        punctures = regular_punctures + irregular_punctures

        sols = get_ramification_points(
            curve=ffr_curve, 
            diff_params=diff_params,
            mt_params=mt_params,
            accuracy=self.accuracy, 
            punctures=punctures,
            g_data=self.g_data,
            logger_name=self.logger_name,
        )

        ffr_ramification_points = []
        for z_i, (x_j, m_x) in sols:
            rp = RamificationPoint(
                z=PSL2C(mt_params, z_i, numerical=True),
                Ciz=z_i, 
                x=x_j, 
                i=m_x, 
                label=('ramification point #{}'
                       .format(len(ffr_ramification_points)))
            )
            ffr_ramification_points.append(rp)

        # TODO: instead of just printing a warning, increase the accuracy 
        # or other relevant parameters to overcome this issue.
        bpzs = n_remove_duplicate(
            [r.z for r in ffr_ramification_points if r.z != oo], 
            self.accuracy,
        )
        spread = spread_of_branch_points(bpzs, min_spread=MIN_SPREAD)
        if spread > MIN_SPREAD:
            pass
        else:
            logger.info(
                '\n\nWarning!\n'
                'Branch points may be too close to each other '
                'for numerical trivialization. Spread is {}'
                '\nBranch points are:\n{}'
                .format(spread, [rp.z for rp in ffr_ramification_points])
            )

        logger.debug('These are the punctures:')
        for pct in punctures:
            logger.debug('{} at z={}'.format(pct.label, pct.z))

        for pi_div in range(max_pi_div + 1):
            if pi_div == 0:
                # we study the case of no rotations at all.
                z_r = sympy.sympify('1') 
                n_r = -1
                logger.info('The z-plane has not been rotated.')
            elif pi_div == 1:
                # there are no nontrivial rotations when this is 1.
                continue
            else: 
                z_r = sympy.sympify('exp(pi* I / {})'.format(pi_div))
                n_r = 1
                logger.info(
                    'Will try rotating z-plane in increments'
                    ' of pi/{}'.format(pi_div)
                )
            
            z_plane_rotation = z_r

            while (rotate_z_plane is True) and n_r < pi_div:
                if pi_div != 0:
                    logger.info(
                        'The z-plane has been rotated {} times.\n'
                        'Current rotation of the z-plane: {}\n'
                        .format(n_r, z_plane_rotation)
                    )

                # TODO: the code should be able to analyze differentials,
                # and determine where are singularities, instead of giving 
                # them by hand. (However, we'll need to supply monodromy 
                # parameters such as masses, and stokes data for irregular 
                # singularities.)
                
                # Now check if the z-plane needs to be rotated
                # NOTE: if we substitute z' = c z in F(x,z)=0,
                # where c is a phase, the position of punctures 
                # and branch points will rotate contravariantly
                # z_pt -> c^{-1} z_pt

                # z-coords of branch points.
                bpzs = n_remove_duplicate(
                    [r.z / complex(z_plane_rotation)
                     for r in ffr_ramification_points if r.z != oo], 
                    self.accuracy,
                )
                # z-coords of punctures.
                pctzs = n_remove_duplicate(
                    [p.z / complex(z_plane_rotation)
                     for p in punctures if p.z != oo],
                    self.accuracy,
                )
                z_list = bpzs + pctzs
                z_r_list = [z.real for z in z_list]
                if len(z_r_list) > 1:
                    min_x_distance = min([
                        abs(x - y) for i, x in enumerate(z_r_list) 
                        for y in z_r_list[i + 1:]
                    ])
                    min_abs_distance = min([
                        abs(x - y) for i, x in enumerate(z_list)
                        for y in z_list[i + 1:]
                    ])
                    
                elif len(z_r_list) == 1:
                    # No need for the rotation.
                    rotate_z_plane = False
                    break

                elif len(z_r_list) == 0:
                    raise Exception(
                        'Could not find any punctures ' 
                        'or branch points'
                    )

                if min_x_distance > min_abs_distance / len(z_list):
                    logger.info(
                        'All branch points and punctures '
                        'are sufficiently separated horizontally.\n'
                        'Will not rotate z-plane any more.\n'
                    )
                    rotate_z_plane = False
                    break
                
                else:
                    logger.info(
                        'Some branch points or punctures '
                        'are vertically aligned.\n'
                        'Need to rotate the z-plane.\n'
                    )
                    n_r += 1
                    z_plane_rotation *= z_r

            if rotate_z_plane is False:
                break

        if n_r == pi_div == max_pi_div:
            raise ValueError(
                'Could not find a suitable rotation for the z-plane.'
            )

        # Apply the z-rotation to numerical attributes.
        self.z_plane_rotation = z_plane_rotation
        self.ffr_curve = ffr_curve.set_z_rotation(z_plane_rotation)
        self.diff = SWDiff(
            'x',
            g_data=self.g_data,
            diff_params=diff_params,
            mt_params=mt_params,
            z_rotation=self.z_plane_rotation,
        )

        for p in punctures + ffr_ramification_points:
            p.set_z_rotation(1/z_plane_rotation)
        self.regular_punctures = regular_punctures
        self.irregular_punctures = irregular_punctures

        self.ffr_ramification_points = ffr_ramification_points
        self.analyze_ffr_ramification_points()

        # Automatically configure various sizes 
        # if not configured manually.
        if min_abs_distance is None:
            min_abs_distance = DEFAULT_LARGE_STEP_SIZE

        if config['size_of_small_step'] is None:
            config['size_of_small_step'] = min_abs_distance / 100.0

        if config['size_of_large_step'] is None:
            config['size_of_large_step'] = min_abs_distance / 10.0

        if config['size_of_bp_neighborhood'] is None:
            config['size_of_bp_neighborhood'] = min_abs_distance / 2.0

        if config['size_of_puncture_cutoff'] is None:
            config['size_of_puncture_cutoff'] = min_abs_distance / 100.0

    def set_z_rotation(self, z_rotation):
        for p in (
            self.regular_punctures + self.irregular_punctures +
            self.ffr_ramification_points
        ):
            p.set_z_rotation(z_rotation)

        if self.ffr_curve is not None:
            self.ffr_curve.set_z_rotation(1/z_rotation)

        if self.curve is not None:
            self.curve.set_z_rotation(1/z_rotation)

        if self.diff is not None:
            self.diff.set_z_rotation(1/z_rotation)

        self.z_plane_rotation /= z_rotation

    def get_json_data(self):
        json_data = {
            'g_data': self.g_data.get_json_data(),
            'regular_punctures': [p.get_json_data() 
                                  for p in self.regular_punctures],
            'irregular_punctures': [p.get_json_data() 
                                    for p in self.irregular_punctures],
            'ffr_ramification_points': [
                rp.get_json_data() for rp in self.ffr_ramification_points
            ],
            'z_plane_rotation': str(self.z_plane_rotation),
            'accuracy': self.accuracy,
        }

        return json_data

    def set_from_json_data(self, json_data):
        logger = logging.getLogger(self.logger_name)

        self.g_data = GData(json_data=json_data['g_data'],
                            logger_name=self.logger_name,)
        # XXX: Remove the following check after deprecating
        # using older data.
        try:
            self.regular_punctures = [
                Puncture(json_data=data)
                for data in json_data['regular_punctures'] 
            ]
            self.irregular_punctures = [
                Puncture(json_data=data)
                for data in json_data['irregular_punctures'] 
            ]
        except KeyError:
            logger.warning(
                'Loading a JSON data of an older version: '
                'no (ir)regular_punctures data, '
                'use punctures data instead.'
            )
            self.regular_punctures = []
            self.irregular_punctures = [
                Puncture(json_data=data)
                for data in json_data['punctures'] 
            ]

        self.ffr_ramification_points = [
            RamificationPoint(json_data=data)
            for data in json_data['ffr_ramification_points']
        ]
        self.z_plane_rotation = sympy.sympify(json_data['z_plane_rotation'])
        for p in (
            self.regular_punctures + self.irregular_punctures +
            self.ffr_ramification_points
        ):
            p.set_z_rotation(self.z_plane_rotation)
        self.accuracy = json_data['accuracy']

    def get_aligned_xs(self, z_0, near_degenerate_branch_locus=False):
        """
        Returns (aligned_ffr_xs, aligned_xs), where each element is
        a numpy array of x-coordinates of the fibers over z.
        The order of x's is compatible with the order of the weights
        in g_data.weights, in the sense that they are linearly related.
        """
        logger = logging.getLogger(self.logger_name)

        algebra_type = self.g_data.type
        algebra_rank = self.g_data.rank
        fund_rep_index = self.g_data.fundamental_representation_index

        # First order x's of the first fundamental cover
        # according to the order of ''g_data.weights'',
        # then construct the list of x's for the given representation
        # ordered according to weights.
        # NOTE: does this correspond to the 7th irrep for E7?
        #       We need that one for that algebra
        ffr_xs = self.ffr_curve.get_xs(z_0) 

        if algebra_type == 'A':
            # Can consider ffr_xs to be aligned according to ffr_weights,
            # [(1, 0, 0, 0, 0, 0),
            #  (0, 1, 0, 0, 0, 0),
            #  ...,
            #  (0, 0, 0, 0, 0, 1)].
            aligned_ffr_xs = ffr_xs

            if fund_rep_index == 1:
                xs = aligned_ffr_xs
            else:
                xs = self.get_xs_of_weights_from_ffr_xs(aligned_ffr_xs)

        elif algebra_type == 'D':
            # Align ffr_xs according to ffr_weights,
            # [(1, 0, 0, 0, 0),
            #  (0, 1, 0, 0, 0),
            #  ...,
            #  (0, 0, 0, 0, 1),
            #  (-1, 0, 0, 0, 0),
            #  ...
            #  (0, 0, 0, 0, -1)]      
            
            zero_xs = [x for x in ffr_xs if abs(x) <= self.accuracy]
            non_zero_xs = [x for x in ffr_xs if abs(x) > self.accuracy]
            n_zero_xs = len(zero_xs)

            if n_zero_xs == 0:
                sorted_ffr_xs = sorted(
                    ffr_xs, key=lambda x: phase(x),  # reverse=True,
                )

                # Pick x's corresponding to the positive weights.
                # The order among the positive x's is arbitrary.
                positive_xs = sorted_ffr_xs[:algebra_rank]
                # Then pick an x corresponding to each negative weight
                # aligned according to the positive x's.
                unsorted_negative_xs = sorted_ffr_xs[algebra_rank:]
                negative_xs = list(numpy.zeros_like(positive_xs))
            
            # Handle also degenerate cases, such as 
            # when two sheets are identically zero
            else:
                if n_zero_xs != 2 and near_degenerate_branch_locus is False:
                    logger.info(
                        'At z ={} found the following sheets \n{}'.format(
                            z_0, ffr_xs
                        ))
                    raise Exception('Zero sheets must be none or two.')
                else:
                    sorted_ffr_xs = sorted(
                        non_zero_xs, key=lambda x: phase(x),  # reverse=True,
                    )
                    positive_xs = (
                        sorted_ffr_xs[:(algebra_rank - 1)] + [zero_xs[0]]
                    )
                    unsorted_negative_xs = (
                        x for x in ffr_xs if x not in positive_xs
                    )
                    negative_xs = list(numpy.zeros_like(positive_xs))

            for nx in unsorted_negative_xs:
                difference = numpy.fromiter(
                    (abs(px - (-nx)) for px in positive_xs),
                    dtype=float,
                )
                j = difference.argsort()[0]
                # Check the pairing of positive and negative x's.
                px_j = positive_xs[j]
                if numpy.isclose(px_j, -nx) is False:
                    logger.warn(
                        "get_ordered_xs(): No pairing of x's in the D-type, "
                        "({}, {}) != (x, -x).".format(px_j, nx)
                    )
                    logger.info('positive xs : {}'.format(positive_xs))
                if numpy.isclose(
                    px_j, -nx, atol=SHEET_NULL_TOLERANCE
                ) is False:
                    logger.warn(
                        "get_ordered_xs(): No pairing of x's in the D-type,"
                        " ({}, {}) != (x, -x).".format(px_j, nx)
                    )
                    logger.info('positive xs : {}'.format(positive_xs))
                else:
                    # Put the negative x at the same index
                    # as its positive pair.
                    negative_xs[j] = nx
            aligned_ffr_xs = list(
                numpy.concatenate((positive_xs, negative_xs))
            )

            if fund_rep_index == 1:
                xs = aligned_ffr_xs
            else:
                xs = self.get_xs_of_weights_from_ffr_xs(aligned_ffr_xs)

        elif algebra_type == 'E':
            if algebra_rank == 6:
                ffr_weights_list = list(self.g_data.ffr_weights)
                aligned_ffr_xs = align_sheets_for_e_6_ffr(
                    ffr_xs, 
                    ffr_weights_list,
                    near_degenerate_branch_locus=near_degenerate_branch_locus,
                    logger_name=self.logger_name,
                )
                if fund_rep_index == 1:
                    xs = aligned_ffr_xs
                else:
                    # xs = self.get_xs_of_weights_from_ffr_xs(aligned_ffr_xs)
                    raise NotImplementedError
            elif algebra_rank == 7:
                raise NotImplementedError
     
        return (aligned_ffr_xs, xs)

    def get_xs_of_weights_from_ffr_xs(self, ffr_xs):
        g_data = self.g_data
        fund_rep_index = g_data.fundamental_representation_index

        if fund_rep_index == 1:
            return ffr_xs
        else:
            # xs = numpy.zeros(len(g_data.weights), dtype=complex)
            xs = [0.0j for i in range(len(g_data.weights))]

            if g_data.type == 'A' or g_data.type == 'D':
                for i, cs in enumerate(g_data.weight_coefficients):
                    for j, c_j in enumerate(cs):
                        xs[i] += c_j * ffr_xs[j]
            else:
                raise NotImplementedError

            return list(xs)

    def analyze_ffr_ramification_points(self):
        logger = logging.getLogger(self.logger_name)
        rp_type = None
        num_eq = self.ffr_curve.num_eq

        # use dz = z - rp.z & dx = x - rp.x
        dz, dx = sympy.symbols('dz, dx')
        for rp in self.ffr_ramification_points:
            logger.info(
                "Analyzing a ramification point at z = {}, x={}."
                .format(rp.z, rp.x)
            )

            zero_threshold = self.accuracy * 100
            g_data = self.g_data

            if g_data.type == 'A' or g_data.type == 'D':
                local_curve = (
                    num_eq.subs(x, rp.x + dx).subs(z, rp.z + dz)
                    .series(dx, 0, rp.i + 1).removeO()
                    .series(dz, 0, 2).removeO()
                )
            elif g_data.type == 'E' and g_data.rank == 6:
                local_curve = (
                    num_eq.subs(x, rp.x + dx).subs(z, rp.z + dz)
                    .series(dx, 0, rp.i + 1).removeO()
                    .series(dz, 0, 3).removeO()
                )
            else:
                raise NotImplementedError
            logger.debug('\nlocal curve = {}\n'.format(local_curve))
                
            # Classify which type of ramification point
            # type_I: ADE type with x_0 != 0
            #   i.e. F ~ a z + b x^k    (in local coordinates)
            # type_II: D-type with x_0 = 0, but nonedgenerate
            #   i.e. F ~ a z + b x^2r   with r=rank(g)
            # type_III: D-type with x_0 = 0, degenerate
            #   i.e. F ~ x^2 (a z + b x^(2r-2))
            # type_IV: E6-type with x_0 = 0, degenerate
            #   i.e. F ~ x^3 (b x^{24} + a_1 x^{12} z + a_2 z^2)
            # type V: Other case.
            # More cases may be added in the future, in particular 
            # for degenerations of E_6 or E_7 curves.
            
            # XXX: Temporary case for D-type AD theories.
            if (
                g_data.type == 'D' and rp.i == 4 and abs(rp.x) < zero_threshold
            ):
                # No need to grow S-walls from this ramification point,
                # there will be another ramification point that gives
                # the same S-wall. This ramification point is a placeholder.
                rp.ramification_type = 'type_AD'
                rp.sw_diff_coeff = None 
                continue
            elif (
                g_data.type == 'A' or (
                    (g_data.type == 'D' or g_data.type == 'E') and 
                    abs(rp.x) > zero_threshold
                ) or (
                    g_data.type == 'D' and rp.i == 2 
                    and abs(rp.x) < zero_threshold
                ) or (
                    g_data.type == 'E' and (rp.i == 2 or rp.i == 3) 
                    and abs(rp.x) < zero_threshold
                )
            ):
                rp_type = 'type_I'
            elif (
                g_data.type == 'D' and abs(rp.x) < zero_threshold
                and 2 * g_data.rank == rp.i
                and abs(local_curve.n().subs(dx, 0).coeff(dz)) > zero_threshold
            ):
                rp_type = 'type_II'
            elif (
                g_data.type == 'D' and 2 * g_data.rank == rp.i
                and abs(local_curve.n().subs(dx, 0).coeff(dz)) < zero_threshold
            ):
                rp_type = 'type_III'
            elif (
                g_data.type == 'E' and g_data.rank == 6
                and abs(local_curve.n().subs(dx, 0).coeff(dz)) < zero_threshold
            ):
                rp_type = 'type_IV'
            else:
                rp_type = 'type_V'
                logger.info(
                    'Lie algebra {}'.format(g_data.type, g_data.rank)
                )
                logger.info('ramification index {}'.format(rp.i))
                logger.info(
                    'local curve {}'
                    .format(abs(local_curve.n().subs(dx, 0).coeff(dz)))
                )
                raise Exception(
                    'Cannot handle this type of ramification point'.format(
                        local_curve
                    )
                )

            # dx_dz = dx(dz) is the local form of x (the local 
            # coordinate around the ramification point) as a function of z
            # (also intended as a local coordinate near a ramification point)
            if rp_type == 'type_I' or rp_type == 'type_II':
                a = local_curve.n().subs(dx, 0).coeff(dz)
                b = local_curve.n().subs(dz, 0).coeff(dx ** rp.i)
                sw_diff_coeffs_a_b = [complex(a), complex(b)]
                dx_dz = (-1.0 * (a / b) * dz) ** sympy.Rational(1, rp.i)

            elif rp_type == 'type_III':
                a = local_curve.n().coeff(dz).coeff(dx, 2)
                b = local_curve.n().subs(dz, 0).coeff(dx ** rp.i)
                sw_diff_coeffs_a_b = [complex(a), complex(b)]
                dx_dz = (-1.0 * (a / b) * dz) ** sympy.Rational(1, rp.i - 2)

            elif rp_type == 'type_IV':
                a_1 = local_curve.n().coeff(dz).coeff(dx, 15)
                a_2 = local_curve.n().coeff(dz, 2).coeff(dx, 3)
                b = local_curve.n().subs(dz, 0).coeff(dx ** rp.i)
                sw_diff_coeffs_a_b = [complex(a_1), complex(a_2), complex(b)]
                dx_dz_plus = (((
                    -1.0 * a_1 
                    + (a_1 ** 2 - 4.0 * b * a_2) ** sympy.Rational(1, 2)
                ) / (2.0 * b)) * dz) ** sympy.Rational(1, 12)
                dx_dz_minus = (((
                    -1.0 * a_1 
                    - (a_1 ** 2 - 4.0 * b * a_2) ** sympy.Rational(1, 2)
                ) / (2.0 * b)) * dz) ** sympy.Rational(1, 12)

            # FIXME: temporary patch until type_AD is removed
            elif rp_type == 'type_AD':
                a = local_curve.n().subs(dx, 0).coeff(dz)
                b = local_curve.n().subs(dz, 0).coeff(dx ** rp.i)
                sw_diff_coeffs_a_b = [complex(a), complex(b)]
                dx_dz = (-1.0 * (a / b) * dz) ** sympy.Rational(1, rp.i)
            
            logger.debug(
                '\nThe ramification point at (z,x)={} is of {}'.format(
                    [rp.z, rp.x], rp_type
                )
            )

            rp.ramification_type = rp_type    
            rp.sw_diff_coeffs_a_b = sw_diff_coeffs_a_b

            # Now we compute the SW differential actual coefficient.
            # The relation of this to a, b will depend on the 
            # ramification point type. For types I, II, III it should be 
            # sw_diff_coeff = (-a / b)^{1/k} * (self.diff.jac)^{+/-1}
            # for a degree-k ramification point
            # because F ~ a z + b x^k so 
            # \lambda ~ x dz ~ (-a/b)^{1/k} (dz/dz') dz'
            # For type IV there will be two coefficients to consider
            # because (after dropping an overall x^3) 
            # F ~ b x^{24} + a_1 x^{12} z + a_2 z^2
            # so there are two solutions
            # \lambda ~ x dz ~ (c_{+/-})^{1/12} (dz/dz') dz'
            # where 
            # c_{+/-} = (-a_1 +/- \sqrt{(a_1)^2 - 4 b a_2})/(2 b)

            # here num_v is essentially: x * (dz'/dz), where the last factor 
            # is the jacobian from z-plane rotations or PSL2C transformations.
            num_v = self.diff.num_v  

            # now we plug this into num_v, in a neighborhood of x_0
            # we have x = x_0 + dx_dz.
            if (
                rp_type == 'type_I' or rp_type == 'type_II' or 
                rp_type == 'type_III' or rp_type == 'type_AD'
            ):
                local_diff = (
                    num_v.subs(x, rp.x + dx_dz).subs(z, rp.z + dz)
                    .series(dz, 0, 1).removeO()
                )
                # get the coefficient and the exponent of the leading term
                (diff_c, diff_e) = local_diff.leadterm(dz)
                if diff_e == 0:
                    # remove the constant term from the local_diff
                    local_diff -= local_diff.subs(dz, 0)
                    (diff_c, diff_e) = local_diff.leadterm(dz)

                rp.sw_diff_coeff = complex(diff_c.n())

            elif rp_type == 'type_IV':
                local_diff_plus = (
                    num_v.subs(x, rp.x + dx_dz_plus).subs(z, rp.z + dz)
                    .series(dz, 0, 1).removeO()
                )
                local_diff_minus = (
                    num_v.subs(x, rp.x + dx_dz_minus).subs(z, rp.z + dz)
                    .series(dz, 0, 1).removeO()
                )

                # get the coefficient and the exponent of the leading term
                (diff_c_plus, diff_e_plus) = local_diff_plus.leadterm(dz)
                if diff_e_plus == 0:
                    # remove the constant term from the local_diff
                    local_diff_plus -= local_diff_plus.subs(dz, 0)
                    (diff_c_plus, diff_e_plus) = local_diff_plus.leadterm(dz)
                (diff_c_minus, diff_e_minus) = local_diff_minus.leadterm(dz)
                if diff_e_minus == 0:
                    # remove the constant term from the local_diff
                    local_diff_minus -= local_diff_minus.subs(dz, 0)
                    (diff_c_minus, diff_e_minus) = (
                        local_diff_minus.leadterm(dz)
                    )

                rp.sw_diff_coeff = (
                    [complex(diff_c_plus.n()), complex(diff_c_minus.n())]
                )
            else:
                raise NotImplementedError

            


def get_punctures_from_config(
    config_punctures_string, label_prefix,
    diff_params, mt_params,
):
    punctures = []

    if (
        config_punctures_string is None or
        config_punctures_string == 'None' or
        config_punctures_string == ''
    ):
        return punctures

    punctures_str = [
        p_raw_str.strip() for p_raw_str 
        in config_punctures_string.lstrip('[').rstrip(']').split(',')
    ]
    for p_n, p_str in enumerate(punctures_str):
        if len(p_str) == 0:
            continue

        Cipz = sympy.sympify(p_str.strip()).subs(diff_params)
        pz = PSL2C(mt_params, Cipz)

        if pz == oo:
            npz = oo
        else:
            npz = complex(pz.evalf(n=ROOT_FINDING_PRECISION, chop=True))

        punctures.append(
            Puncture(
                z=npz, Ciz=Cipz,
                label=(label_prefix + ' #{}'.format(p_n))
            )
        )
    return punctures


def get_ramification_points_using_system_of_eqs(
    curve=None, 
    diff_params=None, 
    mt_params=None,
    accuracy=None, 
    punctures=None,
    g_data=None,
    logger_name='loom',
):
    logger = logging.getLogger(logger_name)

    sols = []
    f = (curve.sym_eq.subs(diff_params)
         .evalf(n=ROOT_FINDING_PRECISION, chop=True))
    # Make f into the form of f_n/f_d
    f_n, f_d = sympy.cancel(f).as_numer_denom()
    eq_1 = f_n

    # Check the curve if it has the D-type factorization.
    # TODO: checking factorization of an E-type curve using SymPy
    # is too slow. Consider using SAGE.
    if g_data.type != 'E':
        num_factor, f_n_factors = sympy.factor_list(f_n)
        if len(f_n_factors) > 1:
            # TODO: check Casimir differentials too?
            if (
                g_data.type == 'D' and
                len(f_n_factors) == 2 and
                (x, 2) in f_n_factors
            ):
                eq_1 = sympy.simplify(eq_1 / x ** 2)
            else:
                logger.warning('The curve to find ramification points '
                               'has an unknown factorization: {} = {}.'
                               .format(f_n, f_n.factor()))

    eq_2 = eq_1.diff(x)

    # NOTE: solve_poly_system vs. solve
    # sols = sympy.solve_poly_system([f, f.diff(x)], z, x)
    # sols = sympy.solve([f, f.diff(x)], z, x)
    z_x_s = sage_subprocess.solve_system_of_eqs(
        [eq_1, eq_2],
        precision=ROOT_FINDING_PRECISION,
        logger_name=logger_name,
    )

 
    logger.info('Analyze solutions from SAGE.')
    # TODO: Consider calculating the discriminant D(z)
    # and double-check if all the z_i's are found.
    for z_i, x_i in z_x_s:
        #logger.info('Analyze a solution (z, x) = ({}, {})'.format(z_i, x_i))
        # Check if z_i is one of the punctures.
        is_puncture = False
        for p in punctures:
            if abs(z_i - p.Ciz) < accuracy:
                is_puncture = True
        if is_puncture:
            continue

        # Calculate the multiplicity of x_i 
        # by finding the maximum k that satisfies
        # (d_x)^k f_n(x, z_i)|_{x = x_i} = 0.
        f_n_i = f_n.subs(z, z_i)
        m_x = 1
        while (
            abs(f_n_i.diff(x, m_x).subs(x, x_i)
                .evalf(n=ROOT_FINDING_PRECISION)) < accuracy
        ):
            m_x += 1

        if m_x == 1:
            dfdx = f_n_i.diff(x).subs(x, x_i).evalf(n=ROOT_FINDING_PRECISION)
            logger.warning(
                'multiplicity of x is 1 at z = {}, x = {}: '
                'df/dx = {}.'.format(z_i, x_i, dfdx)
            )

        sols.append(
            [complex(mpmath.chop(z_i)), 
             (complex(mpmath.chop(x_i)), m_x)]
        )

    return sols      


def get_ramification_points_using_discriminant(
    curve=None, 
    diff_params=None, 
    mt_params=None,
    accuracy=None, 
    punctures=None,
    g_data=None,
    logger_name='loom',
):
    logger = logging.getLogger(logger_name)
    sols = []    
    f = curve.sym_eq
    # Make f into the form of rf = f_n/f_d
    subs_dict = copy.deepcopy(diff_params)
    rf = sympy.simplify(sympy.cancel(sympy.simplify(f.subs(subs_dict))))
    f_n, f_d = rf.as_numer_denom()
    
    # BEGINNING of bifurcation here: must decide which of 
    # these methods to employ in the long run

    # FIRST METHOD: Use sympy to compute discriminants
    # D_z = sympy.discriminant(f_n, x)

    # if D_z == 0:
    #     logger.info(
    #         'The discriminant of F(x,z) is identically zero. '
    #         'Will work with an effective discriminant.'
    #     )
    #     if g_data.type == 'A':
    #         D_z = sympy.discriminant(f_n / x, x)
    #     if g_data.type == 'D':
    #         D_z = sympy.discriminant(f_n / (x ** 2), x)
    #     # TODO: think through possible generalizations here,
    #     # this is only handling certain special cases with E_6
    #     # i.e. the maximally degenerate branch-point, occurring 
    #     # at the origin of the coulomb branch of pure E_6 SYM
    #     if g_data.type == 'E':
    #         logger.info(
    #             'will work with renormalized curve \n{}'.format(
    #                 f_n / (x ** 3)
    #             )
    #         )
    #         logger.debug(
    #             'after simplification \n{}'.format(
    #                 sympy.simplify(f_n / (x ** 3))
    #             )
    #         )
    #         # FIXME: this is a temporary fix, the computation of the 
    #         # discriminant with sage or sympy is just stuck. Mathematica 
    #         # is able to do this in a fraction of a second though
    #         if (
    #             sympy.expand(f_n / (x ** 3)) == (
    #                 -108 * x ** 24 * z ** 26 - 540 * I * x ** 12 * z ** 15 
    #                 + 540 * I * x ** 12 * z ** 13 - z ** 4 + 2 * z ** 2 - 1
    #             )
    #         ):
    #             D_z = z ** 598 * (z ** 2 - 1) ** 46
    #         else:
    #             D_z = sympy.discriminant(
    #                 sympy.simplify(f_n / (x ** 3)), x
    #             )

    #     logger.debug(
    #         'Will work with the effective discriminant:\n{}'.format(D_z)
    #     )

    # SECOND METHOD: Use SAGE to compute discriminants
    D_z = sage_subprocess.compute_discriminant(
        sympy.expand(f_n)
    )

    if D_z == 0:
        logger.info(
            'The discriminant of F(x,z) is identically zero. '
            'Will work with an effective discriminant.'
        )
        if g_data.type == 'A':
            D_z = sympy.discriminant(f_n / x, x)
        if g_data.type == 'D':
            D_z = sympy.discriminant(f_n / (x ** 2), x)
        # NOTE: think through possible generalizations here,
        # this is only handling certain special cases with E_6
        # i.e. the maximally degenerate branch-point, occurring 
        # at the origin of the coulomb branch of pure E_6 SYM
        if g_data.type == 'E':
            logger.info(
                'will work with renormalized curve \n{}'.format(
                    f_n / (x ** 3)
                )
            )
            logger.debug(
                'after simplification \n{}'.format(
                    sympy.simplify(f_n / (x ** 3))
                )
            )

            # FIXME: this is a temporary fix, the computation of the 
            # discriminant with sage or sympy is just stuck. Mathematica 
            # is able to do this in a fraction of a second though
            if (
                sympy.expand(f_n / (x ** 3)) == (
                    -108 * x**24 * z**26 - 540 * I * x**12 * z**15 
                    + 540 * I * x**12 * z**13 - z**4 + 2 * z**2 - 1
                )
            ):
                D_z = z**598 * (z**2 - 1)**46
            else:
                D_z = sage_subprocess.compute_discriminant(
                    sympy.expand(f_n / (x ** 3))
                )

        logger.debug(
            'Will work with the effective discriminant:\n{}'.format(D_z)
        )

    # END of bifurcation.

    # TODO: perform the cancellation of common factors with SAGE
    # because sympy is not great at factorizing polynomials
    # and sometimes it doesn't simplify entirely,
    # leading to fake branch points where the (unsimplfied) 
    # denominator blows up.
    D_z_n, D_z_d = sympy.cancel(D_z).as_numer_denom()
    factors = sympy.factor_list(sympy.expand(sympy.Poly(D_z_n, z)))[1]
    for fact in factors:
        logger.debug('stuyding roots of factor {}'.format(fact))
        # separate the factor itself and the multiplicity
        f_multiplicity = fact[1]
        f_P = sympy.Poly(fact[0], z)
        # f_m = fact[1]
        cs = [
            c_sym.evalf(
                subs=subs_dict, n=ROOT_FINDING_PRECISION
            ).as_real_imag()
            for c_sym in f_P.all_coeffs()
        ]
        f_P_coeffs = [mpmath.mpc(*c) for c in cs]

        f_roots = None

        fact_eq = f_P.as_expr().evalf(n=ROOT_FINDING_PRECISION, chop=True) 
        f_roots = sage_subprocess.solve_single_eq_single_var(
            fact_eq,
            var='z',
            precision=ROOT_FINDING_PRECISION,
            logger_name=logger_name,
        )

        is_same_z = lambda a, b: abs(a - b) < accuracy    
        
        gathered_f_roots = gather(f_roots, is_same_z) 

        # Find the roots of f(x, z=z_i) for the roots {z_i} of D(z).
        for z_i, zs in gathered_f_roots.iteritems():
            # m_z is the multiplicity of z_i.
            # m_z = len(zs)
            # Check if z_i is one of the punctures.
            is_puncture = False
            for p in punctures:
                if abs(z_i - p.Ciz) < accuracy:
                    is_puncture = True
            if is_puncture:
                continue

            subs_dict[z] = z_i
            f_x_eq = f.subs(subs_dict).evalf(n=ROOT_FINDING_PRECISION)
            f_x_roots = sage_subprocess.solve_single_eq_single_var(
                f_x_eq,
                var='x',
                precision=ROOT_FINDING_PRECISION,
                logger_name=logger_name,
            )

            # In general x-roots have worse errors.
            is_same_x = lambda a, b: abs(a - b) < accuracy / 1e-3
            gathered_f_x_roots = gather(f_x_roots, is_same_x)
            for x_j, xs in gathered_f_x_roots.iteritems():
                # m_x is the multiplicity of x_j.
                # m_x = len(xs) * f_multiplicity
                m_x = len(xs)
                if m_x == 1:
                    continue
                else:
                    sols.append([complex(z_i), (complex(x_j), m_x)])
    return sols


def find_xs_at_z_0(
        sw_data, z_0, x_0=None, num_x=1, ffr=False, use_sage=False
):
    """
    Get x's above z_0 and return the num_x of them 
    which are nearest to x_0.
    NOTE: these are NOT sorted according to the weight/sheet 
    dictionary.
    """
    if ffr is True:
        xs_at_z_0 = sw_data.ffr_curve.get_xs(z_0, use_sage=use_sage) 
    else:
        raise NotImplementedError

    # xs_at_z_0 = sw_data.get_sheets_at_z(z_0, ffr=ffr).values()
    if x_0 is None:
        return xs_at_z_0
    else:
        return sorted(
            xs_at_z_0, lambda x1, x2: cmp(abs(x1 - x_0), abs(x2 - x_0))
        )[:num_x]
   

def null_weight_triples(weights):
    null_vec = numpy.array([0 for i in range(len(list(weights[0])))])
    null_triples = []
    for i, j, k in combinations(range(len(weights)), 3):
        w_i = weights[i]
        w_j = weights[j]
        w_k = weights[k]
        # FIXME: weights are in general arrays of floats, so
        # there may be a numerical issue in the following comparison.
        if (w_i + w_j + w_k == null_vec).all():
            null_triples.append([i, j, k])

    return sorted(null_triples)


def null_sheet_triples(sheets):
    null_triples = []
    delta_x = SHEET_NULL_TOLERANCE
    max_attemps = 100
    n_attempt = 0

    while len(null_triples) != N_NULL_TRIPLES and n_attempt < max_attemps:
        for x_i, x_j, x_k in combinations(sheets, 3):
            if (abs(x_i + x_j + x_k) < delta_x):
                null_triples.append([x_i, x_j, x_k])

        n_attempt += 1
        if len(null_triples) != N_NULL_TRIPLES:
            if len(null_triples) < N_NULL_TRIPLES:
                delta_x = delta_x * 1.2
            elif len(null_triples) > N_NULL_TRIPLES:
                delta_x = delta_x / 2
            null_triples = []

    if len(null_triples) != N_NULL_TRIPLES:
        raise ValueError(
            'Wrong number of overall sheet triples after {} attempts'.format(
                n_attempt
            )
        )

    return sorted(null_triples)


def get_quintets(e, null_triples):
    """ Return the quintets of triples for the given weight/sheet."""
    return [t for t in null_triples if e in t]


def quintet_contained(e, quintet):
    """
    Check whether a weight/sheet e
    is contained in a quintet of triples
    Return 0 for 'No' and '1' for 'Yes'
    """
    ans = 0
    for t in quintet:
        if e in t:
            ans = 1
            break
    return ans


def align_sheets_for_e_6_ffr(
    sheets, 
    weights, 
    near_degenerate_branch_locus=None,
    logger_name='loom',
):
    """
    Return the list of sheets aligned according to the list of weights.
    The output is a list such that sheets[i] will correspond to weights[i].
    """
    logger = logging.getLogger(logger_name)
    if len(sheets) != 27:
        raise Exception('Missing some sheets of the E_6 cover.')

    n_sheets_at_origin = len(
        [x for x in sheets if abs(x) < SHEET_NULL_TOLERANCE]
    )
    
    if n_sheets_at_origin == 27 and near_degenerate_branch_locus is True:
        sorted_sheets = sheets

    elif n_sheets_at_origin == 3:
        # In the degenerate E6 curve (e.g. SYM at origin 
        # of coulomb branch), there are 3 sheets at the origin
        # while other sheets are arranged into two circles
        # in two groups of 12.

        have_same_r = (
            lambda a, b: abs(abs(complex(a)) - abs(complex(b))) 
            < SHEET_NULL_TOLERANCE
        )
        gathered_sheets = gather(sheets, have_same_r)  

        if len(gathered_sheets) != 3:
            logger.info('The following sheets appear: ')
            logger.info(sheets)
            raise ValueError(
                'In a degenerate E_6 curve, sheets should arrange into '
                'three rings in the complex plane. '
                '(One or more rings may shrink to the origin)'
            )

        # radii of the three circles
        r_0, r_1, r_2 = sorted(map(abs, gathered_sheets.keys()))

        # build the three groups of sheets
        g_0 = gathered_sheets.values()[0]
        g_1 = gathered_sheets.values()[1]
        g_2 = gathered_sheets.values()[2]
        
        # normalize the phase to run from 0 to 2 \pi
        norm_phase = lambda w: phase(w) % (2 * pi)
        
        # sort sheets within each ring according to their
        # phase, counter-clockwise starting from the real axis
        g_0_sorted = g_0
        g_1_sorted = sorted(g_1, key=norm_phase)
        g_2_sorted = sorted(g_2, key=norm_phase)

        # Groups of sorted weights, according to the Coxeter 
        # projection.
        # The weights are represented by integer labels, 
        # these are given in the paper on ADE networks and run 
        # from 0 to 26. 
        # Each group of weights is ordered clockwise as they 
        # appear in the Coxeter diagram, starting from the real axis.
        # Note: using the Coxeter diagram from the program cproj
        # gives a shift by 1 in the labels of all weights,
        # but otherwise they coincide precisely with the weights
        # used by loom.
        g_0_weights = [8, 13, 17]
        g_1_weights = [6, 3, 9, 5, 11, 16, 19, 20, 22, 21, 10, 15]
        g_2_weights = [4, 2, 1, 0, 7, 14, 18, 26, 25, 24, 23, 12]

        # Now start sorting the sheets
        sorted_sheets = [None for w in weights]
        # The center ring g_0 is easy
        for i in range(len(g_0_weights)):
            sorted_sheets[g_0_weights[i]] = g_0_sorted[i]
        # In the Coxeter diagram, the phase of the 
        # first weight from g_1 is larger than the phase of 
        # the first weight from g_2 (this one lies precisely 
        # on the positive real axis).
        # Therefore, we have to match sheets to the weights accordingly.
        if norm_phase(g_1_sorted[0]) > norm_phase(g_2_sorted[0]):
            for i in range(len(g_1_weights)):
                sorted_sheets[g_1_weights[i]] = g_1_sorted[i]
            for i in range(len(g_2_weights)):
                sorted_sheets[g_2_weights[i]] = g_2_sorted[i]
        else:
            # In this case we start from the 2nd sheet, not the first one
            # we handle this by shifting cyclically the argument of
            # g_1_sorted
            
            for i in range(len(g_1_weights)):
                sorted_sheets[g_1_weights[i]] = (
                    g_1_sorted[(i + 1) % len(g_1_weights)]
                )
            for i in range(len(g_2_weights)):
                sorted_sheets[g_2_weights[i]] = g_2_sorted[i]

    else:
        n_w_triples = null_weight_triples(weights)
        n_s_triples = null_sheet_triples(sheets)

        # Pick the first sheet and the first weight, 
        # we will declare them to match
        # These must be numbers between 0 and 26,
        # any choice should be equivalent.
        sheet_0_index = 0
        weight_0_index = 0

        x_0 = sheets[sheet_0_index]
        # w_0 = weights[weight_0_index]
        # sorted_sheets[weight_0_index] = x_0

        # The quintet of triples of w_0
        # Note: these are not actual weights, 
        # but rather their integer labels
        q_0 = get_quintets(weight_0_index, n_w_triples)

        # The quintet of SHEET triples of x_0
        # Note: these are the actual values of sheet coordinates
        s_q_0 = get_quintets(x_0, n_s_triples)
        if len(s_q_0) != NULL_TRIPLES_INDIVIDUAL:
            logger.info('the sheets\n{}'.format(sheets))
            logger.info('the choice of x_0\n{}'.format(x_0))
            logger.info('the triples of x_0\n{}'.format(s_q_0))
            raise ValueError(
                'Wrong number of sheet triples: {} instead of {}'.format(
                    len(s_q_0), NULL_TRIPLES_INDIVIDUAL
                )
            )

        # Get the ordered list of sheets appearing in the quintet s_q_0
        known_sheets = [x_0]
        for i in range(NULL_TRIPLES_INDIVIDUAL):
            # Get the (ordered) pair of sheets [x_i, x_j] 
            # from each triple [x_0, x_i, x_j]
            s_pair = [s for s in s_q_0[i] if s != x_0]
            known_sheets.append(s_pair[0])
            known_sheets.append(s_pair[1])

        # Get the ordered list of weights appearing in the quintet q_0
        known_weights = [weight_0_index]
        for t in q_0:
            for i in t:
                if i in known_weights:
                    continue
                else:
                    known_weights.append(i)

        missing_weights = [
            i for i in range(len(weights)) if i not in known_weights
        ]
        missing_sheets = [
            x for x in sheets if x not in known_sheets
        ]

        # The ordering of the known weights is now the following:
        # [i_0, j_1, j_2, k_1, k_2, ...]
        # where i_0 + j_1 + j_2 = 0 = i_0 + k_1 + k_2 = ...
        # (intended as a sum of the actual weights, not their labels)
        # This ordering matters, as we will use it to catalogue the 
        # remaining unknown weights, and we will do the same with the sheets.
        # The overall ordering of the pairs [j_1, j_2], [k_1, k_2], ...
        # is ALMOST free, because of the Weyl symmetry.
        # More precisely, there are 5 pairs and we have a W(D_5)
        # symmetry. So pairs can be permuted and an EVEN number of 
        # 'flips' can be performed. I.e. W(D_5) ~ S_5 x (Z_2)^4
        # This leaves us with two inequiavlent choices:
        # eiher
        # [i_0, j_1, j_2, k_1, k_2, ..., n_1, n_2]
        # or
        # [i_0, j_1, j_2, k_1, k_2, ..., n_2, n_1]
        # We have to try both cases.

        last_pair = known_weights[-2:]
        last_pair_r = [last_pair[1], last_pair[0]]
        known_weights_1 = known_weights
        known_weights_2 = [k_s for k_s in known_weights[:-2]] + last_pair_r

        # List all the combos of which SHEET quintets 
        # must/must not contain all missing SHEETS
        sheet_combos = [
            [
                quintet_contained(x, get_quintets(y, n_s_triples)) 
                for y in known_sheets
            ] for x in missing_sheets
        ]

        for known_weights_i in [known_weights_1, known_weights_2]:
            # Will reorder the sheets according to the weights 
            # they correspond to
            # e.g. if 
            # weights = [w_0, w_1, ...]
            # then we aim for
            # [x_0, x_1, ...] --> [x'_0, x'_1, ...]
            # were on LHS is the list of sheets, and on RHS
            # is the list of sorted_sheets.
            sorted_sheets = [None for w in weights]

            # List all the combos of which WEIGHT quintets 
            # must/must not contain all missing WEIGHTS
            weight_combos = [
                [
                    quintet_contained(j, get_quintets(i, n_w_triples)) 
                    for i in known_weights_i
                ] for j in missing_weights
            ]

            # Now place the known sheets in the corresponding position
            # as dictated by the corresponding known weight.
            for i in range(len(known_weights_i)):
                k_w = known_weights_i[i]
                k_s = known_sheets[i]
                sorted_sheets[k_w] = k_s

            # Now match the patterns of inclusion in the quintets
            # between missing weights and missing sheets.
            # When the patterns match, assign the corresponding 
            # sheet to the list of sorted ones.
            for i in range(len(missing_sheets)):
                s_combo = sheet_combos[i]
                for j in range(len(missing_weights)):
                    w_combo = weight_combos[j]
                    if s_combo == w_combo:
                        sorted_sheets[missing_weights[j]] = missing_sheets[i]
                    else:
                        pass
            
            if None in sorted_sheets:
                # this would mean that sorted_weights_i 
                # is not the correct ordering
                pass
            else:
                break

        if None in sorted_sheets:
            raise ValueError('Cannot match all sheets with weights') 
        elif len(sorted_sheets) != len(delete_duplicates(sorted_sheets)):
            raise ValueError('Duplicate identification of sheets and weights')
            
    return sorted_sheets
