import sympy
import numpy
import logging
import warnings
import copy
import pdb
import sympy.mpmath as mpmath

from sympy import oo
from sympy.mpmath import mp
from sympy.mpmath.libmp.libhyper import NoConvergence
from itertools import combinations
from cmath import phase, pi
from matplotlib import cm as mpl_color_map

import sage_subprocess
from misc import (ctor2, r2toc, PSL2C,
                  delete_duplicates, gather, parse_sym_dict_str,
                  n_remove_duplicate)

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
        self.set_root_color_map()

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

    def weyl_monodromy(self, root, br_loc, direction, reverse=False):
        """
        Returns a new root of the segment of an S-wall
        when it crossed a branch cut from a brancing locus
        which could be a branch point or an irregular singularity.
        """
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

    def set_root_color_map(self):
        """
        Create a mapping between a positive root and a color.
        A negative root will have the same color as the corresponding
        positive root.
        """
        p_roots = self.positive_roots
        num_p_roots = len(p_roots)
        p_root_colors = []
        for i in range(num_p_roots):
            #r, g, b, alpha = mpl_color_map.jet((i / float(n_rts)), bytes=True)
            try:
                color_map = mpl_color_map.viridis
            except AttributeError:
                color_map = mpl_color_map.jet
            r, g, b, alpha = color_map(
                (i / float(num_p_roots)), bytes=True
            )
            p_root_colors.append('#{:02x}{:02x}{:02x}'.format(r, g, b))

        self.root_color_map = p_root_colors

    def get_root_color(self, root):
        logger = logging.getLogger(self.logger_name)
        if self.root_color_map is None:
            self.set_root_color_map()

        for i, p_root in enumerate(self.positive_roots):
            if (numpy.array_equal(p_root, root) or
                numpy.array_equal(-p_root, root)):
                return self.root_color_map[i]

        logger.warning('No color mapped for the root {}'
                       .format(root.tolist()))
        return None


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
            self.is_puncture = is_puncture
        else:
            self.set_from_json_data(json_data)

    def __str__(self):
        return 'z = {}, x = {}, i = {}'.format(self.z, self.x, self.i)

    def __eq__(self, other):
        return self.label == other.label

    def get_json_data(self):
        if self.sw_diff_coeff is None:
            sw_diff_coeff = None
        else:
            sw_diff_coeff = ctor2(self.sw_diff_coeff)
        json_data = {
            'z': ctor2(self.z),
            'Ciz': str(self.Ciz),
            'x': ctor2(self.x),
            'i': ctor2(self.i),
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
        self.i = r2toc(json_data['i'])
        self.label = json_data['label']
        self.ramification_type = json_data['ramification_type']

        if json_data['sw_diff_coeff'] is not None:
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
                 z_rotation=None, ffr=False):
        self.sym_eq = None
        self.num_eq = None

        if ffr is True:
            # Build a cover in the first fundamental representation.
            ffr_eq_str = get_ffr_curve_string(
                casimir_differentials, g_data.type, g_data.rank
            )
            try:
                self.sym_eq = sympy.sympify(ffr_eq_str)
            except:
                raise ValueError('syntax error in the Casimir differentials.')
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

    def get_xs(self, z_0):
        """
        Return a numpy array of x-coordinates over z = z_0.
        """
        if self.num_eq is None:
            raise NotImplementedError

        fx = sympy.simplify(self.num_eq.subs(z, z_0))
        sym_poly = sympy.Poly(fx, x, domain='CC')
        coeff_list = map(complex, sym_poly.all_coeffs())
        return numpy.roots(coeff_list)


class SWDiff:
    def __init__(
            self, v_str, g_data=None, diff_params=None, mt_params=None,
            z_rotation=None,):
        # sym_v is a SymPy expression. 
        self.sym_v = sympy.sympify(v_str)
        # num_v is from sym_v with its parameters 
        # substituted with numerical values.
        # NOTE: We apply PSL2C only to the numerical curve
        # for the simplicity of analysis.
        Ciz = PSL2C(mt_params, z_rotation * z, inverse=True) 
        dCiz = Ciz.diff(z)
        self.num_v = (
            (self.sym_v.subs(z, Ciz) * dCiz).subs(diff_params)
            .evalf(n=ROOT_FINDING_PRECISION, chop=True)
        )


class SWDataBase(object):
    """
    A class containing the geometric data of a Seiberg-Witten curve
    in the first fundamental representation,
        x^N + \sum_{k=2}^N a_k(z) dz^k x^{N-k}, 
    where x is the Seiberg-Witten differential of the form 
        x = x dz.

    This class also includes ramification points of the curve, and 
    another curve in the representation given in the configuration.

    This is the base class of SWDataWithTrivialization, where 
    the trivialization data of the curve is contained.
    """
    def __init__(self, config, logger_name='loom', json_data=None,):
        self.logger_name = logger_name
        logger = logging.getLogger(self.logger_name)

        self.g_data = None
        self.regular_punctures = []
        self.irregular_punctures = []
        self.ffr_ramification_points = None
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

        logger.info(
            'Seiberg-Witten curve in the 1st fundamental '
            'representation:\n{} = 0\n(numerically\n{}=0\n)'
            .format(sympy.latex(self.ffr_curve.sym_eq),
                    sympy.latex(self.ffr_curve.num_eq))
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

        self.diff = SWDiff(
            'x',
            g_data=self.g_data,
            diff_params=diff_params,
            mt_params=mt_params,
            z_rotation=self.z_plane_rotation,
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
        self.accuracy = json_data['accuracy']

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
                
                regular_punctures = get_punctures_from_config(
                    config['regular_punctures'], 'Regular puncture',
                    diff_params, mt_params, z_plane_rotation,
                )
                irregular_punctures = get_punctures_from_config(
                    config['irregular_punctures'], 'Irregular puncture',
                    diff_params, mt_params, z_plane_rotation,
                )

                ffr_curve = SWCurve(
                    casimir_differentials=casimir_differentials, 
                    g_data=self.g_data,
                    diff_params=diff_params,
                    mt_params=mt_params,
                    z_rotation=z_plane_rotation,
                    ffr=True,
                )

                logger.info(
                    'Calculating ramification points of '
                    'the Seiberg-Witten curve '
                    'in the first fundamental rep.'
                )
            
                punctures = regular_punctures + irregular_punctures

                ffr_ramification_points = []
                # print '\nthis is the curve'
                # print ffr_curve.num_eq
                # print '\nthese are the roots'
                # print self.g_data.roots.tolist()
                sols = get_ramification_points(
                    curve=ffr_curve, 
                    diff_params=diff_params,
                    mt_params=mt_params,
                    accuracy=self.accuracy, 
                    punctures=punctures,
                    g_data=self.g_data,
                    logger_name=self.logger_name,
                )
                
                for z_i, (x_j, m_x) in sols:
                    rp = RamificationPoint(
                        # Note: if we substitute z' = c z in F(x,z)=0,
                        # where c is a phase, the position of punctures 
                        # and branch points will rotate contravariantly
                        # z_pt -> c^{-1} z_pt
                        z=(PSL2C(mt_params, z_i, numerical=True) /
                           complex(z_plane_rotation)),
                        Ciz=z_i, 
                        x=x_j, 
                        i=m_x, 
                        label=('ramification point #{}'
                               .format(len(ffr_ramification_points)))
                    )
                    ffr_ramification_points.append(rp)


                logger.debug('These are the punctures:')
                for pct in punctures:
                    logger.debug('{} at z={}'.format(pct.label, pct.z))

                # Now check if the z-plane needs to be rotated

                # z-coords of branch points.
                bpzs = n_remove_duplicate(
                    [r.z for r in ffr_ramification_points if r.z != oo],
                    self.accuracy,
                )
                # z-coords of punctures.
                pctzs = n_remove_duplicate(
                    [p.z for p in punctures if p.z != oo],
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
                    # This will block rotation of the z-plane
                    # enable only for debugging purposes
                    # rotate_z_plane = False
                    #
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

        self.z_plane_rotation = z_plane_rotation
        self.regular_punctures = regular_punctures
        self.irregular_punctures = irregular_punctures
        self.ffr_ramification_points = ffr_ramification_points
        self.ffr_curve = ffr_curve

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
                aligned_ffr_xs = sort_sheets_for_e_6_ffr(
                    ffr_xs, ffr_weights_list
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


def get_punctures_from_config(
    config_punctures_string, label_prefix,
    diff_params, mt_params, z_plane_rotation,
):
    punctures = []

    if config_punctures_string is None:
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
            # Note: if we substitute z' = c z in F(x,z)=0,
            # where c is a phase, the position of punctures 
            # and branch points will rotate contravariantly
            # z_pt -> c^{-1} z_pt
            npz = complex(
                (pz / z_plane_rotation)
                .evalf(n=ROOT_FINDING_PRECISION, chop=True)
            )
        punctures.append(
            Puncture(
                z=npz, Ciz=Cipz,
                label=(label_prefix + ' #{}'.format(p_n))
            )
        )
    return punctures

# E_6 curve strings
#
# tau_str = 'z + 1/z + ({u_6})'
# q_1_str = (
#     '270*x^(15) + 342*(({u_1}))*x^(13) + 162*(({u_1}))^2*x^(11)'  
#     '- 252*(({u_2}))*x^(10) + (26*(({u_1}))^3 + 18*(({u_3})))*x^9' 
#     '- 162*(({u_1}))*(({u_2}))*x^8 + (6*(({u_1}))*(({u_3})) '
#     '- 27*(({u_4})))*x^7' 
#     '- (30*(({u_1}))^2*(({u_2})) - 36*(({u_5})))*x^6' 
#     '+ (27*(({u_2}))^2 - 9*(({u_1}))*(({u_4})))*x^5' 
#     '- (3*(({u_2}))*(({u_3})) - 6*(({u_1}))*(({u_5})))*x^4' 
#     '- 3*(({u_1}))*(({u_2}))^2*x^3 - 3*(({u_2}))*(({u_5}))*x '
#     '- (({u_2}))^3'
# )
# q_2_str = '1/(2*x^3)*(({q_1})^2 - ({p_1})^2*({p_2}))'
# p_1_str = (
#     '78*x^10 + 60*(({u_1}))*x^8 + 14*(({u_1}))^2*x^6 '
#     '- 33*(({u_2}))*x^5' 
#     '+ 2*(({u_3}))*x^4 - 5*(({u_1}))*(({u_2}))*x^3 - (({u_4}))*x^2 '
#     '- (({u_5}))*x - (({u_2}))^2'
# )
# p_2_str = (
#     '12*x^10 + 12*(({u_1}))*x^8 + 4*(({u_1}))^2*x^6 '
#     '- 12*(({u_2}))*x^5 + (({u_3}))*x^4' 
#     '- 4*(({u_1}))*(({u_2}))*x^3 - 2*(({u_4}))*x^2 + 4*(({u_5}))*x '
#     '+ (({u_2}))^2'
# )


# Another version of the SW curve, obtained after replacing 
# x -> - z x / 2
# so that now x ~ \lambda_{SW}
#
phi_12_str = '({u_6})'
q_1_str = (
    '-(13/256) * x^9 * ({u_1})^3 * z^9 '
    '- 15/32 * x^6 * ({u_1})^2 * ({u_2}) * z^6'
    '-(81 * x^(11) * ({u_1})^2 * z^(11))/(1024) '
    '+ (3/8) * x^3 * ({u_1}) * ({u_2})^2 * z^3 '
    '- (81/128) * x^8 * ({u_1}) * ({u_2}) * z^8'
    '- (3/64) * x^7 * ({u_1}) * ({u_3}) * z^7 '
    '+ (9 / 32) * x^5 * ({u_1}) * ({u_4}) * z^5'
    '+ (3/8) * x^4 * ({u_1}) * ({u_5}) * z^4 '
    '- (171 * x^(13) * ({u_1}) * z^(13))/(4096)'
    '- ({u_2})^3 - (27/32) * x^5 * ({u_2})^2 * z^5 '
    '- (3/16) * x^4 * ({u_2}) * ({u_3}) * z^4'
    '+ (3/2) * x * ({u_2}) * ({u_5}) * z '
    '- (63/256) * x^(10) * ({u_2}) * z^(10) '
    '- (9/256) * x^9 * ({u_3}) * z^9 '
    '+ (27/128) * x^7 * ({u_4}) * z^7 '
    '+ (9/16) * x^6 * ({u_5}) * z^6 '
    '- (135 * x^(15) * z^(15))/(16384)'
)
q_2_str = '1/2 * ((-2) / (x * z))^3 *(({q_1})^2 - ({p_1})^2*({p_2}))'
p_1_str = (
    '(7/32) * x^6 * ({u_1})^2 * z^6 + (5/8) * x^3 * ({u_1}) * ({u_2}) * z^3' 
    '+ (15/64) * x^8 * ({u_1}) * z^8 - ({u_2})^2 '
    '+ (33/32) * x^5 * ({u_2}) * z^5 + (1/8) * x^4 * ({u_3}) * z^4 '
    '- (1/4) * x^2 * ({u_4}) * z^2 + (x * ({u_5}) * z)/2 '
    '+ (39 * x^(10) * z^(10))/512'
)
p_2_str = (
    '(1/256) * x * z * (16 * x^5 * ({u_1})^2 * z^5 '
    '+ 12 * x^7 * ({u_1}) * z^7 + 16 * x^3 * ({u_3}) * z^3 '
    '- 128 * x * ({u_4}) * z - 512 * ({u_5}) + 3 * x^9 * z^9) '
    '+ (1/8) * ({u_2}) * (4 * x^3 * ({u_1}) * z^3 + 3 * x^5 * z^5) '
    '+ ({u_2})^2'
)


def get_ffr_curve_string(casimir_differentials, g_type, g_rank):
    """
    Construct a Seiberg-Witten curve in the 1st fundamental representation
    using Casimir differentials.
    """
    cs = []
    # FIXME: make checks for the casimirs given by the user
    # For example, I was able to insert the 0th Casimir with no complaint
    # Moreover I didn get the curve I wanted because below we add +1 to phi_0
    # We should prevent the wrong casimirs from being inserted at all.
    if g_type == 'A':
        N = g_rank + 1
        for k, phi_k in casimir_differentials.iteritems():
            cs.append([k, phi_k])

    elif g_type == 'D':
        N = 2 * g_rank
        for k, phi_k in casimir_differentials.iteritems():
            if k % 2 != 0:
                raise ValueError(
                    'phi_{} = {} dz^k'.format(k, phi_k) + 
                    'is an incorrect casimir differential '
                    'for a D-type curve.' 
                )
            else:
                cs.append([k, phi_k])

    elif g_type == 'E':
        phi = casimir_differentials
        # u_(1, 2, 3, 4, 5) = phi[2, 5, 6, 8, 9, 12]
        if g_rank == 6:
            # The following goes with the other presentation of the SW curve.
            # 
            phi_12 = phi_12_str.format(u_6=phi[12])
            q_1 = q_1_str.format(u_1=phi[2], u_2=phi[5], u_3=phi[6],
                                u_4=phi[8], u_5=phi[9])
            p_1 = p_1_str.format(u_1=phi[2], u_2=phi[5], u_3=phi[6],
                                u_4=phi[8], u_5=phi[9])
            p_2 = p_2_str.format(u_1=phi[2], u_2=phi[5], u_3=phi[6],
                                u_4=phi[8], u_5=phi[9])
            q_2 = q_2_str.format(q_1=q_1, p_1=p_1, p_2=p_2)
            curve_str = (
               '(1/2)*(-(1/2) *z *x)^3*({phi_12})^2 '
               '- ({q_1})*({phi_12}) + ({q_2})'
               .format(phi_12=phi_12, q_1=q_1, q_2=q_2)
            )
            
            # tau = tau_str.format(u_6=phi[12])
            # q_1 = q_1_str.format(u_1=phi[2], u_2=phi[5], u_3=phi[6],
            #                     u_4=phi[8], u_5=phi[9])
            # p_1 = p_1_str.format(u_1=phi[2], u_2=phi[5], u_3=phi[6],
            #                     u_4=phi[8], u_5=phi[9])
            # p_2 = p_2_str.format(u_1=phi[2], u_2=phi[5], u_3=phi[6],
            #                     u_4=phi[8], u_5=phi[9])
            # q_2 = q_2_str.format(q_1=q_1, p_1=p_1, p_2=p_2)
            # curve_str = (
            #    '(1/2)*(x^3)*({tau})^2 - ({q_1})*({tau}) + ({q_2})'
            #    .format(tau=tau, q_1=q_1, q_2=q_2)
            # )

            # print '\nthe curve string'
            # print curve_str

            return curve_str

    else:
        raise NotImplemented(
            'get_ffr_curve_string(): construction of a Seiberg-Witten curve '
            'of g = {}_{} is not implemented.'.format(g_type, g_rank)
        )

    # g_type == 'A' or g_type == 'D':
    curve_str = 'x^{} '.format(N)
    for k, c_k in cs:
        curve_str += ' + ({}) '.format(c_k)
        if k != N:
            curve_str += ' * x^{}'.format(N - k)

    return curve_str


#def get_ramification_points(
#    curve=None, 
#    diff_params=None, 
#    mt_params=None,
#    z_rotation=None,
#    accuracy=None, 
#    punctures=None,
#    method=None,
#    g_data=None,
#    logger_name='loom',
#):
#    logger = logging.getLogger(logger_name)
#
#    if curve.sym_eq is None:
#        raise NotImplementedError
#
#    if method == 'discriminant':
#        sols = get_ramification_points_using_discriminant(
#            curve=curve, 
#            diff_params=diff_params, 
#            mt_params=mt_params,
#            accuracy=accuracy, 
#            punctures=punctures,
#            g_data=g_data,
#            logger_name=logger_name,
#        )
#    else:
#        if method != 'system_of_eqs':
#            logger.warning(
#                'Unknown or no method set to find ramification points.\n'
#                'Use system_of_eqs by default.'
#            )
#        sols = get_ramification_points_using_system_of_eqs(
#            curve=curve, 
#            diff_params=diff_params, 
#            mt_params=mt_params,
#            accuracy=accuracy, 
#            punctures=punctures,
#            logger_name=logger_name,
#        )
#
#    ramification_points = []
#
#    for z_i, (x_j, m_x) in sols:
#        rp = RamificationPoint(
#            # Note: if we substitute z' = c z in F(x,z)=0,
#            # where c is a phase, the position of punctures 
#            # and branch points will rotate contravariantly
#            # z_pt -> c^{-1} z_pt
#            z=PSL2C(mt_params, z_i, numerical=True) / complex(z_rotation),
#            Ciz=z_i, 
#            x=x_j, 
#            i=m_x, 
#            label=('ramification point #{}'
#                   .format(len(ramification_points)))
#        )
#        ramification_points.append(rp)
#
#    return ramification_points


def get_ramification_points_using_system_of_eqs(
    curve=None, 
    diff_params=None, 
    mt_params=None,
    accuracy=None, 
    punctures=None,
    g_data=None,
    logger_name='loom',
):
    sols = []
    f = curve.sym_eq
    # Make f into the form of f_n/f_d
    f_n, f_d = sympy.cancel(f).as_numer_denom()
    eq_1 = f_n.subs(diff_params).evalf(n=ROOT_FINDING_PRECISION, chop=True)

    d_x_f_n, d_x_f_d = sympy.cancel(f.diff(x)).as_numer_denom()
    eq_2 = d_x_f_n.subs(diff_params).evalf(n=ROOT_FINDING_PRECISION, chop=True) 

    # NOTE: solve_poly_system vs. solve
    # sols = sympy.solve_poly_system([f, f.diff(x)], z, x)
    # sols = sympy.solve([f, f.diff(x)], z, x)
    z_x_s = sage_subprocess.solve_system_of_eqs(
        [eq_1, eq_2],
        precision=ROOT_FINDING_PRECISION,
        logger_name=logger_name,
    )
    # TODO: Consider calculating the discriminant D(z)
    # and double-check if all the z_i's are found.
    for z_i, x_i in z_x_s:
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
        f_n_i = f_n.subs(diff_params).subs(z, z_i)
        m_x = 1
        while (
            abs(f_n_i.diff(x, m_x).subs(x, x_i)
                .evalf(n=ROOT_FINDING_PRECISION)) < accuracy
        ):
            m_x += 1

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
    
    # Old way -- keep until testing is complete:
    # f = curve.sym_eq
    # Make f into the form of rf = f_n/f_d
    # rf = sympy.cancel(f)
    # rf = sympy.cancel(f.subs(subs_dict))
    # f_n, f_d = rf.as_numer_denom()
    # subs_dict = copy.deepcopy(diff_params)
    
    # New way:
    f = curve.sym_eq
    # Make f into the form of rf = f_n/f_d
    subs_dict = copy.deepcopy(diff_params)
    rf = sympy.simplify(sympy.cancel(sympy.simplify(f.subs(subs_dict))))
    f_n, f_d = rf.as_numer_denom()
    
    # # Find the roots of D(z), the discriminant of f(x, z)
    # # as a polynomial of x. 
    # # TODO: test if sage's discriminant works well also with
    # # A and D-types curves. Then evaluate whether to get rid of
    # # the sympy discriminant for them.
    # if g_data.type == 'E':
    #     D_z = sage_subprocess.compute_discriminant(f_n.subs(subs_dict))
    # else:
    #     D_z = sympy.discriminant(f_n.subs(subs_dict), x)

    D_z = sympy.discriminant(f_n.subs(subs_dict), x)

    if D_z == 0:
        logger.info(
            'The discriminant of F(x,z) is identically zero. '
            'Will work with an effective discriminant.'
        )
        if g_data.type == 'A':
            D_z = sympy.discriminant(f_n.subs(subs_dict) / x, x)
        if g_data.type == 'D':
            D_z = sympy.discriminant(f_n.subs(subs_dict) / (x ** 2), x)
        # NOTE: think through possible generalizations here,
        # this is only handling certain special cases with E_6
        # i.e. the maximally degenerate branch-point, occurring 
        # at the origin of the coulomb branch of pure E_6 SYM
        if g_data.type == 'E':
            logger.info(
                'will work with renormalized curve {}'.format(
                    f_n.subs(subs_dict) / (x ** 3)
                )
            )
            logger.debug(
                'after simplification {}'.format(
                    sympy.simplify(f_n.subs(subs_dict) / (x ** 3))
                )
            )
            D_z = sympy.discriminant(
                sympy.simplify(f_n.subs(subs_dict) / (x ** 3)), x
            )
        logger.debug(
            'Will work with the effective discriminant:\n{}'.format(D_z)
        )

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

        # OLD METHOD
        #
        # Increase maxsteps & extraprec when root-finding fails.
        # polyroots_maxsteps = ROOT_FINDING_MAX_STEPS
        # polyroots_extra_precision = ROOT_FINDING_PRECISION
        # while f_roots is None:
        #     try:
        #         f_roots = mpmath.polyroots(
        #             f_P_coeffs, 
        #             maxsteps=polyroots_maxsteps,
        #             extraprec=polyroots_extra_precision,
        #         )
        #     except NoConvergence:
        #         logger.warning(
        #             'mpmath.polyroots failed; increase maxsteps & extraprec '
        #             'by 10.'
        #         )
        #         polyroots_maxsteps += 10
        #         polyroots_extra_precision += 10
                
        fact_eq = f_P.as_expr().evalf(n=ROOT_FINDING_PRECISION, chop=True) 
        f_roots = sage_subprocess.solve_single_eq(
            [fact_eq],
            precision=ROOT_FINDING_PRECISION,
            logger_name=logger_name,
        )

        # # TODO: numerical errors in the determination of 
        # # roots of the discriminant induces an artificial 
        # # separation in the roots sometimes.
        # # This is especially the case with E_6, hence the following
        # # criterion. Find something better.
        # if g_data.type == 'E':
        #     is_same_z = lambda a, b: abs(a - b) < 0.001
        # else:
        #     is_same_z = lambda a, b: abs(a - b) < accuracy   

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
            
            # OLD METHOD
            #      
            # subs_dict[z] = z_i
            # f_x_cs = [c.evalf(subs=subs_dict, n=ROOT_FINDING_PRECISION) 
            #           for c in sympy.Poly(f, x).all_coeffs()]
            # f_x_coeffs = [mpmath.mpc(*c.as_real_imag()) for c in f_x_cs]
            # f_x_roots = None
            #
            # while f_x_roots is None:
            #     try:
            #         f_x_roots = mpmath.polyroots(
            #             f_x_coeffs,
            #             maxsteps=polyroots_maxsteps,
            #             extraprec=polyroots_extra_precision
            #         )
            #     except NoConvergence:
            #         logger.warning(
            #             'mpmath.polyroots failed; increase maxsteps & '
            #             'extraprec by 10.'
            #         )
            #         polyroots_maxsteps += 10
            #         polyroots_extra_precision += 10

            subs_dict[z] = z_i
            f_x_eq = f.subs(subs_dict).evalf(n=ROOT_FINDING_PRECISION)
            # print '\n this is the f_x equation: {}'.format(f_x_eq)
            f_x_roots = sage_subprocess.solve_single_eq_x(
                [f_x_eq],
                precision=ROOT_FINDING_PRECISION,
                logger_name=logger_name,
            )
            # print 'this is the equation for sage : {}'.format(f_x_eq)
            # print 'this is the solution {}'.format(f_x_roots)

            # In general x-roots have worse errors.
            is_same_x = lambda a, b: abs(a - b) < accuracy / 1e-2
            gathered_f_x_roots = gather(f_x_roots, is_same_x)

            # print '\nfound potential ramification point'
            # print [complex(z_i), gathered_f_x_roots]

            for x_j, xs in gathered_f_x_roots.iteritems():
                # m_x is the multiplicity of x_j.
                m_x = len(xs) * f_multiplicity
                if m_x == 1:
                    continue

                sols.append([complex(z_i), (complex(x_j), m_x)])

    return sols


def find_xs_at_z_0(sw_data, z_0, x_0=None, num_x=1, ffr=False):
    """
    Get x's above z_0 and return the num_x of them 
    which are nearest to x_0.
    NOTE: these are NOT sorted according to the weight/sheet 
    dictionary.
    """
    if ffr is True:
        xs_at_z_0 = sw_data.ffr_curve.get_xs(z_0) 
    else:
        raise NotImplementedError

    # xs_at_z_0 = sw_data.get_sheets_at_z(z_0, ffr=ffr).values()
    if x_0 is None:
        return xs_at_z_0
    else:
        return sorted(xs_at_z_0,
                      lambda x1, x2: cmp(abs(x1 - x_0), abs(x2 - x_0)))[:num_x]
   

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
            print '\n{} triples'.format(len(null_triples))
            print 'delta_x = {}'.format(delta_x)
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


def sort_sheets_for_e_6_ffr(sheets, weights):
    """
    Return the list of sheets sorted according to the list of weights.
    The output is a list such that sheets[i] will correspond to weights[i].
    """

    n_sheets_at_origin = len(
        [x for x in sheets if abs(x) < SHEET_NULL_TOLERANCE]
    )
    
    if n_sheets_at_origin >= 3:
        # print "There are {} sheets at the origin. ".format(n_sheets_at_origin)
        have_same_r = (
            lambda a, b: abs(abs(complex(a)) - abs(complex(b))) 
            < SHEET_NULL_TOLERANCE
        )
        gathered_sheets = gather(sheets, have_same_r)  
        # print 'the sheets'
        # print sheets
        # print 'gathered sheets {}'.format(gathered_sheets)
        # print len(gathered_sheets)
        # print 'radii'
        # print map(abs, gathered_sheets.keys())

        if len(gathered_sheets)!=3:
            raise ValueError(
                'In a degenerate E_6 curve, sheets should arrange into '
                'three rings in the complex plane. '
                '(One ring may shrink to the origin)'
            )

        # radii of the three circles
        r_0, r_1, r_2 = map(abs, gathered_sheets.keys())

        # build groups of sheets
        g_0 = [x for x in sheets if abs(abs(x) - r_0) < SHEET_NULL_TOLERANCE]
        g_1 = [x for x in sheets if abs(abs(x) - r_1) < SHEET_NULL_TOLERANCE]
        g_2 = [x for x in sheets if abs(abs(x) - r_2) < SHEET_NULL_TOLERANCE]

        # print 'groups of sheets'
        # print g_0
        # print g_1
        # print g_2
        
        # normalize the phase to run from 0 to 2 \pi
        norm_phase = lambda w: phase(w) % (2 * pi)
        # sort sheets within each ring according to their
        # phase, counter-clockwise starting from the real axis
        g_0_sorted = g_0
        g_1_sorted = sorted(g_1, key=norm_phase)
        g_2_sorted = sorted(g_2, key=norm_phase)

        # print 'sorted groups of sheets'
        # print g_0_sorted
        # print g_1_sorted
        # print g_2_sorted

        # Groups of sorted weights, according to the Coxeter 
        # projection.
        # The weights are represented by integer labels, 
        # these are given in the paper on ADE networks and run 
        # from 0 to 26. 
        # Each group of weights is ordered clockwise as they 
        # appear in the Coxeter diuagram, starting from the real axis.
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
        if norm_phase(g_1[0]) > norm_phase(g_2[0]):
            for i in range(len(g_1_weights)):
                sorted_sheets[g_1_weights[i]] = g_1_sorted[i]
            for i in range(len(g_2_weights)):
                sorted_sheets[g_2_weights[i]] = g_2_sorted[i]
        else:
            # In this case we start from the 2nd sheet, not the first one
            # we handle this by shifting cyclically the argument of
            # g_1_sorted
            for i in range(len(g_1_weights)):
                sorted_sheets[g_1_weights[i]] = g_1_sorted[(i + 1) % len(g_1_weights)]
            for i in range(len(g_2_weights)):
                sorted_sheets[g_2_weights[i]] = g_2_sorted[i]


        # raise Exception("This is a degenerate case, cannot handle it yet!")
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
            print '\nthe sheets'
            print sheets
            print '\nthe choice of x_0'
            print x_0
            print '\nthe triples of x_0'
            print s_q_0
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

