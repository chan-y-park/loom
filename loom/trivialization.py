import numpy
import logging
import sympy

from cmath import exp, pi, phase
from sympy import oo
from numpy.linalg import matrix_rank
from itertools import combinations
from pprint import pprint
from heapq import nsmallest

from geometry import SWDataBase
from misc import (
    n_remove_duplicate, ctor2, r2toc, delete_duplicates, is_weyl_monodromy
)

x, z = sympy.symbols('x z')


# TODO: set both of the following automatically: e.g. using the
# minimal_distance attribute of the SW fibration

# number of steps used to track the sheets along a leg
# the path used to trivialize the cover at any given point
N_PATH_TO_PT = 100

# number of steps for each SEGMENT of the path around a
# branching point (either branch-point, or irregular singularity)
# N_PATH_AROUND_PT = 60
N_PATH_AROUND_PT = 100

# Number of times the tracking of sheets is allowed to automatically zoom in.
# Usual values
MAX_ZOOM_LEVEL = 3
ZOOM_FACTOR = 10
# Tuned for E_6
# MAX_ZOOM_LEVEL = 0
# ZOOM_FACTOR = 10

# Tolerance for recognizing colliding sheets at a branch-point
BP_PROXIMITY_THRESHOLD = 0.05


class BranchPoint:
    """
    The BranchPoint class.

    This class is strictly related to the
    cover corresponding to the first fundamental
    representation.

    Attributes
    ----------

    z :
        The position of the branch point on the z-plane

    trivialization :
        The trivialization of the cover to which the
        branch point is associated.

    groups :
        A list of groups of sheets which collide together
        at the branch point.

    singles :
        The list of sheets which do not collide with any
        other sheet.

    enum_sh :
        The enumerated sheets at the branch point.
        A list of pairs [i, x] where i is the sheet
        identifier referring to the reference sheets
        of the trivialization class; x is the corresponding
        coordinate in the fiber above the branch point.

    path_to_bp - UNAVAILABLE :
        A path running from the basepoint of the trivialization
        to the branch point without crossing any branch cut.

    sheet_tracks_to_bp - UNAVAILABLE :
        A list of sheet tracks, i.e. the x-values of each
        sheet as it is tracked along a path that runs to
        the branch point, to determine collision structure
        of the various sheets.

    positive_roots :
        A minimal list of positive roots characterizing the
        groups of colliding sheets at the branch point.

    path_around_bp - UNAVAILABLE :
        A path encircling the branch point and no one else,
        used to compute the monodromy.

    sheet_tracks_around_bp - UNAVAILABLE :
        A list of sheet tracks, i.e. the x-values of each
        sheet as it is tracked along a path that runs around
        the branch point, to determine the monodromy.

    monodromy :
        The monodromy matrix acting on the column vector
        of sheets (hence, acting FROM the left).
        Sheets are ordered according to the reference
        sheets of the trivialization.

    order :
        At a branch point, the dual of the higgs field
        lies on the boundary of a Weyl chamber.
        In general, it will li at the intersection of
        k of the walls delimiting the chamber.
        The order of the branch point is then k + 1.

    ffr_ramification_points :
        A list of all ramification point objects, which
        lie in the fiber above the branch point.

    """
    def __init__(
        self, z=None, json_data=None, ffr_ramification_points=None,
        logger_name='loom',
    ):
        self.logger_name = logger_name
        if json_data is None:
            self.z = z
            self.label = None
            self.monodromy = None

            self.groups = None
            self.positive_roots = None
            self.order = None

            self.ffr_ramification_points = None
            self.seeds = []
        else:
            self.set_from_json_data(json_data, ffr_ramification_points,)

        self.data_attributes = [
            'z', 'label', 'monodromy', 'groups', 'positive_roots',
            'order', 'ffr_ramification_points'
        ]

    def set_z_rotation(self, z_rotation):
        if self.z != oo:
            self.z /= complex(z_rotation)

    def get_json_data(self):
        json_data = {
            'z': ctor2(self.z),
            'label': self.label,
            'monodromy': self.monodromy.tolist(),
            'groups': self.groups,
            'positive_roots': self.positive_roots.tolist(),
            'order': self.order,
            'ffr_ramification_points': [
                rp.label for rp in self.ffr_ramification_points
            ],
            # 'seeds': map(ctor2, self.seeds),
        }
        return json_data

    def set_from_json_data(self, json_data, ffr_ramification_points,):
        self.z = r2toc(json_data['z'])
        self.label = json_data['label']
        self.monodromy = numpy.array(json_data['monodromy'])
        self.groups = json_data['groups']
        self.positive_roots = numpy.array(json_data['positive_roots'])
        self.order = json_data['order']
        self.ffr_ramification_points = []
        # self.seeds = map(r2toc, json_data['seeds'])
        for rp_label in json_data['ffr_ramification_points']:
            for rp in ffr_ramification_points:
                if rp_label == rp.label:
                    self.ffr_ramification_points.append(rp)

    def print_info(self):
        logger = logging.getLogger(self.logger_name)
        logger.info(
            "---------------------------------------------------------\n"
            "Branch Point at z = {}\n"
            "---------------------------------------------------------"
            .format(self.z)
        )
        for key, value in vars(self).iteritems():
            logger.info("{}:".format(key))
            pprint(value)


class IrregularSingularity:
    """
    The IrregularSingularity class.
    Just a container of information.
    Strictly related to the first fundamental representation cover.
    """
    def __init__(self, z=None, label=None, json_data=None):
        if json_data is None:
            self.z = z
            self.label = label
            self.monodromy = None
        else:
            self.z = r2toc(json_data['z'])
            self.label = json_data['label']
            self.monodromy = numpy.array(json_data['monodromy'])

        self.data_attributes = ['z', 'label', 'monodromy']

    def set_z_rotation(self, z_rotation):
        if self.z != oo:
            self.z /= complex(z_rotation)

    def get_json_data(self):
        json_data = {
            'z': ctor2(self.z),
            'label': self.label,
            'monodromy': self.monodromy.tolist(),
        }
        return json_data

    def print_info(self):
        logger.info(
            "---------------------------------------------------------\n"
            "Irregular singularity at z = {}\n"
            "---------------------------------------------------------"
            .format(self.z)
        )
        for key, value in vars(self).iteritems():
            logger.info("{}:".format(key))
            pprint(value)


# TODO: Use g_data.weights at the base point as labels of sheets,
# instead of integer indicies. Just a conceptual issue, because
# the integer indicies are labeling g_data.weights in the same order.
# PL: not sure if we want to do that: for a human interface,
# labeling by integers is much more readable.
class SWDataWithTrivialization(SWDataBase):
    """
    All branch cuts are assumed to run vertically, emanating
    upwards from branch points and irregular singularities.

    Arguments
    ---------

    sw_data :
        an object of the type SWData, whose attribute 'curve'
        should correspond to the curve in the FIRST fundamental
        representation of the Lie algebra

    ffr_ramification_points :
        a list of objects of the type RamificationPoint, corresponding
        to the given Seiberg-Witten curve in the first fundamental rep.

    Attributes & Methods
    --------------------

    base_point :
        the base point of the trivialization

    reference_ffr_xs :
        a list of x's
            [x_0, x_1, ..., x_i, ...]
        where 'i' is an integer label for the sheet,
        and 'x' is its position in the fiber of T^*C 
        over the basepoint. This is aligned with 
        g_data.ffr_weights.

    get_sheets_at_z(z) :
        this method returns the set of sheets and their integer label 
        identifier at any point 'z' on the C-plane.
        The labels are consistent with those at the basepoint.
        To get the corresponding weights, of the first fundamental 
        representation, use g_data.weights[i].
        The output looks like this
        {0 : x_0, ... , i : x_i, ...}
    """
    # NOTE: I am assuming that branch points NOR irregular singularities 
    # overlap vertically.
    # This should be guaranteed by the automatic rotation of 
    # the z-plane which is performed before calling this class.
    # TODO: irregular singularities ARE branch points.
    # PL: in what sense? 
    def __init__(self, config, logger_name='loom', json_data=None,):
        super(SWDataWithTrivialization, self).__init__(
            config, logger_name=logger_name, json_data=json_data,
        )

        self.branch_points = []
        self.irregular_singularities = []

        self.min_distance = None 
        self.min_horizontal_distance = None 
        self.base_point = None 
        self.reference_ffr_xs = None
        self.reference_xs = None
        # self.branch_cut_angle = None

        if json_data is None:
            self.set_trivialization()
        else:
            self.set_trivialization_from_json_data(json_data)

        self.data_attributes += [
            'branch_points', 'irregular_singularities', 'min_distance',
            'min_horizontal_distance', 'base_point', 'reference_ffr_xs',
            'reference_xs'
        ]

    def set_z_rotation(self, z_rotation):
        super(SWDataWithTrivialization, self).set_z_rotation(z_rotation)

        for p in (self.branch_points + self.irregular_singularities):
            p.set_z_rotation(z_rotation)

        self.base_point /= complex(z_rotation)

    def get_json_data(self):
        json_data = super(SWDataWithTrivialization, self).get_json_data()
        json_data['branch_points'] = [
            bp.get_json_data()
            for bp in self.branch_points
        ]
        json_data['irregular_singularities'] = [
            irs.get_json_data()
            for irs in self.irregular_singularities
        ]
        json_data['min_distance'] = self.min_distance
        json_data['min_horizontal_distance'] = self.min_horizontal_distance
        json_data['base_point'] = ctor2(self.base_point)
        json_data['reference_ffr_xs'] = [
            ctor2(x) for x in self.reference_ffr_xs
        ]
        json_data['reference_xs'] = [
            ctor2(x) for x in self.reference_xs
        ]
        json_data['farthest_branching_locus'] = ctor2(
            self.farthest_branching_locus
        )

        return json_data

    def set_trivialization_from_json_data(self, json_data,):
        self.branch_points = [
            BranchPoint(json_data=data, 
                        ffr_ramification_points=self.ffr_ramification_points,)
            for data in json_data['branch_points']
        ]
        self.irregular_singularities = [
            IrregularSingularity(json_data=data)
            for data in json_data['irregular_singularities']
        ]
        self.min_distance = json_data['min_distance'] 
        self.min_horizontal_distance = json_data['min_horizontal_distance']
        self.base_point = r2toc(json_data['base_point']) 
        self.reference_ffr_xs = [
            r2toc(x) for x in json_data['reference_ffr_xs']
        ]
        self.reference_xs = [
            r2toc(x) for x in json_data['reference_xs']
        ]
        try:
            self.farthest_branching_locus = r2toc(
                json_data['farthest_branching_locus']
            )
        except KeyError:
            logger.warning('No JSON data for farthest_branching_locus, '
                           'setting to None.')
            self.farthest_branching_locus = None

    def set_trivialization(self):
        logger = logging.getLogger(self.logger_name)

        # z-coords of branch points.
        bpzs = n_remove_duplicate(
            [r.z for r in self.ffr_ramification_points if r.z != oo],
            self.accuracy,
        )

        logger.info(
            'Ramification points arrange into {} branch points '
            'at positions {}'.format(
                len(bpzs), bpzs
            )
        )

        # z-coords of irregular punctures.
        ipzs = n_remove_duplicate(
            [p.z for p in self.irregular_punctures if p.z != oo],
            self.accuracy,
        )
        
        # Automatically choose a basepoint, based on the positions of
        # both branch points and irregular singularities
        all_points_z = bpzs + ipzs
        n_critical_loci = len(all_points_z)
        
        if n_critical_loci > 1:
            all_distances = [
                abs(x - y) for x in all_points_z for y in all_points_z
            ]
            max_distance = max(all_distances)
            # Minimun mutual distance among all the
            # branch points/punctures.
            non_zero_distances = [
                x for x in all_distances if abs(x) > self.accuracy
            ]
            self.min_distance = min(non_zero_distances)
            horizontal_distances = [
                abs(x.real - y.real) for x in all_points_z 
                for y in all_points_z
            ]
            self.min_horizontal_distance = min(
                [x for x in horizontal_distances if abs(x) > self.accuracy]
            )
            self.farthest_branching_locus = max(
                [abs(x) for x in all_points_z]
            )
            
        elif n_critical_loci == 1:
            # If there is only one branching locus, we still
            # need to set distance scales, as these will be used to 
            # circumvent the branch locus when constructing paths 
            # to trivializae the cover, as well as to fix a basepoint
            # for the trivialization
            max_distance = 3.0
            self.min_distance = max_distance
            self.min_horizontal_distance = max_distance
            self.farthest_branching_locus = max_distance

        elif n_critical_loci == 0:
            raise Exception('Must have at least one critical locus on C.')

        center = sum([z_pt for z_pt in all_points_z]) / n_critical_loci
        self.base_point = center - 1j * max_distance

        # Fix reference x's at the basepoint.
        # These sheets are aligned in the order of
        # sw.g_data.weights, i.e. reference_sheets[i]
        # is the value of x corresponding to 
        # sw.g_data.weights[i].
        logger.info(
            "Getting aligned x's at the base point z = {}."
            .format(self.base_point)
        )
        self.reference_ffr_xs, self.reference_xs = self.get_aligned_xs(
            self.base_point,
        )

        # Construct the list of branch points. 
        # Also assign them the respective ramification points 
        # and analyze them
        for i, z_bp in enumerate(bpzs):
            bp = BranchPoint(z=z_bp)
            bp.label = 'Branch point #{}'.format(i)

            self.analyze_branch_point(bp)
            if bp.order > 1:
                # only add if there are any positive roots associated
                # otherwise may be an accidental BP
                # FIXME: Must handle also accidental BP; for example 
                # a point like F~ z + x^2(1+x^2) can happen in D-type
                # and will have no obvious monodromy. Need to deal with it. 
                self.branch_points.append(bp)

        # Construct the list of irregular singularities
        for j, z_irr_sing in enumerate(ipzs):
            irr_sing = IrregularSingularity(
                z=z_irr_sing, label='Irr.Sing. #{}'.format(j)
            )
            self.analyze_irregular_singularity(irr_sing)
            self.irregular_singularities.append(irr_sing)

        # perform extensive global checks on monodromies
        logger.info('Doing global checks on monodromies...')
        self.global_monodromy_checks()

    def global_monodromy_checks(self):
        """
        A set of global checks on monodromies.
        """
        logger = logging.getLogger(self.logger_name)
        rep_d = len(self.g_data.weights)
        ram_loci = self.branch_points + self.irregular_singularities
        ordered_monodromies = [rl.monodromy for rl in sorted(
            ram_loci, key=lambda rl: rl.z.real
        )]

        # If there is no irregular singularity at infinity, 
        # we make a check that the numerically computed monodromy is trivial
        m_oo = None
        for p in self.irregular_punctures: 
            if p.z == oo:
                m_oo = True
                logger.debug(
                    "There is branching at infinity, "
                    "due to an irregular singularity."
                )

        path_around_oo = get_path_around(oo, self.base_point, self)
        if m_oo is None:
            numerical_m_oo = self.get_sheet_monodromy(path_around_oo)
        elif m_oo is True:
            numerical_m_oo = self.get_sheet_monodromy(
                path_around_oo, is_irr_sing=True
            )

        if m_oo is None:
            if numpy.array_equal(
                numerical_m_oo, numpy.identity(rep_d)
            ) is True:
                m_oo = numpy.identity(rep_d)
                logger.info(
                    "Direct computation confirms trivial monodromy at oo"
                )
            else:
                raise Exception(
                    "Monodromy at infinity is nontrivial, but no branching "
                    "has been declared at infinity."
                )
        elif m_oo is True:
            m_oo = numerical_m_oo

        # Check that product of all monodromies equals 
        # the monodromy at infinity    
        m_tot = numpy.identity(rep_d)
        for m_i in ordered_monodromies:
            # NOTE: since monodromies are counter-clockwise and 
            # ramification loci are ordered according to ascending 
            # real part of their coordinate, multiplication must 
            # be to the right.
            m_tot = numpy.dot(m_tot, m_i)

        if numpy.array_equal(m_tot, m_oo) is True:
            logger.debug(
                "Monodromy at infinity agrees with the product "
                "of all monodromies, as expected."
            )
        else:
            logger.info(
                '\n\nWARNING!\n'
                'Monodromy at infinity does not agree with the product '
                'of all monodromies. Probable error in trivialization.'
            )
            # raise Exception(
            #     'Monodromy at infinity does not agree with the product '
            #     'of all monodromies. Probable error in trivialization.'
            # )

        # Now we check consistency for various combinations of monodromies
        # e.g. M_1 M_2 ... M_k = M_{k+1} ... M_N (M_{oo})^{-1}
        # for all possible sequential pairwise splittings.
        for i in range(len(ordered_monodromies)):
            m_1 = numpy.identity(rep_d)
            m_2 = numpy.identity(rep_d)
            for m_i in ordered_monodromies[:(i + 1)]:
                # NOTE: since monodromies are counter-clockwise and 
                # ramification loci are ordered according to ascending 
                # real part of their coordinate, multiplication must 
                # be to the right.
                m_1 = numpy.dot(m_1, m_i)
            for m_i in ordered_monodromies[(i + 1):]:
                # NOTE: since monodromies are counter-clockwise and 
                # ramification loci are ordered according to ascending 
                # real part of their coordinate, multiplication must 
                # be to the right.
                m_2 = numpy.dot(m_2, m_i)
            if numpy.array_equal(numpy.dot(m_1, m_2), m_oo):
                continue
            else:
                logger.info(
                    '\n\nWARNING!\n'
                    'The product of the first {} monodromies '
                    'does not match the (inverse of) the remaining ones'
                    .format(i + 1)
                )
                # raise Exception(
                #     'The product of the first {} monodromies '
                #     'does not match the (inverse of) the remaining ones'
                #     .format(i + 1)
                # )


    def get_sheets_along_path(
        self, z_path, is_path_to_bp=False, ffr=False,
        ffr_xs_0=None, zoom_level=MAX_ZOOM_LEVEL,
        accuracy=None, ffr_sheets_along_path=None,
    ):
        """
        Tracks the sheets along a path.
        It checks at each step that tracking is successful,
        meaning that all sheets can be distinguished correctly.
        This would fail if we choose a path ending on a branch-point.
        For tracking roots as we run into a branch point, one should
        set the variable 'is_path_to_bp=True', and the check for 
        sheets becoming too similar will be ignored altogether.
        If tracking fails, an attempt will be made to 'zoom in',
        up to a certain number of times.
        If zooming also fails, the tracking will take into account 
        the first derivative of sheets and match them according to it.
        """       
        logger = logging.getLogger(self.logger_name)
        g_data = self.g_data
        if accuracy is None:
            accuracy = self.accuracy
        # If the initial sheets are unspecified, 
        # the initial point should be the basepoint of the trivialization 
        if ffr_xs_0 is None:
            if abs(z_path[0] - self.base_point) < self.accuracy:
                ffr_xs_0 = self.reference_ffr_xs
                # xs_0 = self.reference_xs
            else:
                raise Exception('Must specify initial sheets for tracking.')
        logger.debug(
            'Zooming level: {}/{}'.format(zoom_level, MAX_ZOOM_LEVEL)
        )
        # Each element is a sheet, which is a list of x's along the path.
        # Initialized with reference_xs.
        # TODO: set each element to an integer rather than a float.
        
        if ffr_sheets_along_path is None:
            ffr_sheets_along_path = [[x] for x in ffr_xs_0]
        
        for i, z in enumerate(z_path):
            near_degenerate_branch_locus = False
            if is_path_to_bp is True and abs(z - z_path[-1]) < self.accuracy:
                near_degenerate_branch_locus = True
            # NOTE: Commenting this, since we don't need to bother
            # with alignment for now. This saves a lot of problems 
            # with numerics and make the code much faster.
            # Delete the following commend after EXTENSIVE testing.
            # ffr_xs_1, xs_1 = self.get_aligned_xs(
            #     z, 
            #     near_degenerate_branch_locus=near_degenerate_branch_locus
            # )
            ffr_xs_1 = self.ffr_curve.get_xs(z) 

            # if it's not a path to branch point, check tracking
            if is_path_to_bp is False:
                sorted_ffr_xs = get_sorted_xs(
                    ffr_xs_0, ffr_xs_1, 
                    accuracy=accuracy,
                    check_tracking=True, index=i,
                    z_0=z_path[i - 1], z_1=z_path[i],
                    g_data=g_data,
                    logger_name=self.logger_name,
                    sw_curve=self.curve,
                )
            # if it's a path to branch point, but we are far from it,
            # still check tracking
            elif near_degenerate_branch_locus is False:
                sorted_ffr_xs = get_sorted_xs(
                    ffr_xs_0, ffr_xs_1,
                    accuracy=accuracy,
                    check_tracking=True,
                    z_0=z_path[i - 1], z_1=z_path[i],
                    g_data=g_data,
                    logger_name=self.logger_name,
                    sw_curve=self.curve,
                )
            # if it's a path to a branch point and we are getting close to it, 
            # don't check tracking anymore
            else:
                sorted_ffr_xs = get_sorted_xs(
                    ffr_xs_0, ffr_xs_1,
                    accuracy=accuracy,
                    check_tracking=False,
                    g_data=g_data,
                    logger_name=self.logger_name,
                    sw_curve=self.curve,
                )
            if sorted_ffr_xs == 'sorting failed':
                logger.warning(
                    'Encountered a problem with sheet tracking.'
                )
                if zoom_level > 0:
                    if i > 0:
                        delta_z = (z_path[i] - z_path[i - 1]) / ZOOM_FACTOR
                        zoomed_path = [
                            z_path[i - 1] + j * delta_z  
                            for j in range(ZOOM_FACTOR + 1)
                        ]
                    else:
                        delta_z = (z_path[1] - z_path[0]) / ZOOM_FACTOR
                        zoomed_path = [
                            z_path[0] + j * delta_z  
                            for j in range(ZOOM_FACTOR + 1)
                        ]
                    sheets_along_zoomed_path = self.get_sheets_along_path(
                        zoomed_path, 
                        is_path_to_bp=near_degenerate_branch_locus,
                        ffr=True,
                        ffr_xs_0=ffr_xs_0,
                        zoom_level=(zoom_level - 1),
                        accuracy=(accuracy / ZOOM_FACTOR),
                        ffr_sheets_along_path=ffr_sheets_along_path,
                    )
                    # Just keep the last tracked value for each sheet 
                    sorted_ffr_xs = [
                        zoom_s[-1] for zoom_s in sheets_along_zoomed_path
                    ]
                    # for j, s_j in enumerate(ffr_sheets_along_path):
                    #     s_j += sorted_ffr_xs[j]

                else:
                    old_ffr_xs = [s[-2] for s in ffr_sheets_along_path]
                    delta_xs = [
                        ffr_xs_0[j] - old_ffr_xs[j] 
                        for j in range(len(ffr_xs_0))
                    ]
                    sorted_ffr_xs = sort_xs_by_derivative(
                        ffr_xs_0, ffr_xs_1, delta_xs, self.accuracy,
                        logger_name=self.logger_name,
                    )
                    if sorted_ffr_xs == 'sorting failed':
                        ffr_xs_1_s = self.ffr_curve.get_xs(z, use_sage=True)
                        sorted_ffr_xs = sort_xs_by_derivative(
                            ffr_xs_0, ffr_xs_1_s, delta_xs, self.accuracy,
                            logger_name=self.logger_name,
                        )
                        if sorted_ffr_xs == 'sorting failed':
                            logger.info(
                                'Studying sheets near z = {} found sheets'
                                '\n ffr_xs_0 = {} \n ffr_xs_1 = {}'
                                .format(z, ffr_xs_0, ffr_xs_1_s)
                            )
                            raise Exception(
                                '\nCannot track the sheets!\n'
                                'Probably passing too close to a branch point'
                                ' or a puncture. Try increasing N_PATH_TO_PT '
                                'or N_PATH_AROUND_PT, or MAX_ZOOM_LEVEL.'
                            )

            # this is just the ordinary step, where we add the 
            # latest value of ordered sheets
            for j, s_j in enumerate(ffr_sheets_along_path):
                s_j.append(sorted_ffr_xs[j])

            ffr_xs_0 = sorted_ffr_xs
            
        # the result is of the form [sheet_path_1, sheet_path_2, ...]
        # where sheet_path_i = [x_0, x_1, ...] are the fiber coordinates
        # of the sheet along the path
        if ffr is True:
            return ffr_sheets_along_path
        elif ffr is False:
            rep_dim = len(g_data.weights)
            # Warning: the following is offset by 1 from len(z_path)!
            n_steps = len(ffr_sheets_along_path[0])
            sheets_along_path = [
                [0 for j in range(n_steps)] for i in range(rep_dim)
            ]

            for i in range(n_steps):
                # Note: here it is crucial that the sheets are 
                # ordered according to the weights of the rep
                # to quich they correspond
                ffr_xs_at_step_i = [s[i] for s in ffr_sheets_along_path]
                xs_at_step_i = self.get_xs_of_weights_from_ffr_xs(
                    ffr_xs_at_step_i
                )
                for j, sheet_track_j in enumerate(sheets_along_path):
                    sheet_track_j[i] = xs_at_step_i[j]

            return sheets_along_path

    def get_sheets_at_z(self, z_pt, g_data=None, ffr=False):
        """
        Returns a dict of (sheet_index, x) at a point ''z_pt'', 
        which cannot be a branch point or a singularity.
        """
        z_path = get_path_to(z_pt, self)
        sheets = self.get_sheets_along_path(z_path, ffr=ffr)
        final_xs = [s_i[-1] for s_i in sheets]
        final_sheets = {i: x for i, x in enumerate(final_xs)}
        return final_sheets

    # TODO: Review this method.
    def get_sheet_monodromy(
        self, z_path, is_higher_bp=False, higher_bp_type=None,
        is_irr_sing=False
    ):
        """
        Compares the x-coordinates of sheets at the 
        beginning and at the end of a CLOSED path.
        Returns a permutation matrix, expressed in 
        the basis of reference sheets, such that
        new_sheets = M . old_sheets

        There are subtleties associarted with representations
        of D and E-type algebras, where two or more sheets
        can be identically zero. Then one needs to 
        figure out how to permute them, if at all.
        For this reason, it is important to specify whenever 
        we are computing the monodromy around an irregular singularity.
        """
        logger = logging.getLogger(self.logger_name)
        logger.debug(
            "Analyzing the monodromy around a closed path "
            "of length {}.".format(len(z_path))
        )
        initial_xs = self.reference_xs
        initial_sheets = [[i, x] for i, x in enumerate(initial_xs)]
        sheet_tracks = self.get_sheets_along_path(z_path)
        final_xs = [s_t[-1] for s_t in sheet_tracks]
        final_sheets = [[i, x] for i, x in enumerate(final_xs)]

        # Now we compare the initial and final sheets 
        # to extract the monodromy permutation
        # recall that each entry of initial_sheets and final_sheets
        # is of the form [i, x] with i the integer label
        # and x the actual position of the sheet in the fiber 
        # above the basepoint.
        sorted_sheets = []
        for s_1 in initial_sheets:
            closest_candidate = final_sheets[0]
            min_d = abs(s_1[1] - closest_candidate[1])
            for s_2 in final_sheets:
                if abs(s_2[1] - s_1[1]) < min_d:
                    min_d = abs(s_2[1] - s_1[1])
                    closest_candidate = s_2
            sorted_sheets.append(closest_candidate)
        
        # Now we check that sheet tracking is not making a mistake.
        # NOTE: cannot use the function 'delete_duplicates' with this 
        # data structure.
        seen = set()
        uniq = []
        for s in sorted_sheets:
            if s[1] not in seen:
                uniq.append(s[1])
                seen.add(s[1])

        if len(uniq) < len(sorted_sheets):
            # When studying D-type covers there may be situations
            # where two sheets collide at x=0 everywhere
            # Do not raise an error in this case.
            # Likewise for E-type
            if (
                min(map(abs, [s[1] for s in sorted_sheets])) < self.accuracy 
            ) and (
                self.g_data.type == 'D' and len(sorted_sheets) - len(uniq) == 1
            ) or (
                self.g_data.type == 'E' and len(sorted_sheets) - len(uniq) == 2
            ):
                # If two or more sheets are equal (and both zero) then the 
                # (integer) labels they got assigned in sorting above may 
                # be the same, this would lead to a singular permutation matrix
                # and must be corrected, as follows.
                int_labels = [s[0] for s in sorted_sheets]
                uniq_labels = list(set(int_labels))
                labels_multiplicities = [
                    len([i for i, x in enumerate(int_labels) if x == u]) 
                    for u in uniq_labels
                ]
                multiple_labels = []
                for i, u in enumerate(uniq_labels):
                    if labels_multiplicities[i] > 1:
                        if (
                            labels_multiplicities[i] == 2 
                            and self.g_data.type == 'D'
                        ) or (
                            labels_multiplicities[i] == 3 
                            and self.g_data.type == 'E'
                        ):
                            multiple_labels.append(u)
                        else:
                            logger.debug(
                                'int labels = {}'.format(int_labels)
                            )
                            logger.debug(
                                'multiple labels = {}'
                                .format(multiple_labels)
                            )
                            raise Exception('Too many degenerate sheets')
                if len(multiple_labels) != 1:
                    raise Exception(
                        'Cannot determine which sheets are' +
                        'degenerate, tracking will fail.'
                    )

                if self.g_data.type == 'D':
                    double_sheets = [
                        i for i, s in enumerate(sorted_sheets) 
                        if s[0] == multiple_labels[0]
                    ]

                    corrected_sheets = sorted_sheets 
                    if is_higher_bp is False and is_irr_sing is False:
                        corrected_sheets[double_sheets[0]] = (
                            initial_sheets[double_sheets[0]]
                        )
                        corrected_sheets[double_sheets[1]] = (
                            initial_sheets[double_sheets[1]]
                        )
                    elif (
                        is_higher_bp is True and (
                            higher_bp_type == 'type_II' 
                            or higher_bp_type == 'type_III'
                            # FIXME: remove this after removing 'type_AD' 
                            # placeholder
                            or higher_bp_type == 'type_AD'

                        ) or is_irr_sing is True
                    ):
                        # Should decide case-by-case whether to employ 
                        # (0,1) -> (0,1) or (0,1) -> (1,0)
                        # One way to do so would be to pick each of them, and 
                        # construct the whole monodromy, then see if applying 
                        # the monodromy to every root gives back a root
                        # By direct checks, only one of the two options works
                        # so there should be no ambiguity left.

                        # TODO: print a warning if BOTH options give a Weyl 
                        # sheet matrix, because in that case there may be 
                        # ambiguity

                        # UPDATE: Disabling the option 0. Because in these 
                        # types  of ramification points there should be no 
                        # sheet that remains fixed by the permutation.
                        # Keep it in comment in case we encounter trouble

                        # corrected_sheets_0 = [s for s in corrected_sheets]
                        corrected_sheets_1 = [s for s in corrected_sheets]

                        # corrected_sheets_0[double_sheets[0]] = (
                        #     initial_sheets[double_sheets[0]]
                        # )
                        # corrected_sheets_0[double_sheets[1]] = (
                        #     initial_sheets[double_sheets[1]]
                        # )

                        corrected_sheets_1[double_sheets[0]] = (
                            initial_sheets[double_sheets[1]]
                        )
                        corrected_sheets_1[double_sheets[1]] = (
                            initial_sheets[double_sheets[0]]
                        )

                        # m_0 = build_monodromy_matrix(
                        #     initial_sheets, corrected_sheets_0
                        # )
                        m_1 = build_monodromy_matrix(
                            initial_sheets, corrected_sheets_1
                        )

                        # if is_weyl_monodromy(m_0, self.g_data):
                        #     corrected_sheets = corrected_sheets_0
                        if is_weyl_monodromy(m_1, self.g_data):
                            corrected_sheets = corrected_sheets_1
                        else:
                            raise Exception(
                                'Failed to assign a Weyl-type monodromy.'
                            )
                        
                    else:
                        raise Exception(
                            'higher-type ramification points for D-type '
                            'theories can only be of types II or III. '
                            'Found {} instead'.format(higher_bp_type)
                        )
                    sorted_sheets = corrected_sheets
                    pass

                elif self.g_data.type == 'E':
                    triple_sheets = [
                        i for i, s in enumerate(sorted_sheets) 
                        if s[0] == multiple_labels[0]
                    ]

                    corrected_sheets = sorted_sheets 
                    if is_higher_bp is False and is_irr_sing is False:
                        corrected_sheets[triple_sheets[0]] = (
                            initial_sheets[triple_sheets[0]]
                        )
                        corrected_sheets[triple_sheets[1]] = (
                            initial_sheets[triple_sheets[1]]
                        )
                        corrected_sheets[triple_sheets[2]] = (
                            initial_sheets[triple_sheets[2]]
                        )

                    elif (
                        (is_higher_bp is True and higher_bp_type == 'type_IV')
                        or is_irr_sing is True
                    ):
                        # Should decide case-by-case whether to employ 
                        # (0,1,2) -> (1,2,0) or (0,1,2) -> (2,0,1)
                        # One way to do so would be to pick each of them, and 
                        # construct the whole monodromy, then see if applying 
                        # the monodromy to every root gives back a root
                        # By direct checks, only one of the two options works
                        # so there should be no ambiguity left.
                        
                        # TODO: print a warning if BOTH options give a Weyl 
                        # sheet matrix, because in that case there may be 
                        # ambiguity

                        # UPDATE: Disabling the option 0. 
                        # Because in these types of ramification points
                        # there should be no sheet that remains fixed by 
                        # the permutation.
                        # Keep this in comment in case we encounter trouble

                        # corrected_sheets_0 = [s for s in corrected_sheets]
                        corrected_sheets_1 = [s for s in corrected_sheets]
                        corrected_sheets_2 = [s for s in corrected_sheets]

                        # corrected_sheets_0[triple_sheets[0]] = (
                        #     initial_sheets[triple_sheets[0]]
                        # )
                        # corrected_sheets_0[triple_sheets[1]] = (
                        #     initial_sheets[triple_sheets[1]]
                        # )
                        # corrected_sheets_0[triple_sheets[2]] = (
                        #     initial_sheets[triple_sheets[2]]
                        # )

                        corrected_sheets_1[triple_sheets[0]] = (
                            initial_sheets[triple_sheets[1]]
                        )
                        corrected_sheets_1[triple_sheets[1]] = (
                            initial_sheets[triple_sheets[2]]
                        )
                        corrected_sheets_1[triple_sheets[2]] = (
                            initial_sheets[triple_sheets[0]]
                        )

                        corrected_sheets_2[triple_sheets[0]] = (
                            initial_sheets[triple_sheets[2]]
                        )
                        corrected_sheets_2[triple_sheets[1]] = (
                            initial_sheets[triple_sheets[0]]
                        )
                        corrected_sheets_2[triple_sheets[2]] = (
                            initial_sheets[triple_sheets[1]]
                        )

                        # m_0 = build_monodromy_matrix(
                        #     initial_sheets, corrected_sheets_0
                        # )
                        m_1 = build_monodromy_matrix(
                            initial_sheets, corrected_sheets_1
                        )
                        m_2 = build_monodromy_matrix(
                            initial_sheets, corrected_sheets_2
                        )

                        # if is_weyl_monodromy(m_0, self.g_data):
                        #     corrected_sheets = corrected_sheets_0
                        if is_weyl_monodromy(m_1, self.g_data):
                            corrected_sheets = corrected_sheets_1
                        elif is_weyl_monodromy(m_2, self.g_data):
                            corrected_sheets = corrected_sheets_2
                        else:
                            raise Exception(
                                'Failed to assign a Weyl-type monodromy.'
                            )
                            
                    else:
                        raise Exception(
                            'higher-type ramification points for E-type '
                            'theories can only be of type IV. Found {}'
                            .format(higher_bp_type)
                        )
                    sorted_sheets = corrected_sheets
                    pass

            else:
                raise ValueError(
                    '\nError in determination of monodromy!\n' +
                    'Cannot match uniquely the initial sheets' + 
                    ' to the final ones.'
                )
        else:
            pass

        logger.debug(
            'Sorted sheets around locus {}'.format(sorted_sheets)
        )

        perm_matrix = build_monodromy_matrix(initial_sheets, sorted_sheets)
        if is_weyl_monodromy(perm_matrix, self.g_data) is False:
            raise Exception('Failed to assign a Weyl-type monodromy.')
        else:
            logger.debug('Sheet monodromy is of Weyl type')
        logger.debug(
            'Permutation matrix {}'.format(perm_matrix)
        )
        return perm_matrix

    def analyze_branch_point(self, bp):
        logger = logging.getLogger(self.logger_name)
        logger.info(
            "Analyzing a branch point at z = {}."
            .format(bp.z)
        )
        logger.debug(
            "Studying groups of colliding sheets to determine root(s) type."
        )
        path_to_bp = get_path_to(bp.z, self)
        # it's important that we set ffr=False, as the sheets here
        # will be used along with the algorithm for matching them to
        # roots, and there the sheets will be assumed to be in 1-1
        # correspondence (and ordered) with the weights of the representation
        sheets_along_path = self.get_sheets_along_path(
            path_to_bp, is_path_to_bp=True, ffr=False,
        )
        xs_at_bp = [s_i[-1] for s_i in sheets_along_path]
        enum_sh = [[i, x_i] for i, x_i in enumerate(xs_at_bp)]
        
        clusters = []
        for i, x in enum_sh:
            is_single = True
            for c_index, c in enumerate(clusters):
                x_belongs_to_c = belongs_to_cluster(x, c, enum_sh)
                if x_belongs_to_c is True:
                    clusters[c_index].append(i)
                    is_single = False
                    break
            if is_single is True:
                clusters.append([i])

        # bp.enum_sh = enum_sh
        bp.groups = [c for c in clusters if len(c) > 1]
        # bp.singles = [c[0] for c in clusters if len(c) == 1]

        bp.positive_roots = get_positive_roots_of_branch_point(
            bp, self.g_data, self.logger_name,
        )
        bp.order = len(bp.positive_roots) + 1

        bp.ffr_ramification_points = [
            rp for rp in self.ffr_ramification_points
            if abs(rp.z - bp.z) < self.accuracy
        ]

        ramification_types = ([
            rp.ramification_type for rp in bp.ffr_ramification_points 
            if rp.ramification_type != 'type_I'
        ])

        logger.debug("Computing the monodromy.")
        path_around_bp = get_path_around(bp.z, self.base_point, self,)

        if len(ramification_types) == 0:
            bp.monodromy = self.get_sheet_monodromy(path_around_bp)
        else:
            if len(delete_duplicates(ramification_types)) == 1:
                bp.monodromy = self.get_sheet_monodromy(
                    path_around_bp, 
                    is_higher_bp=True, higher_bp_type=ramification_types[0]
                )
            else:
                raise Exception(
                    'Multiple ramification types for BP at z = {}'.format(bp.z)
                )

        # If the order of the branch point is 1, it means that it does 
        # not correspond to a root, and will be dropped later on.
        # Since all branch points should be of weyl-type, there
        # should be no monodromy around this branch point.
        # Here we check there is no monodromy at all around this branch point.
        if bp.order == 1:
            rep_d = len(self.g_data.weights)
            if numpy.array_equal(bp.monodromy, numpy.identity(rep_d)) is True:
                logger.debug(
                    "This branch point has trivial monodromy, in agreement"
                    "with the fact that it does not correspond to any root."
                )
            else:
                logger.info(
                    "Warning! While the branch point does not correspond to "
                    "any root, its monodromy is nontrivial. "
                    "Probable error in trivialization of the covering."
                )

        # TODO: it would be good to make a check here, e.g. on the 
        # relation between vanishing roots and the monodromy.

        # XXX: Temporary analysis for D-type AD theories.
        # TODO: Delete this, 'for' cycle, and implement this fix 
        # inside the routine get_sheet_monodromy as for all the other types.
        for rp in bp.ffr_ramification_points:
            if (
                rp.ramification_type == 'type_AD' 
                # or rp.ramification_type == 'type_III'
                # NOTE: no need to do this for type_III, already handled
                # in the routine get_sheet_monodromy
                # But, double-check me with explicit examples.
            ):
                logger.warning('Temporary fix-up of the monodromy '
                               'for a D-type AD theory.')
                # Find two diagonal elements and permute them.
                indices = []
                for i, row_i in enumerate(bp.monodromy):
                    if row_i[i] == 1:
                        indices.append(i)
                if len(indices) == 2:
                    i, j = indices
                    bp.monodromy[i][i] = 0
                    bp.monodromy[j][j] = 0
                    bp.monodromy[i][j] = 1
                    bp.monodromy[j][i] = 1
                else:
                    RuntimeError('Unknown form of the monodromy at {}.'
                                 .format(bp.label))

    def analyze_irregular_singularity(self, irr_sing):
        logger = logging.getLogger(self.logger_name)
        logger.info(
            "Analyzing an irregular singularity at z = {}."
            .format(irr_sing.z)
        )
        path_around_z = get_path_around(irr_sing.z, self.base_point, self)
        irr_sing.monodromy = self.get_sheet_monodromy(
            path_around_z, is_irr_sing=True
        )
    

def get_path_to(z_pt, sw_data, logger_name='loom'):
    """
    Return a rectangular path from the base point to z_pt.
    If the path has to pass too close to a branch point, 
    we avoid the latter by drawing an arc around it.
    """
    logger = logging.getLogger(logger_name)
    base_pt = sw_data.base_point
    closest_bp = None
    # if n_loci==None:
    #     n_loci = len(sw_data.branch_points + sw_data.irregular_singularities)
    # radius = sw_data.min_distance / n_loci
    radius = sw_data.min_horizontal_distance / 2.0

    logger.debug("Constructing a path [{}, {}]".format(base_pt, z_pt))

    # Determine if the path will need to pass 
    # close to a branch point.
    for bp in sw_data.branch_points:
        delta_z = z_pt - bp.z
        # NOTE we only check for one possible nearby point
        # based on the fact that the radius is always less
        # than the minimal horizontal separation of them
        if abs(delta_z.real) < radius and delta_z.imag > 0:
            closest_bp = bp
            break

    # If there the path does not pass near a branch point:
    if closest_bp is None:
        z_0 = base_pt
        z_1 = 1j * base_pt.imag + z_pt.real
        z_2 = z_pt
        half_steps = int(N_PATH_TO_PT / 2)
        return (
            [z_0 + ((z_1 - z_0) / half_steps) * i 
             for i in range(half_steps + 1)] + 
            [z_1 + ((z_2 - z_1) / half_steps) * i 
             for i in range(half_steps + 1)]                
        )

    # If there the path needs to pass near a branch point:
    else:
        z_0 = base_pt
        z_1 = 1j * base_pt.imag + closest_bp.z.real
        z_2 = 1j * (closest_bp.z.imag - radius) + closest_bp.z.real
        z_3 = closest_bp.z + radius * exp(1j * phase(z_pt - closest_bp.z))
        z_4 = z_pt
        
        if (z_pt - closest_bp.z).real > 0:
            # way_around = 'ccw'
            sign = 1.0
            delta_theta = phase(z_pt - closest_bp.z) + pi / 2
        else:
            # way_around = 'cw'
            sign = -1.0
            delta_theta = 3 * pi / 2 - phase(z_pt - closest_bp.z) 

        steps = int(N_PATH_TO_PT / 5)

        path_segment_1 = [z_0 + ((z_1 - z_0) / steps) * i
                          for i in range(steps + 1)]
        path_segment_2 = [z_1 + ((z_2 - z_1) / steps) * i 
                          for i in range(steps + 1)]
        path_segment_3 = [
            closest_bp.z + radius * (-1j) * exp(
                sign * 1j * (delta_theta) * (float(i) / float(steps))
            ) for i in range(steps + 1)
        ]
        path_segment_4 = [
            z_3 + ((z_4 - z_3) / steps) * i for i in range(steps + 1)
        ]
        
        return (
            path_segment_1 + path_segment_2 + path_segment_3 + path_segment_4
        )
    

def get_path_around(z_pt, base_pt, sw):
    logger = logging.getLogger(sw.logger_name)
    logger.debug("Constructing a closed path around z = {}".format(z_pt))
    z_0 = base_pt
    
    if z_pt != oo:
        z_1 = 1j * base_pt.imag + z_pt.real
        # if n_loci==None:
        #     n_loci = len(sw.branch_points + sw.irregular_singularities)
        # radius = min_distance / n_loci
        radius = sw.min_horizontal_distance / 2.0
        z_2 = z_pt - 1j * radius

        steps = N_PATH_AROUND_PT
        path_segment_1 = [
            z_0 + ((z_1 - z_0) / steps) * i for i in range(steps + 1)
        ]
        path_segment_2 = [
            z_1 + ((z_2 - z_1) / steps) * i for i in range(steps + 1)
        ]
        path_segment_3 = [
            z_pt + radius * (-1j) * exp(i * 2.0 * pi * 1j / steps) 
            for i in range(steps + 1)
        ]
        path_segment_4 = path_segment_2[::-1]
        path_segment_5 = path_segment_1[::-1]
        return (
            path_segment_1 + path_segment_2 + path_segment_3 + 
            path_segment_4 + path_segment_5
        )
    else:
        radius = sw.farthest_branching_locus * 2.0
        steps = N_PATH_AROUND_PT
        # need more steps to get good tracking on a large circle
        large_steps = 10 * steps

        z_1 = 1j * base_pt.imag
        z_2 = -1j * radius
        path_segment_1 = [
            z_0 + ((z_1 - z_0) / steps) * i for i in range(steps + 1)
        ]
        path_segment_2 = [
            z_1 + ((z_2 - z_1) / steps) * i for i in range(steps + 1)
        ]
        path_segment_3 = [
            radius * (-1j) * exp(i * 2.0 * pi * 1j / large_steps) 
            for i in range(large_steps + 1)
        ]
        path_segment_4 = path_segment_2[::-1]
        path_segment_5 = path_segment_1[::-1]
        return (
            path_segment_1 + path_segment_2 + path_segment_3 + 
            path_segment_4 + path_segment_5
        )



# TODO: Try using numba.
# TODO: Make smarter checks based on the types
# of ramification points above the branch point.
def get_sorted_xs(ref_xs, new_xs, accuracy=None, check_tracking=True, 
                  index=None, z_0=None, z_1=None, g_data=None,
                  logger_name='loom', sw_curve=None):
    """
    Returns a sorted version of 'new_xs'
    based on matching the closest points with 
    'ref_xs'
    """
    logger = logging.getLogger(logger_name)
    sorted_xs = []

    # It may happen that numerical inaccuracy of numpy
    # leads to believe that sheet tracking is going wrong.
    # For example in E_6 we expect to have three sheets with x=0
    # but sometimes numpy returns as many as 14 sheets with x=0
    # then the above check will fail because the unique sheets are
    # far less than 27 - (3-1), to begin with.
    # In such case, one can use SAGE to get a better accuracy for 
    # the values of sheets x's. But sage is slow, so keep it for 
    # neceessary cases.
    # Maybe introduce a similar check for D-type, but it seems 
    # not necessary
    if g_data.type == 'E' and check_tracking is True:
        if (27 - len(n_remove_duplicate(new_xs, accuracy))) > 2:
            new_xs = sw_curve.get_xs(z_1, use_sage=True)
            if (27 - len(n_remove_duplicate(new_xs, accuracy))) > 2:
                logging.info(
                    'Warning! at z={} sheets seem too degenerate:\n{}'
                    .format(z_1, new_xs)
                )

    for s_1 in ref_xs:
        closest_candidate = nsmallest(1, new_xs, key=lambda x: abs(x - s_1))[0]
        sorted_xs.append(closest_candidate)
    
    if check_tracking is True:
        # Now we check that sheet tracking is not making a mistake.
        unique_sorted_xs = n_remove_duplicate(sorted_xs, accuracy)

        max_dx = max(
            [abs(x_i - x_j) for x_i in sorted_xs for x_j in sorted_xs]
        )
        max_dx_dt = max(
            [abs(sorted_xs[i] - ref_xs[i]) for i in range(len(ref_xs))]
        )

        if len(unique_sorted_xs) < len(sorted_xs):
            # When studying D-type covers there may be situations
            # where two sheets collide at x=0 everywhere
            # Do not raise an error in this case.
            # The same is true for E-type covers at the origin of the
            # Coulomb branch
            if (
                g_data.type == 'D' and min(map(abs, sorted_xs)) < accuracy
                and len(sorted_xs) - len(unique_sorted_xs) == 1
            ) or (
                g_data.type == 'E' and min(map(abs, sorted_xs)) < accuracy
                and len(sorted_xs) - len(unique_sorted_xs) == 2
            ):
                return sorted_xs
            else:
                logger.debug(
                    "At step %s, between %s and %s " % (index, z_0, z_1)
                )
                logger.debug("ref_xs:\n{}".format(ref_xs)) 
                logger.debug("new_xs:\n{}".format(new_xs)) 
                logger.debug("sorted_xs:\n{}".format(sorted_xs)) 
                logger.debug("unique_sorted_xs:\n{}".format(unique_sorted_xs)) 
                logger.debug('Having trouble tracking sheets, will zoom in.')
                return 'sorting failed'
        # don't trust tracking if dx/dt for a single step is larger than
        # the maximum difference between any two xs
        elif max_dx_dt > max_dx:
            return 'sorting failed'
        else:
            return sorted_xs
    else:
        # If the path is one ending on a branch-point, 
        # the check that tracking is correct is disabled
        # because it would produce an error, since by definition
        # sheets will be indistinguishable at the very end.
        return sorted_xs


def sort_xs_by_derivative(ref_xs, new_xs, delta_xs, accuracy, 
                          logger_name='loom'):
    # will only work if there are at most two sheets being 
    # too close two each other, not three or more.
    # Unless there are three or more sheets all equal to zero
    # In this case we assume it's a degenerate curve and we 
    # sort those sheets accordingly.
    # TODO: generalize to handle more general cases (if we need it at all)
    logger = logging.getLogger(logger_name)
    logger.debug('Resorting to tracking sheets by their derivatives')

    # first, identify the problematic sheets
    ys = []
    for s_1 in ref_xs:
        closest_candidate = nsmallest(1, new_xs, key=lambda x: abs(x - s_1))[0]
        ys.append(closest_candidate)
        
    # the list of ys corresponds to the ref_xs as
    # ref_xs = [x_1, x_2, x_3, ...]
    #     ys = [y_1, y_2, y_1, ...]
    # and will contain doubles. 
    # We use them to identify the pairs that give trouble
    correct_xy_pairs = {}   # (a dictionary)
    trouble_xs = []
    trouble_ys = []
    for i in range(len(ys)):
        if ys.count(ys[i]) == 1:
            # add this key and value to the dictioanry
            correct_xy_pairs.update({ref_xs[i]: ys[i]})
        else:
            trouble_ys.append(ys[i])

    trouble_ys = n_remove_duplicate(trouble_ys, 0.0)
    for y_t in trouble_ys:
        # get all positions of the troubling y
        y_positions = [j for j, y in enumerate(ys) if y == y_t]
        # then get all x's which are mapped to it
        trouble_xs.append([ref_xs[j] for j in y_positions])

    for x_pair in trouble_xs:
        if len(x_pair) != 2:
            # Check if we are dealing with a set of 
            # sheets which are all (approximately) 0.0
            if max(map(abs, x_pair)) < accuracy:
                for x_pair_i in x_pair:
                    correct_xy_pairs.update({x_pair_i: x_pair_i})
            else:
                raise Exception('Cannot handle this kind of sheet degeneracy')
        else:
            closest_ys_0 = nsmallest(
                2, new_xs, key=lambda x: abs(x - x_pair[0])
            )
            closest_ys_1 = nsmallest(
                2, new_xs, key=lambda x: abs(x - x_pair[1])
            )
            # a check
            if (
                closest_ys_0 != closest_ys_1 and 
                closest_ys_0 != closest_ys_1.reverse()
            ):
                logger.warning((
                    'the closest sheets to the reference pair {}'
                    '\ndont match: they are respectively:\n{}\n{}'
                ).format(x_pair, closest_ys_0, closest_ys_1))
            
            # compute the differences of the various combinations
            dx_00 = closest_ys_0[0] - x_pair[0]
            dx_01 = closest_ys_0[1] - x_pair[0]
            dx_10 = closest_ys_1[0] - x_pair[1]
            dx_11 = closest_ys_1[1] - x_pair[1]
            # pick for each x in the x_pair its companion based on 
            # the phase of the displacement, choosing the closest to 
            # the previous step in the tracking
            i_0 = ref_xs.index(x_pair[0])
            i_1 = ref_xs.index(x_pair[1])
            ref_dx_0 = delta_xs[i_0]
            ref_dx_1 = delta_xs[i_1]
            # first find the companion for x_pair[0]
            if abs(phase(dx_00 / ref_dx_0)) < abs(phase(dx_01 / ref_dx_0)):
                correct_xy_pairs.update({x_pair[0]: closest_ys_0[0]})
            else:
                correct_xy_pairs.update({x_pair[0]: closest_ys_0[1]})
            # then repeat for x_pair[1]
            if abs(phase(dx_10 / ref_dx_1)) < abs(phase(dx_11 / ref_dx_1)):
                correct_xy_pairs.update({x_pair[1]: closest_ys_0[0]})
            else:
                correct_xy_pairs.update({x_pair[1]: closest_ys_0[1]})

    # at this point, we should have sorted all the new_xs
    # we check if the sorting was successful

    sorted_xs = [correct_xy_pairs[x] for x in ref_xs]
    unique_sorted_xs = n_remove_duplicate(sorted_xs, 0.0)

    if len(sorted_xs) == len(unique_sorted_xs):
        return sorted_xs
    else:
        return 'sorting failed'


def kr_delta(i, j):
    if i == j:
        return 1
    else:
        return 0


def get_positive_roots_of_branch_point(bp, g_data, logger_name='loom'):
    """
    Determines the positive roots associated with 
    a branch point's 'structure', i.e. how the sheets
    collide at the branch point.
    It will return a minimal list, i.e. it will drop
    any redundant roots that can be obtained as linear
    combinations of others.
    """
    logger = logging.getLogger(logger_name)
    vanishing_positive_roots = []
    positive_roots = g_data.positive_roots
    # Note that bp.groups carries indices, which can be used
    # to map each x at the reference point to the weights, i.e.
    # reference_xs[i] <-> weights[i].
    weights = g_data.weights

    for g in bp.groups:
        # Within each group of colliding sheets/weights,
        # consider all possible pairs, and compute 
        # the corresponding difference.
        # Then add it to the vanishing positive roots.
        for s_1, s_2 in combinations(g, 2):
            v_1 = weights[s_1]
            v_2 = weights[s_2]
            if any(numpy.allclose(v_1 - v_2, x) for x in positive_roots):
                vanishing_positive_roots.append(v_1 - v_2)

            elif any(numpy.allclose(v_2 - v_1, x) for x in positive_roots):
                vanishing_positive_roots.append(v_2 - v_1)

            else:
                continue
    if vanishing_positive_roots == []:
        # TODO: Check if the situation is what we expect and is normal.
        # If that's the case, don't print a message. 
        # Otherwise print a warning instead of info.
        logger.debug(
            "Branch point doesn't correspond "
            "to a positive root. May be an accidental branch point."
        )
        return []

    # Finally, cleanup the duplicates, 
    # as well as the roots which are not linearly independent
    return numpy.array(
        keep_linearly_independent_vectors(vanishing_positive_roots)
    )


def belongs_to_cluster(x, c, enum_sh):
    """
    Given a cluster of sheets, c = [i_0, i_1, ...]
    specified by means of their integer labels,
    it determines whether a sheet with coordinate 'x'
    is close enough to ANY of the sheets in 'c'
    to be considered as part of it.
    The positions of sheets in the cluster are extracted
    from enum_sh = [...[i_k, x_k]...]
    """
    test = False
    for i in c:
        # pick the coordinate of the sheet with label 'i'
        y_i = [y for j, y in enum_sh if j == i][0]
        if abs(y_i - x) < BP_PROXIMITY_THRESHOLD:
            test = True
            break

    if test is False:
        return False
    if test is True:
        return True


def keep_linearly_independent_vectors(vector_list):
    """
    Takes a list of numpy arrays and returns a 
    subset of linearly independent ones.
    """
    
    first_vector = vector_list[0]
    independent_list = [first_vector]

    m_rank = 1
    m = numpy.matrix([first_vector])
    for v in vector_list:
        # add the vector as a row to the matrix, 
        # then compute the rank
        new_m = numpy.vstack([m, v])
        new_m_rank = matrix_rank(new_m)
        if new_m_rank > m_rank:
            m = new_m
            m_rank = new_m_rank
            independent_list.append(v)

    return independent_list


def build_monodromy_matrix(initial_sheets, sorted_sheets, logger_name='loom'):
    # Now we have three lists:
    # initial_sheets = [[0, x_0], [1, x_1], ...]
    # final_sheets = [[0, x'_0], [1, x'_1], ...]
    # sorted_sheets = [[i_0, x_0], [i_1, x_1], ...]
    # therefore the monodromy permutation corresponds
    # to 0 -> i_0, 1 -> i_1, etc.
    logger = logging.getLogger(logger_name)
    n_sheets = len(initial_sheets)

    # NOTE: in the following basis vectors, i = 0 , ... , n-1
    def basis_e(i):
        return numpy.array([kr_delta(j, i) for j in range(n_sheets)])

    perm_list = []
    for i in range(n_sheets):
        new_sheet_index = sorted_sheets[i][0]
        perm_list.append(basis_e(new_sheet_index))

    # perm_matrix = numpy.array(perm_list).transpose()
    perm_matrix = numpy.array(perm_list)
    logger.debug('Permutation matrix {}'.format(perm_matrix))
    return perm_matrix

