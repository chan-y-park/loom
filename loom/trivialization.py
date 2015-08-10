import sympy
import matplotlib.pyplot as plt
import cmath
import numpy as np

from sympy import Poly
from cmath import exp, pi
from numpy.linalg import matrix_rank
from itertools import combinations

from geometry import SWData
from misc import delete_duplicates

### number of steps used to track the sheets along a leg 
### the path used to trivialize the cover at any given point
N_PATH_TO_PT = 100

### number of steps for each SEGMENT of the path around a 
### branching point (either branch-point, or irregular singularity)
N_PATH_AROUND_PT = 30

### Tolerance for recognizing colliding sheets at a branch-point
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

    path_to_bp :
        A path running from the basepoint of the trivialization
        to the branch point without crossing any branch cut.
    
    sheet_tracks_to_bp :
        A list of sheet tracks, i.e. the x-values of each
        sheet as it is tracked along a path that runs to
        the branch point, to determine collision structure 
        of the various sheets.
    
    positive_roots :
        A minimal list of positive roots characterizing the 
        groups of colliding sheets at the branch point.
    
    path_around_bp :
        A path encircling the branch point and no one else,
        used to compute the monodromy.
    
    sheet_tracks_around_bp :
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

    """
    def __init__(self, z=None):
        self.z = z
        self.groups = None 
        self.singles = None
        self.enum_sh = None
        #self.path_to_z = None
        self.positive_roots = None 
        self.order = None
        self.monodromy = None

        #self.data = trivialization.analyze_branch_point(self.z)
        #self.trivialization = trivialization
        
        #bp_data = trivialization.analyze_branch_point(self.z)
        #self.sheet_tracks_to_bp = bp_data['tracked_sheets']
        #self.path_to_bp = bp_data['path_to_branch_point']

        #self.sheet_tracks_around_bp = (
        #    trivialization.track_sheets_along_path(path_around_z)
        #)

        #path_around_z = trivialization.path_around_pt(z)
        #self.monodromy = trivialization.sheet_monodromy(path_around_z)

    def print_info(self):
        print "\n---------------------------------------------------------\
               \nBranch Point at z = %s\
               \n---------------------------------------------------------"\
               % self.z
        print "this is the branch point structure"
        print "groups = {}".format(self.groups)
        print "singles = {}".format(self.singles)
        print "positive roots = {}".format(self.positive_roots)
        print "order = {}".format(self.order)
        print "sheets at the branch point = {}".format(self.enum_sh)
        print "sheet monodromy permutation matrix = \n{}".format(self.monodromy)        

class IrregularSingularity:
    """
    The IrregularSingularity class.
    Just a container of information.
    Strictly related to the first fundamental representation cover.
    """
    def __init__(self, z=None):
        self.z = z
        self.monodromy = None
        
        #path_around_z = trivialization.path_around(z)
        #self.sheet_tracks_around_irr_sing = (
        #    trivialization.track_sheets_along_path(path_around_z)
        #)
        #self.monodromy = (
        #    trivialization.sheet_monodromy(path_around_z)
        #)
        
    def print_info(self):
        print "\n---------------------------------------------------------\
               \nIrregular singularity at z = %s\
               \n---------------------------------------------------------"\
               % self.z
        print "sheet monodromy permutation matrix = \n{}".format(self.monodromy)        


### TODO: Use g_data.weights at the base point as labels of sheets,
### instead of integer indicies. Just a conceptual issue, because
### the integer indicies are labeling g_data.weights in the same order.
class SWDataWithTrivialization(SWData):
    """
    All branch cuts are assumed to run vertically, emanating
    upwards from branch points and irregular singularities.

    Arguments
    ---------
    
    sw_data : 
        an object of the type SWData, whose attribute 'curve' 
        should correspond to the curve in the FIRST fundamental
        representation of the Lie algebra
    
    ramification_points : 
        a list of objects of the type RamificationPoint, corresponding
        to the given sw_data.

    Attributes & Methods
    --------------------

    base_point : 
        the base point of the trivialization

    reference_sheets :
        a list x's 
            [x_0, x_1, ..., x_i, ...]
        where 'i' is an integer label for the sheet,
        and 'x' is its position in the fiber of T^*C 
        over the basepoint. This is aligned with 
        g_data.ffr_weights.

    sheets_at_z(z) :
        this method returns the set of sheets and their integer label 
        identifier at any point 'z' on the C-plane.
        The labels are consistent with those at the basepoint.
        To get the corresponding weights, of the firt fundamental 
        representation, use g_data.weights[i].
        The output looks like this
        {0 : x_0, ... , i : x_i, ...}
    """
    ### NOTE: I am assuming that branch points do not overlap vertically
    ### this should be guaranteed by introducing an automatic rotation of 
    ### the z-plane before calling this class.
    ### NOTE: I am restricting to square-root type branch points.
    ### Although I am not printing any explicit warning/error message 
    ### and the computation will go through for higher-type, but give a 
    ### wrong answer!
    def __init__(self, config,):
        super(SWDataWithTrivialization, self).__init__(config)

        self.branch_points = []
        self.irregular_singularities = []

        # z-coords of branch points.
        bpzs = delete_duplicates(
            [r.z for r in self.ramification_points if not r.is_puncture]
        )
        # z-coords of irregular singularities.
        iszs = delete_duplicates(
            [r.z for r in self.ramification_points if r.is_puncture]
        )
        
        ### Automatically choose a basepoint, based on the positions of
        ### both branch points and irregular singularities
        all_points_z = bpzs + iszs
        all_distances = [abs(x - y) for x in all_points_z
                         for y in all_points_z]
        max_distance = max(all_distances)
        center = sum([z_pt for z_pt in all_points_z]) / len(all_points_z)
        self.base_point = center - 1j * max_distance
        
        ### Minimun distance between the base point and 
        ### branch points/punctures.
        non_zero_distances = [x for x in all_distances if x!=0.0]
        self.min_distance = min(non_zero_distances)

        ### Fix reference x's at the basepoints.
        ### These sheets are aligned in the order of
        ### sw.g_data.weights, i.e. reference_sheets[i]
        ### is the value of x corresponding to 
        ### sw.g_data.weights[i].
        self.reference_ffr_xs, self.reference_xs = self.get_aligned_xs(
            self.base_point,
        )
        #self.reference_sheets = {i: x for i, x in enumerate(self.reference_xs)}
        
        ### XXX: sheet_weight_dictionary -> sw.g_data.weights
        #self.sheet_weight_dictionary = self.build_dictionary()

        ### Construct the list of branch points
        for z_bp in bpzs:
            bp = BranchPoint(z=z_bp)
            self.analyze_branch_point(bp)
            self.branch_points.append(bp)

        ### Construct the list of irregular singularities
        for z_irr_sing in iszs:
            irr_sing = IrregularSingularity(z=z_irr_sing)
            self.analyze_irregular_singularity(irr_sing)
            self.irregular_singularities.append(irr_sing)

        
    ### TODO: Need to implement tracking without using aligned x's?
    def get_sheets_along_path(self, z_path, is_path_to_bp=False,
                              from_ffr=True):
        """
        Tracks the sheets along a path.
        It checks at each step that tracking is successful,
        meaning that all sheets can be distinguished correctly.
        This would fail if we choose a path ending on a branch-point.
        For tracking roots as we run into a branch point, one should
        set the variable 'is_path_to_bp=True', and the check for 
        sheets becoming too similar will be ignored altogether.
        """
        ### NOTE: instead of disabling the check, we could 
        ### still perform it, and suppress the check only towards the 
        ### end of the tracking.

        g_data = self.g_data
        ffr_xs_0 = self.reference_ffr_xs
        xs_0 = self.reference_xs
        ### Each element is a sheet, which is a list of x's along the path.
        ### Initialized with reference_xs.
        ### TODO: set each element to an integer rather than a float.
        sheets_along_path = [[x] for x in xs_0]
        
        for i, z in enumerate(z_path):
            #xs_1 = self.get_xs(z)
            ffr_xs_1, xs_1 = self.get_aligned_xs(z)
            if is_path_to_bp == False:
                sorted_ffr_xs = get_sorted_xs(
                    ffr_xs_0, ffr_xs_1, check_tracking=True, 
                    index=1, z_0=z_path[i-1], z_1=z_path[i]
                )
            else:
                sorted_ffr_xs = get_sorted_xs(ffr_xs_0, ffr_xs_1,
                                              check_tracking=False)
            if g_data.fundamental_representation_index == 1:
                sorted_xs = sorted_ffr_xs
            else:
                sorted_xs = self.get_xs_of_weights_from_ffr_xs(sorted_ffr_xs)
            for j, s_j in enumerate(sheets_along_path):
                s_j.append(sorted_xs[j])
            ffr_xs_0 = sorted_ffr_xs

        ### the result is of the form [sheet_path_1, sheet_path_2, ...]
        ### where sheet_path_i = [x_0, x_1, ...] are the fiber coordinates
        ### of the sheet along the path
        return sheets_along_path


    def get_sheets_at_z(self, z_pt, g_data=None):
        """
        Returns a dict of (sheet_index, x) at a point ''z_pt'', 
        which cannot be a branch point or a singularity.
        """
        z_path = get_path_to(z_pt, self.base_point)
        sheets = self.get_sheets_along_path(z_path)
        final_xs = [s_i[-1] for s_i in sheets]
        final_sheets = {i : x for i, x in enumerate(final_x)}
        return final_sheets

    
    ### TODO: Review this method.
    def get_sheet_monodromy(self, z_path):
        """
        Compares the x-coordinates of sheets at the 
        beginning and at the end of a CLOSED path.
        Returns a permutation matrix, expressed in 
        the basis of reference sheets, such that
        new_sheets = M . old_sheets
        """
        initial_xs = self.reference_xs
        initial_sheets = [[i, x] for i, x in enumerate(initial_xs)]
        final_xs = [sheet_i[-1] 
                    for sheet_i in self.get_sheets_along_path(z_path)]
        final_sheets = [[i, x] for i, x in enumerate(final_xs)]

        ### Now we compare the initial and final sheets 
        ### to extract the monodromy permutation
        ### recall that each entry of initial_sheets and final_sheets
        ### is of the form [i, x] with i the integer label
        ### and x the actual position of the sheet in the fiber 
        ### above the basepoint.
        sorted_sheets = []
        for s_1 in initial_sheets:
            closest_candidate = final_sheets[0]
            min_d = abs(s_1[1] - closest_candidate[1])
            for s_2 in final_sheets:
                if abs(s_2[1] - s_1[1]) < min_d:
                    min_d = abs(s_2[1] - s_1[1])
                    closest_candidate = s_2
            sorted_sheets.append(closest_candidate)
        
        ### Now we check that sheet tracking is not making a mistake.
        ### NOTE: cannot use the function 'delete_duplicates' with this 
        ### data structure.
        seen = set()
        uniq = []
        for s in sorted_sheets:
            if s[1] not in seen:
                uniq.append(s[1])
                seen.add(s[1])
        if len(uniq) < len(sorted_sheets):
            raise ValueError(
                '\nError in determination of monodromy!\n' +
                'Cannot match uniquely the initial sheets to the final ones.'
            )
        else:
            pass

        ### Now we have tree lists:
        ### initial_sheets = [[0, x_0], [1, x_1], ...]
        ### final_sheets = [[0, x'_0], [1, x'_1], ...]
        ### sorted_sheets = [[i_0, x_0], [i_1, x_1], ...]
        ### therefore the monodromy permutation corresponds
        ### to 0 -> i_0, 1 -> i_1, etc.

        n_sheets = len(initial_sheets)
        
        ### NOTE: in the following basis vectors, i = 0 , ... , n-1
        def basis_e(i):
            return np.array([kr_delta(j, i) for j in range(n_sheets)])

        perm_list = []
        for i in range(n_sheets):
            new_sheet_index = sorted_sheets[i][0]
            perm_list.append(basis_e(new_sheet_index))

        perm_matrix = np.matrix(perm_list).transpose()

        return perm_matrix


    def analyze_branch_point(self, bp):
        g_data = self.g_data
        path_to_bp = get_path_to(bp.z, self.base_point)
        sheets_along_path = self.get_sheets_along_path(
            path_to_bp, is_path_to_bp=True
        )
        xs_at_bp = [s_i[-1] for s_i in sheets_along_path]
        enum_sh = [[i, x_i] for i, x_i in enumerate(xs_at_bp)]
        
        clusters = []
        for i, x in enum_sh:
            is_single = True
            for c_index, c in enumerate(clusters):
                x_belongs_to_c = belongs_to_cluster(x, c, enum_sh)
                if x_belongs_to_c == True:
                    clusters[c_index].append(i)
                    is_single = False
                    break
            if is_single == True:
                clusters.append([i])

        bp.enum_sh = enum_sh
        bp.groups = [c for c in clusters if len(c) > 1]
        bp.singles = [c[0] for c in clusters if len(c) == 1]

        bp.positive_roots = get_positive_roots_of_branch_point(
            bp, self.g_data,  
        )
        bp.order = len(bp.positive_roots) + 1

        path_around_bp = get_path_around(bp.z, self.base_point,
                                         self.min_distance)
        bp.monodromy = self.get_sheet_monodromy(path_around_bp)


    def analyze_irregular_singularity(self, irr_sing):
        path_around_z = get_path_around(z)
        self.monodromy = (
            self.get_sheet_monodromy(path_around_z)
        )
        

def get_path_to(z_pt, base_pt):
    """
    Return a rectangular path from the base point to z_pt.
    """
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


def get_path_around(z_pt, base_pt, min_distance):
    z_0 = base_pt
    z_1 = 1j * base_pt.imag + z_pt.real
    radius = min_distance / 2.0
    z_2 = z_pt - 1j * radius

    steps = N_PATH_AROUND_PT
    path_segment_1 = [z_0 + ((z_1 - z_0) / steps) * i
                      for i in range(steps + 1)]
    path_segment_2 = [z_1 + ((z_2 - z_1) / steps) * i 
                      for i in range(steps + 1)]
    path_segment_3 = [z_pt + radius * (-1j) * exp(i * 2.0 * pi * 1j
                                                  / steps) 
                      for i in range(steps +1)]
    path_segment_4 = path_segment_2[::-1]
    path_segment_5 = path_segment_1[::-1]

    return (path_segment_1 + path_segment_2 + path_segment_3 +
            path_segment_4 + path_segment_5)


### TODO: Try using numba.
def get_sorted_xs(ref_xs, new_xs, check_tracking=True, 
                  index=None, z_0=None, z_1=None):
    """
    Returns a sorted version of 'new_xs'
    based on matching the closest points with 
    'ref_xs'
    """
    sorted_xs = []
    for s_1 in ref_xs:
        closest_candidate = new_xs[0]
        min_d = abs(s_1 - closest_candidate)
        for s_2 in new_xs:
            if abs(s_2 - s_1) < min_d:
                min_d = abs(s_2 - s_1)
                closest_candidate = s_2
        sorted_xs.append(closest_candidate)
    
    if check_tracking == True:
        ### Now we check that sheet tracking is not making a mistake.
        unique_sorted_xs = delete_duplicates(sorted_xs)
        if len(unique_sorted_xs) < len(sorted_xs):
            print "\nAt step %s, between %s and %s " % (index, z_0, z_1)
            print "old xs" 
            print ref_xs
            print "new xs"
            print new_xs
            raise ValueError(
                '\nCannot track the sheets!\n'
                'Probably passing too close to a branch point.'
            )
        else:
            return sorted_xs
    else:
        ### If the path is one ending on a branch-point, 
        ### the check that tracking is correct is disabled
        ### because it would produce an error, since by definition
        ### sheets will be indistinguishable at the very end.
        return sorted_xs


def kr_delta(i, j):
    if i == j:
        return 1
    else:
        return 0


def get_positive_roots_of_branch_point(bp, g_data):
    """
    Determines the positive roots associated with 
    a branch point's 'structure', i.e. how the sheets
    collide at the branch point.
    It will return a minimal list, i.e. it will drop
    any redundant roots that can be obtained as linear
    combinations of others.
    """
    vanishing_positive_roots = []
    positive_roots = g_data.positive_roots
    ### Note that bp.groups carries indicies, which can be used
    ### to map each x at the reference point to the weights, i.e.
    ### reference_xs[i] <-> weights[i].
    weights = g_data.weights


    for g in bp.groups:
        ### Within each group of colliding sheets/weights,
        ### consider all possible pairs, and compute 
        ### the corresponding difference.
        ### Then add it to the vanishing positive roots.
        for s_1, s_2 in combinations(g, 2):
            v_1 = weights[s_1]
            v_2 = weights[s_2]
            if any(np.allclose(v_1 - v_2, x) for x in positive_roots):
                vanishing_positive_roots.append(v_1 - v_2)

            elif any(np.allclose(v_2 - v_1, x) for x in positive_roots):
                vanishing_positive_roots.append(v_2 - v_1)

            else:
                raise ValueError("Branch point doesn't correspond "
                                 "to a positive root.")

    ### Finally, cleanup the duplicates, 
    ### as well as the roots which are not linearly independent
    ### TODO: Check if we really need to remove linearly depedent 
    ### roots. Isn't it part of the information a branch pt carries?
    return keep_linearly_independent_vectors(vanishing_positive_roots)


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
        ### pick the coordinate of the sheet with label 'i'
        y_i = [y for j, y in enum_sh if j==i][0]
        if abs(y_i - x) < BP_PROXIMITY_THRESHOLD:
            test = True
            break

    if test == False:
        return False
    if test == True:
        return True


def keep_linearly_independent_vectors(vector_list):
    """
    Takes a list of numpy arrays and returns a 
    subset of linearly independent ones.
    """
    
    first_vector = vector_list[0]
    independent_list = [first_vector]

    m_rank = 1
    m = np.matrix([first_vector])
    for v in vector_list:
        ### add the vector as a row to the matrix, 
        ### then compute the rank
        new_m = np.vstack([m,v])
        new_m_rank = matrix_rank(new_m)
        if new_m_rank > m_rank:
            m = new_m
            m_rank = new_m_rank
            independent_list.append(v)

    return independent_list
