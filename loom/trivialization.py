from api import load_config
from geometry import SWData, get_ramification_points, RamificationPoint
import sympy
import matplotlib.pyplot as plt
import cmath
import numpy as np
from sympy import Poly
from cmath import exp, pi

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
    Just a container of information.
    """
    def __init__(self, z=None, trivialization=None):
        self.z = z
        self.trivialization = trivialization
        
        bp_data = self.trivialization.analyze_branch_point(self.z)
        self.pairs = bp_data['pairs']
        self.singles = bp_data['singles']
        self.enum_sh = bp_data['enum_sh']
        self.sheets_at_branch_point = bp_data['sheets_at_branch_point']
        self.sheet_tracks_to_bp = bp_data['tracked_sheets']
        self.path_to_bp = bp_data['path_to_branch_point']
                        
        self.positive_root = trivialization.positive_root(self.pairs, self.singles, self.enum_sh)
        
        self.path_around_bp = self.trivialization.path_around_pt(self.z)
        self.sheet_tracks_around_bp = self.trivialization.track_sheets_along_path(self.path_around_bp)
        self.monodromy = trivialization.sheet_monodromy(self.path_around_bp)


class IrregularSingularity:
    """
    The IrregularSingularity class.
    Just a container of information.
    """
    def __init__(self, z=None, trivialization=None):
        self.z = z
        self.trivialization = trivialization
        
        self.path_around_irr_sing = self.trivialization.path_around_pt(self.z)
        self.sheet_tracks_around_irr_sing = self.trivialization.track_sheets_along_path(self.path_around_irr_sing)
        self.monodromy = trivialization.sheet_monodromy(self.path_around_irr_sing)
        
    


class Trivialization:
    """
    The Trivialization class.

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

    lie_algebra : 
    the Lie algebra associated with the cover, expressed
    as a list of a capital letter and the rank, e.g.
    ['A', 3] for the A_3 algebra.


    Attributes & Methods
    --------------------

    basepoint : 
    the base point of the trivialization

    reference_sheets :
    a list of pairs [i, x] where 'i' is an integer label
    for the sheet, and 'x' is its position in the fiber of T^*C 
    over the basepoint

    sheet_weight_dictionary :
    a dictionary between the sheet integer labels and the
    weights of the FIRST fundamental representation
    it is structured as follows
    {i_0 : v_0, ... , i_k : v_k , ...}
    where 'v_k' are numpy arrays corresponding to weights.
    - For g=A_n Lie algebras, the weights are given in IR^{n+1}
        v_0 = (1,0,...,0)
        v_1 = (0,1,0,..) 
        ...
        v_n = (0,...,0,1) 
      In this case, it does not matter how we identify weights
      with sheets, since the Weyl group acts by permuting all of 
      them freely.
    - For g=D_n, the weights are given in IR^{n}
        v_0 = (1,0,...,0)
        v_1 = (0,1,...,0)
        v_{n-1} = (0,...,0,1)
        v_n = (-1,0,...,0)
        v_{n+1} = (0,-1,...,0)
        v_{2n-1} = (0,...,0,-1)
      In this case, we diivde the sheets into positive and negative ones,
      and assign the weights accordingly.
      The assignment of positive sheets is almost arbitrary: from each pair
      of positive/negative sheets one can pick either, as long as one makes
      an even number of "sign mistakes". We don't keep track of this,
      as a result there is an ambiguity in distinguishing one spinor 
      representation from the other

    sheets_at_arbitrary_z(z) :
    this method returns the set of sheets and their integer label 
    identifier at any point 'z' on the C-plane.
    These are the sheets of the FIRST FUNDAMENTAL representation.
    The labels are consistent with those at the basepoint.
    To get the corresponding weights, of the firt fundamental 
    representation, the dictionary should be invoked.
    The output looks like this
    [[0, x_0] ... [i, x_i] ...]

    branch_points :
    A list of all the branch points.

    irregular_singularities :
    A list of all the irregular singularities.

    """
    ### NOTE: I am assuming that branch points do not overlap vertically
    ### this should be guaranteed by introducing an automatic rotation of 
    ### the z-plane before calling this class.
    ### NOTE: I am restricting to square-root type branch points.
    ### Although I am not printing any explicit warning/error message 
    ### and the computation will go through for higher-type, but give a wrong answer!
    def __init__(self, sw_data, ramification_points, lie_algebra):
        self.sw_data = sw_data
        self.algebra = lie_algebra
        
        self.branch_points = []
        self.irregular_singularities = []
        self.reference_sheets = None
        self.sheet_weight_dictionary = None
        self.basepoint = None
        self.min_distance = None
        self.max_distance = None
        self.center = None

        b_points_z = bp_from_ramif(ramification_points)
        irr_sing_z = irr_sing_from_ramif(ramification_points)

        print "\nbranch points"
        print b_points_z
        print "\nirregular singularities"
        print irr_sing_z
        
        ### Automatically choose a basepoint, based on the positions of
        ### both branch points and irregular singularities
        all_points_z = b_points_z + irr_sing_z
        all_distances = [abs(x - y) for x in all_points_z
                                                        for y in all_points_z]
        self.max_distance = max(all_distances)
        non_zero_distances = [x for x in all_distances if x!=0.0]
        self.min_distance = min(non_zero_distances)
        self.center = sum([z_pt for z_pt in all_points_z]) / len(all_points_z)
        self.basepoint = self.center - 1j * self.max_distance

        ### Fix reference sheets at the basepoints, i.e. assign an integer 
        ### label to each sheet
        self.reference_sheets = [[i, x] \
                    for i, x in enumerate(self.sheets_at_z(self.basepoint))]
        
        ### Now that we have determined the reference sheets, 
        ### we build a weight-sheet dictionary
        self.sheet_weight_dictionary = self.build_dictionary()

        ### Construct the list of branch points
        for i, z_bp in enumerate(b_points_z):
            self.branch_points.append(BranchPoint(
                                                    z=z_bp, 
                                                    trivialization=self
                                                ))

        ### Construct the list of irregular singularities
        for z_irr_sing in irr_sing_z:
            self.irregular_singularities.append(IrregularSingularity(
                                                    z=z_irr_sing, 
                                                    trivialization=self
                                                ))

        ### PRESENT THE DATA FOR DEBUGGING PURPOSES
        print "\nSheets of the cover at the baspoint z_0 = {}".format(self.basepoint)
        print self.reference_sheets

        print "\nThe dictionary between sheets and weights:"
        print self.sheet_weight_dictionary

        for i, bp in enumerate(self.branch_points):
            if SHOW_TRACKING_PLOTS == True:
                print "\nroot tracking along the path to branch point #{}".format(i)
                for j, sheet_list in enumerate(bp.sheet_tracks_to_bp):
                    data_plot(sheet_list, 'sheet {} tracked to branch point {}'.format(j, i))
                for j, sheet_list in enumerate(bp.sheet_tracks_around_bp):
                    data_plot(sheet_list, 'sheet {} tracked around branch point {}'.format(j, i))
            print "this is the branch point structure for branch point #{}".format(i)
            print "pairs = {}".format(bp.pairs)
            print "singles = {}".format(bp.singles)
            print "positive root = {}".format(bp.positive_root)
            print "sheets at the branch point = {}".format(bp.sheets_at_branch_point)
            print "sheet monodromy permutation matrix = \n{}".format(bp.monodromy)


        
    def sheets_at_z(self, z_0):
        from sympy.abc import x, z
        sw_curve_fiber = self.sw_data.curve.num_eq.subs(z, z_0)
        sym_poly = Poly(sw_curve_fiber, x, domain='CC')
        coeff_list = map(complex, sym_poly.all_coeffs())
        return map(complex, np.roots(coeff_list))


    def path_to_pt(self, z_pt):
        z_0 = self.basepoint
        z_1 = 1j * self.basepoint.imag + z_pt.real
        z_2 = z_pt
        half_steps = int(N_PATH_TO_PT / 2)
        return [z_0 + ((z_1 - z_0) / half_steps) * i \
                                        for i in range(half_steps + 1)] \
            + [z_1 + ((z_2 - z_1) / half_steps) * i \
                                        for i in range(half_steps + 1)]                


    def path_around_pt(self, z_pt):
        z_0 = self.basepoint
        z_1 = 1j * self.basepoint.imag + z_pt.real
        radius = self.min_distance / 2.0
        z_2 = z_pt - 1j * radius

        steps = N_PATH_AROUND_PT
        path_segment_1 = [z_0 + ((z_1 - z_0) / steps) * i \
                                        for i in range(steps + 1)]
        path_segment_2 = [z_1 + ((z_2 - z_1) / steps) * i \
                                        for i in range(steps + 1)]
        path_segment_3 = [z_pt + radius * (-1j) * \
                                exp(i * 2.0 * pi * 1j / steps) \
                                for i in range(steps +1)]
        path_segment_4 = path_segment_2[::-1]
        path_segment_5 = path_segment_1[::-1]

        return path_segment_1 + path_segment_2 + path_segment_3 \
                + path_segment_4 + path_segment_5


    def track_sheets_along_path(self, z_path, is_path_to_bp=False):
        """
        Tracks the sheets along a path.
        It checks at each step that tracking is successful,
        meaning that all sheets can be distinguished correctly.
        This would fail if we choose a path ending on a branch-point.
        For tracking roots as we run into a branch point, one should
        set the variable 'is_path_to_bp=True', and the check for 
        sheets becoming too similar will be ignored altogether.
        """
        ### NOTE: instead of disabling the chack, we could 
        ### still perform it, and suppress the check only towards the 
        ### end of the tracking.

        sheets_0 = [x for i, x in self.reference_sheets]
        sheets_along_path = [[s] for s in sheets_0]
        
        if is_path_to_bp == False:
            for i, z in enumerate(z_path):
                sheets_1 = self.sheets_at_z(z)
                sheets_0 = self.sort_sheets(sheets_0, sheets_1, check_tracking=True, index=1, z_0=z_path[i-1], z_1=z_path[i])
                for i, s_list in enumerate(sheets_along_path):
                    s_list.append(sheets_0[i])

            return sheets_along_path
            ### the result is of the form [sheet_path_1, sheet_path_2, ...]
            ### where sheet_path_i = [x_0, x_1, ...] are the fiber coordinates
            ### of the sheet along the path
        else:
            for i, z in enumerate(z_path):
                sheets_1 = self.sheets_at_z(z)
                sheets_0 = self.sort_sheets(sheets_0, sheets_1, check_tracking=False)
                for i, s_list in enumerate(sheets_along_path):
                    s_list.append(sheets_0[i])

            return sheets_along_path


    def sheets_at_arbitrary_z(self, z_pt):
        """
        Returns a list of sheets of the fuandamental cover
        with their identifiers at any point 'z'.
        The choice of 'z' cannot be a branch point or a singularity.
        The weights of the fundamental cover, corresponding to the sheets,
        can be obtained by invoking the dictionary of the trivialization
        """

        z_path = self.path_to_pt(z_pt)
        sheet_tracks = self.track_sheets_along_path(z_path)
        final_x = [sheet_list[-1] for sheet_list in sheet_tracks]
        final_sheets = [[i, x] for i, x in enumerate(final_x)]
        return final_sheets



    def sort_sheets(self, ref_sheets, new_sheets, check_tracking=True, index=None, z_0=None, z_1=None):
        """
        Returns a sorted version of 'new_sheets'
        based on matching the closest points with 
        'ref_sheets'
        """
        sorted_sheets = []
        for s_1 in ref_sheets:
            closest_candidate = new_sheets[0]
            min_d = abs(s_1 - closest_candidate)
            for s_2 in new_sheets:
                if abs(s_2 - s_1) < min_d:
                    min_d = abs(s_2 - s_1)
                    closest_candidate = s_2
            sorted_sheets.append(closest_candidate)
        
        if check_tracking == True:
            ### Now we check that sheet tracking is not making 
            ### a mistake.
            seen = set()
            uniq = []
            for x in sorted_sheets:
                if x not in seen:
                    uniq.append(x)
                    seen.add(x)
            if len(uniq) < len(sorted_sheets):
                print "\nAt step %s, between %s and %s " % (index, z_0, z_1)
                print "old sheets" 
                print ref_sheets
                print "new sheets"
                print new_sheets
                raise ValueError('\nCannot track the sheets!\n'+\
                        'Probably passing too close to a branch point.')
            else:
                return sorted_sheets
        else:
            ### If the path is one ending on a branch-point, 
            ### the check that tracking is correct is disabled
            ### because it would produce an error, since by definition
            ### sheets will be indistinguishable at the very end.
            return sorted_sheets


    def analyze_branch_point(self, z_bp):
        z_bp_path = self.path_to_pt(z_bp)
        tracked_sheets = self.track_sheets_along_path(z_bp_path, is_path_to_bp=True)
        sheets_at_bp = [sheet_list[-1] for sheet_list in tracked_sheets]
        enum_sh = [[i, s_i] for i, s_i in enumerate(sheets_at_bp)]
        
        pairs = []
        singles = []

        for i, x in enum_sh:
            if i in flatten(pairs):
                pass
            elif i == len(enum_sh)-1:
                ### this is the last sheet of the list
                ### if it's not already in a pair, then it's a single
                singles.append(i)
            else:
                paired = False
                for j, y in enum_sh[i+1:]:                    
                    if abs(x - y) < BP_PROXIMITY_THRESHOLD:
                        paired = True
                        pairs.append([i, j])
                ### NOTE: we can in principle have multiple pairings, meaning
                ### that three or more sheets could collide together
                ### these will show up as several pairs containing the 
                ### same numbers e.g. [i, j], [j, k], [k, i] would mean
                ### that sheets i, j, k collide all together
                ### Should introduce a check that handles this situations!

                if paired == False:
                    singles.append(i)

        return {'pairs' : pairs, \
                'singles' : singles, \
                'enum_sh' : enum_sh, \
                'tracked_sheets' : tracked_sheets, \
                'sheets_at_branch_point' : sheets_at_bp, \
                'path_to_branch_point' : z_bp_path
                }


    def sheet_monodromy(self, z_path):
        """
        Compares the x-coordinates of sheets at the 
        beginning and at the end of a CLOSED path.
        Returns a permutation matrix, expressed in 
        the basis of reference sheets, such that
        new_sheets = M . old_sheets
        """

        initial_sheets = self.reference_sheets
        final_x = [sheet_list[-1] \
                        for sheet_list in self.track_sheets_along_path(z_path)]
        final_sheets = [[i, x] for i, x in enumerate(final_x)]

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
        
        ### Now we check that sheet tracking is not making 
        ### a mistake.
        seen = set()
        uniq = []
        for s in sorted_sheets:
            if s[1] not in seen:
                uniq.append(s[1])
                seen.add(s[1])
        if len(uniq) < len(sorted_sheets):
            raise ValueError('\nError in determination of monodromy!\n'+\
                'Cannot match uniquely the initial sheets to the final ones.')
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


    def build_dictionary(self):
        algebra = self.algebra
        r = algebra[1]

        if algebra[0] == 'A':
            ### for example, fund_weights(2) will be [0, 1, 0, ...]
            def fund_weights(i):
                return np.array([kr_delta(j, i - 1) for j in range(r+1)])
            
            return {i : fund_weights(i+1) for i, x in self.reference_sheets}

        elif algebra[0] == 'D':
            def pos_fund_weights(i):
                return np.array([kr_delta(j, i - 1) for j in range(r)])
            
            def neg_fund_weights(i):
                return -1 * pos_fund_weights(i)

            positive_sheets = [[i, x] for i, x in self.reference_sheets if d_positivity(x)]
            negative_sheets = [[i, x] for i, x in self.reference_sheets if not d_positivity(x)]
            sorted_negative_sheets = sort_negatives(positive_sheets, negative_sheets)
            print "\nThe positive sheets at z_0:"
            print positive_sheets
            print "The corresponding sorted negative sheets:"
            print sorted_negative_sheets
            pos_dict = {i : pos_fund_weights(j+1) for j, [i, x] in enumerate(positive_sheets)}
            # print "\npositive dictionary"
            # print pos_dict
            neg_dict = {i : neg_fund_weights(j+1) for j, [i, x] in enumerate(sorted_negative_sheets)}
            # print "\nnegative dictionary"
            # print neg_dict
            full_dict = pos_dict.copy()
            full_dict.update(neg_dict)
            return full_dict


        elif algebra[0] == 'E':
            raise ValueError('I am not ready for E-type algebras yet!')


    def positive_root(self, pairs, singles, enum_sh):
        """
        Determines the positive root associated with 
        a branch point's 'structure', i.e. how the sheets
        collide at the branch point
        """
        algebra = self.algebra
        first_pair = pairs[0]
        i_1 = first_pair[0]
        i_2 = first_pair[1]
        v_1 = self.sheet_weight_dictionary[i_1]
        v_2 = self.sheet_weight_dictionary[i_2]
        
        if algebra[0] == 'A':
            if i_1 < i_2:
                return v_1 - v_2
            else:
                return v_2 - v_1

        elif algebra[0] == 'D':
            return d_positive_root(v_1 - v_2)

        elif algebra[0] == 'E':
            raise ValueError('I am not ready for E-type algebras yet!')


def bp_from_ramif(ramification_points):
    b_points = [r.z for r in ramification_points if not r.is_puncture]
    seen = set()
    uniq = []
    for b in b_points:
        if b not in seen:
            uniq.append(b)
            seen.add(b)
    return uniq

def irr_sing_from_ramif(ramification_points):
    irr_sing = [r.z for r in ramification_points if r.is_puncture]
    seen = set()
    uniq = []
    for irr in irr_sing:
        if irr not in seen:
            uniq.append(irr)
            seen.add(irr)
    return uniq



def kr_delta(i, j):
    if i == j:
        return 1
    else:
        return 0

def getkey_real(item):
    return item[1].real

def getkey_imag(item):
    return item[1].imag

def getkey_last(item):
    return item[-1]

def flatten(l):
    return [item for sublist in l for item in sublist]

def d_positivity(x):
    """
    Criterion for establishing whether a sheet of a D-type 
    vector-representation cover is positive or not.
    """
    if x.imag > 0 or (x.imag == 0 and x.real > 0):
        return True
    elif x.imag < 0 or (x.imag == 0 and x.real < 0):
        return False
    else:
        raise ValueError('There is a sheet located at x=0.0 for a D-type cover!'
                        +'\nBy Z_2 symmetry there must be two sheets at x=0.'
                        +'\nThis means that the basepoint is not good.')

def d_positive_root(alpha):
    """
    Given a root 'alpha', it determines whether it's positive or negative.
    If positive, it returns alpha. If negative it returns -alpha
    """
    non_zero_elements = [x for x in alpha if x!=0]
    first_element = non_zero_elements[0]
    if first_element > 0:
        return alpha
    else:
        return -1 * alpha



def sort_negatives(pos, neg):
    sorted_neg = []
    for i, x in pos:
        distances = [[j, y, abs(x - (-y))] for j, y in neg]
        ### sorting them by distance, the closest will be the first one
        sorted_distances = sorted(distances, key=getkey_last)
        closest_sheet = sorted_distances[0][:2]
        sorted_neg.append(closest_sheet)

    return sorted_neg


def data_plot(cmplx_list, title):
    """
    Plot the real and imaginary parts of 
    a list of complex numers
    """

    l = len(cmplx_list)
    r_list = [x.real for x in cmplx_list]
    i_list = [x.imag for x in cmplx_list]

    plt.figure(1)

    plt.subplot(211)
    plt.plot(r_list, "r.")
    r_delta = abs(max(r_list)-min(r_list))
    plt.axis([0.0, float(l), min(r_list) - 0.15 * r_delta, \
                             max(r_list) + 0.15 * r_delta])
    plt.ylabel("real part")
    plt.title(title)

    plt.subplot(212)
    plt.plot(i_list, "r.")
    i_delta = abs(max(i_list)-min(i_list))
    plt.axis([0.0, float(l), min(i_list) - 0.15 * i_delta, \
                             max(i_list) + 0.15 * i_delta])
    plt.ylabel("imaginary part")
    plt.title(title)
    
    plt.show()


### Some Testing

# config = load_config('../default.ini')
# algebra = ['A', 2]

# config = load_config('../config/pure_SO_4.ini')
# algebra = ['D', 2]

config = load_config('../config/coset_D_3.ini')
algebra = ['D', 3]

# config = load_config('../config/coset_D_4.ini')
# algebra = ['D', 4]

sw = SWData(config)
ramification_points = get_ramification_points(sw, config['accuracy'])


print "\nThe SW curve in the FIRST fundamental representation"
print sw.curve.sym_eq
print sw.curve.num_eq

print "\nThe ramification points"
for r in ramification_points:
    print r


SHOW_TRACKING_PLOTS = False
t = Trivialization(sw, ramification_points, algebra)

z_arb = 1.50 + 2.35j
print "\nThe sheets trivializing the FUNDAMENTAL cover at z = %s" % z_arb
print t.sheets_at_arbitrary_z(z_arb)

