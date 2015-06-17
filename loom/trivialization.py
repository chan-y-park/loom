from api import load_config
from geometry import SWData, get_ramification_points, RamificationPoint
import sympy
import matplotlib.pyplot as plt
import cmath
import numpy as np
from sympy import Poly

### number of steps used to track the sheets along a leg 
### of the lefshetz spider
N_PATH_TO_BP = 100

### Tolerance for recognizing colliding sheets at a branch-point
PROXIMITY_THRESHOLD = 0.05



class BranchPoint:
    def __init__(self, z=None, structure=None, trivialization=None):
        self.z = z
        self.pairs = structure[0]
        self.singles = structure[1]
        self.sheets_at_branch_point = structure[2]
        self.trivialization = trivialization
        self.positive_root = trivialization.positive_root(structure)
        self.sheet_tracks_from_basepoint = trivialization.track_sheets_to_bp(z)
        
    


class Trivialization:
    """
    The trivialization class.

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

    Attributes
    ----------

    basepoint : 
    the base point of the trivialization

    reference_sheets :
    a list of pairs [i, x] where 'i' is an integer label
    for the sheet, and 'x' is its position in the fiber of T^*C 
    over the basepoint

    """
    ### NOTE: I am assuming that branch points do not overlap vertically
    ### this should be guaranteed by introducing an automatic rotation of 
    ### the z-plane before calling this class.
    ### NOTE: I am only analyzing the ramification_points that are branch points 
    ### for now, not irregular singularities.
    ### NOTE: I am restricting to square-root type branch points.
    ### Although I am not printing any explicit warning/error message 
    ### and the computation will go through for higher-type, but give a wrong answer!
    def __init__(self, sw_data, ramification_points, lie_algebra):
        self.sw_data = sw_data
        self.algebra = lie_algebra
        self.branch_points = []
        self.sheet_weight_dictionary = None

        b_points = [r for r in ramification_points if not r.is_puncture]
        b_points_z = [b.z for b in b_points]
        max_distance = max([abs(x - y) for x in b_points_z 
                                                    for y in b_points_z])
        
        ### Should generalize to take into account punctures as well?
        b_center = sum(b_points_z)
        self.basepoint = b_center - 1j * max_distance

        self.reference_sheets = [[i, x] \
                    for i, x in enumerate(self.sheets_at_z(self.basepoint))]
        ### Now that we have determined the reference sheets, 
        ### we build a weight-sheet dictionary
        self.sheet_weight_dictionary = self.build_dictionary()

        for i, z_bp in enumerate(b_points_z):
            bp_structure = self.branch_point_structure(z_bp)
            self.branch_points.append(BranchPoint(
                                                    z=z_bp, 
                                                    structure=bp_structure,
                                                    trivialization=self
                                                ))

        ### PRESENT THE DATA FOR DEBUGGING PURPOSES
        print "\nSheets of the cover at z_0 = {}".format(self.basepoint)
        print self.reference_sheets

        print "\nThe dictionary between sheets and weights:"
        print self.sheet_weight_dictionary

        for i, bp in enumerate(self.branch_points):
            print "\nroot tracking along the path to branch point #{}".format(i)
            for j, sheet_list in enumerate(bp.sheet_tracks_from_basepoint):
                data_plot(sheet_list, 'sheet {} tracket to branch point {}'.format(j, i))
            print "this is the branch point structure for branch point #{}".format(i)
            print "pairs = {}".format(bp.pairs)
            print "singles = {}".format(bp.singles)
            print "positive root = {}".format(bp.positive_root)
            print "sheets at the branch point = {}".format(bp.sheets_at_branch_point)

        
    def sheets_at_z(self, z_0):
        from sympy.abc import x, z
        sw_curve_fiber = self.sw_data.curve.num_eq.subs(z, z_0)
        sym_poly = Poly(sw_curve_fiber, x, domain='CC')
        coeff_list = map(complex, sym_poly.all_coeffs())
        return map(complex, np.roots(coeff_list))


    def path_to_bp(self, z_bp):
        z_0 = self.basepoint
        z_1 = 1j * self.basepoint.imag + z_bp.real
        z_2 = z_bp
        half_steps = int(N_PATH_TO_BP / 2)
        return [z_0 + ((z_1 - z_0) / half_steps) * i \
                                        for i in range(half_steps + 1)] \
            + [z_1 + ((z_2 - z_1) / half_steps) * i \
                                        for i in range(half_steps + 1)]

    def track_sheets_to_bp(self, z_bp):
        sheets_0 = self.sheets_at_z(self.basepoint)
        z_path = self.path_to_bp(z_bp)
        sheets_along_path = [[s] for s in sheets_0]

        for z in z_path:
            sheets_1 = self.sheets_at_z(z)
            sheets_0 = self.sort_sheets(sheets_1, sheets_0)
            for i, s_list in enumerate(sheets_along_path):
                s_list.append(sheets_0[i])

        return sheets_along_path
        ### the result is of the form [sheet_path_1, sheet_path_2, ...]
        ### where sheet_path_i = [x_0, x_1, ...] are the fiber coordinates
        ### of the sheet along the path

    def sort_sheets(self, new_sheets, ref_sheets):
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
        ### SHOULD INTRODUCE A CHECK THAT A NEW_SHEET CAN'T BE PICKED TWICE!
        return sorted_sheets


    def branch_point_structure(self, z_bp):
        tracked_sheets = self.track_sheets_to_bp(z_bp)
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
                    if abs(x - y) < PROXIMITY_THRESHOLD:
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

        return pairs, singles, enum_sh

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
            print "The positive sheets:"
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


    def positive_root(self, structure):
        """
        Determines the positive root associated with 
        a branch point's 'structure', i.e. how the sheets
        collide at the branch point
        """
        pairs, singles, enum_sh = structure
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


sw = SWData(config)
ramification_points = get_ramification_points(sw, config['accuracy'])


print "\nThe SW curve in the FIRST fundamental representation"
print sw.curve.sym_eq
print sw.curve.num_eq

print "\nThe ramification points"
for r in ramification_points:
    print r



t = Trivialization(sw, ramification_points, algebra)
