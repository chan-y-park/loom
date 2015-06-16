from api import load_config
from geometry import SWData, get_ramification_points
import sympy
import matplotlib.pyplot as plt
import cmath
import numpy as np
from sympy import Poly

### number of steps used to track the sheets along a leg 
### of the lefshetz spider
N_PATH_TO_BP = 50

### Tolerance for recognizing colliding sheets at a branch-point
PROXIMITY_THRESHOLD = 0.05


config = load_config('../default.ini')
sw = SWData(config)
ramification_points = get_ramification_points(sw, config['accuracy'])


class Trivialization:
    ### NOTE: I am assuming that branch points do not overlap vertically
    ### this should be guaranteed by introducing an automatic rotation of 
    ### the z-plane before calling this class.
    ### NOTE: I am assuming the ramification_points are all branch points 
    ### for now, not irregular singularities.
    def __init__(self, sw_data, ramification_points):
        self.sw_data = sw_data
        self.ramification_points = ramification_points
        r_points_z = [r.z for r in self.ramification_points]
        max_distance = max([abs(x - y) for x in r_points_z 
                                                    for y in r_points_z])
        
        ### Should generalize to take into account punctures as well?
        r_center = sum(r_points_z)
        self.basepoint = r_center - 1j * max_distance

        self.reference_sheets = [[i, x] \
                    for i, x in enumerate(self.sheets_at_z(self.basepoint))]

        print "\nSheets of the cover at z_0 = {}".format(self.basepoint)
        print self.reference_sheets

        print "\npath to first branch point"
        print self.path_to_bp(r_points_z[0])

        for z_bp in r_points_z:
            # sheets_path = self.track_sheets_to_bp(z_bp)
            # for sheet_list in sheets_path:
            #     data_plot(sheet_list, z_bp)
            pass

        pairs, singles  = self.branch_point_structure(r_points_z[0])
        print "\npairs = {}".format(pairs)
        print "singles = {}".format(singles)



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

    # def branch_point_structure(self, z_bp):
    #     tracked_sheets = self.track_sheets_to_bp(z_bp)
    #     sheets_at_bp = [sheet_list[-1] for sheet_list in tracked_sheets]
    #     enum_sh = [[i, s_i] for i, s_i in enumerate(sheets_at_bp)]

    #     print "sheets of the cover at z = {}".format(z_bp)
    #     print enum_sh
    #     ### TO DO:
    #     ### Will sort sheets according to their real part in the fiber C-plane
    #     ### then will only compare adjacent ones
    #     ### But we should check that not more than two roots have the same 
    #     ### real part. If that's the case, should tilt the x-plane
    #     ### Alternative: don't use this sorting at all, and compare all posible 
    #     ### pairs of sheets.
    #     real_sorted_enum_sh = sorted(enum_sh, key=getkey_real)
        
    #     ### These 'pairs' and 'singles' should be probably
    #     ### introduced as attributes of the branch-point class
    #     pairs = []
    #     singles = []

    #     for j, sheet in enumerate(real_sorted_enum_sh):
    #         ### here 'sheet' stands for the sheet identifier
    #         ### and contains sheet = [i, x]
    #         ### while 'j' is the index of the sheets sorted
    #         ### according to the real parts of x
    #         i = sheet[0]
    #         x = sheet[1]
                        
    #         if i in flatten(pairs):
    #             pass
    #         elif j == len(real_sorted_enum_sh)-1:
    #             ### this is the last sheet of the list
    #             ### if it's not already in a pair, then it's a single
    #             singles.append(i)
    #         else:
    #             x_1 = real_sorted_enum_sh[j+1][1]
    #             i_1 = real_sorted_enum_sh[j+1][0]
                
    #             if abs(x - x_1) < PROXIMITY_THRESHOLD:
    #                 pairs.append([i, i_1])
    #             else:
    #                 singles.append(i)

    #     return pairs, singles

    def branch_point_structure(self, z_bp):
        tracked_sheets = self.track_sheets_to_bp(z_bp)
        sheets_at_bp = [sheet_list[-1] for sheet_list in tracked_sheets]
        enum_sh = [[i, s_i] for i, s_i in enumerate(sheets_at_bp)]

        print "sheets of the cover at z = {}".format(z_bp)
        print enum_sh
        
        ### These 'pairs' and 'singles' should be probably
        ### introduced as attributes of the branch-point class
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

        return pairs, singles




def getkey_real(item):
    return item[1].real

def getkey_imag(item):
    return item[1].imag

def flatten(l):
    return [item for sublist in l for item in sublist]


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


print "\nThe SW curve"
print sw.curve.sym_eq
print sw.curve.num_eq

print "\nThe ramification points"
for r in ramification_points:
    print r

t = Trivialization(sw, ramification_points)
