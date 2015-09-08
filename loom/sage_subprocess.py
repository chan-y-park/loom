# Force integer division to give a float, i.e. 1/2 = 0.5.
from __future__ import division

import os
import subprocess
import pdb

base_dir = os.path.dirname(os.path.realpath(__file__))
sage_script_dir = base_dir + '/sage_scripts/'

def solve_poly_system(poly_system):
    """
    Use sage to solve the given system of polynomial equations of x and z.
    """
    sols_str = subprocess.check_output(
        ['sage', sage_script_dir + 'solve_poly_system.sage'] +
        [str(poly) for poly in poly_system]
    )
    sols = eval(sols_str)
    return sols


def get_g_data(root_system, highest_weight):
    g_data_str = subprocess.check_output(
        ['sage', sage_script_dir + 'get_g_data.sage', root_system,
         str(highest_weight)]
    )
    g_data = eval(g_data_str)
    
    return g_data

