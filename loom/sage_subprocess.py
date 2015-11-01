# Force integer division to give a float, i.e. 1/2 = 0.5.
from __future__ import division
import sympy.mpmath as mpmath
from sympy.mpmath import mp, mpc
import os
import subprocess
import pdb

base_dir = os.path.dirname(os.path.realpath(__file__))
sage_script_dir = base_dir + '/sage_scripts/'

def solve_system_of_eqs(eqs, precision=None):
    """
    Use sage to solve the given system of polynomial equations of x and z.
    """
    sols = []
    if precision is not None:
        mp.dps = precision
    else:
        precision = 15
    try:
        sols_str_list_str = subprocess.check_output(
            ['sage', sage_script_dir + 'solve_system_of_eqs.sage'] +
            [str(precision)] +
            [str(eq) for eq in eqs]
        )
    except (KeyboardInterrupt, SystemExit) as e:
        raise

    sols_str_list = eval(sols_str_list_str)
    for sols_str in sols_str_list:
        (z_re, z_im), (x_re, x_im) = sols_str
        sols.append(
            (mpc(z_re, z_im), mpc(x_re, x_im))
        )

    return sols


def get_g_data(root_system, highest_weight):
    try:
        g_data_str = subprocess.check_output(
            ['sage', sage_script_dir + 'get_g_data.sage', root_system,
             str(highest_weight)]
        )
    except (KeyboardInterrupt, SystemExit) as e:
        raise

    g_data = eval(g_data_str)
    
    return g_data

