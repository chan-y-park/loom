# Force integer division to give a float, i.e. 1/2 = 0.5.
from __future__ import division
from sympy import sympify, poly
from sympy.mpmath import mp, mpc
# from sympy.parsing.sympy_parser import parse_expr
import logging
import os
import subprocess
import pdb

base_dir = os.path.dirname(os.path.realpath(__file__))
sage_script_dir = base_dir + '/sage_scripts/'

def solve_system_of_eqs(eqs, precision=None, logger_name='loom',):
    """
    Use sage to solve the given system of polynomial equations of x and z.
    """
    logger = logging.getLogger(logger_name)
    sols = []
    if precision is not None:
        mp.dps = precision
    else:
        precision = 15
    logger.info('Use SAGE to solve {} @ precision = {}.'
                 .format(eqs, precision))
    try:
        rv_str = subprocess.check_output(
            ['sage', sage_script_dir + 'solve_system_of_eqs.sage'] +
            [str(precision)] +
            [str(eq) for eq in eqs]
        )
    except (KeyboardInterrupt, SystemExit):
        raise
    
    rv = eval(rv_str)
    sols_str_list, messages = rv

    for msg in messages:
        logger.warning(msg)

    #sols_str_list = eval(sols_str_list_str)
    for sols_str in sols_str_list:
        (z_re, z_im), (x_re, x_im) = sols_str
        sols.append(
            (mpc(z_re, z_im), mpc(x_re, x_im))
        )

    return sols


# TODO: Unify this and the next function 
# with the above method for multi-equations
# these were written in a rush to finish something else, 
# apologies for the repetition of code :)
# But Careful for the subtle differences, not just relabelings..
# In particular the use of the polynomial ring and of the
# method called .roots()
# Also, cleanup the messages log forwarding in here, it's probably useless.
def solve_single_eq_z(eqs, precision=None, logger_name='loom',):
    """
    Use sage to solve a single polynomial equation in z.
    """
    logger = logging.getLogger(logger_name)
    sols = []
    if precision is not None:
        mp.dps = precision
    else:
        precision = 15
    try:
        rv_str = subprocess.check_output(
            ['sage', sage_script_dir + 'solve_single_eq_z.sage'] +
            [str(precision)] +
            [str(eq) for eq in eqs]
        )
    except (KeyboardInterrupt, SystemExit):
        raise
    
    rv = eval(rv_str)
    sols_str_list, mult_str_list, messages = rv

    for msg in messages:
        logger.warning(msg)
    
    for i, sols_str in enumerate(sols_str_list):
        (z_re, z_im) = sols_str
        for j in range(mult_str_list[i]):
            sols.append(
                mpc(z_re, z_im)
            )

    return sols


def solve_single_eq_x(eqs, precision=None, logger_name='loom',):
    """
    Use sage to solve a single polynomial equation in x.
    """
    logger = logging.getLogger(logger_name)
    sols = []
    if precision is not None:
        mp.dps = precision
    else:
        precision = 15
    try:
        rv_str = subprocess.check_output(
            ['sage', sage_script_dir + 'solve_single_eq_x.sage'] +
            [str(precision)] +
            [str(eq) for eq in eqs]
        )
    except (KeyboardInterrupt, SystemExit):
        raise
    
    rv = eval(rv_str)
    sols_str_list, mult_str_list, messages = rv

    for msg in messages:
        logger.warning(msg)
    
    for i, sols_str in enumerate(sols_str_list):
        (z_re, z_im) = sols_str
        for j in range(mult_str_list[i]):
            sols.append(
                mpc(z_re, z_im)
            )

    return sols


def get_g_data(root_system, highest_weight):
    try:
        g_data_str = subprocess.check_output(
            ['sage', sage_script_dir + 'get_g_data.sage', root_system,
             str(highest_weight)]
        )
    except (KeyboardInterrupt, SystemExit):
        raise

    g_data = eval(g_data_str)
    
    return g_data


def compute_discriminant(f):
    """
    Use SAGE to compute the discriminant of a polynomial f.
    f must be expressed in variables x, z.
    The discriminant will be computed with respect to x.
    """

    # TODO: find a more civilized way, instead of this very dirty trick
    f_list = list(str(f))
    for i, letter in enumerate(f_list):
        if letter == 'I':
            f_list[i] = '1j'
    f_str = "".join(f_list)

    try:
        disc_str = subprocess.check_output(
            ['sage', sage_script_dir + 'compute_discriminant.sage'] +
            [f_str]
        )
    except (KeyboardInterrupt, SystemExit):
        raise
    print 'tha answer of sage'
    print disc_str
    disc_sym = sympify(disc_str)
    if disc_sym == 0:
        return 0
    else:
        return poly(sympify(disc_str))


