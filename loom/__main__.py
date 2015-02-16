"""
Main module. Run 'python -m loom' at the directory that contains the
directory 'loom'.

This module is based on the following libraries:
NumPy 1.8.2
SciPy 0.15.1
SymPy 0.7.6 -> must be installed
"""
import sys
import logging
import os
import getopt
import time

from gui import open_gui
from api import load_spectral_network, generate_spectral_network

shortopts = 'c:gl:p:'
longopts = [
    'gui-mode',
    'logging-level',
    'phase',
    'show-plot',
    'show-plot-on-cylinder',
    'load-data=',
    #'save-data=',
]

CONFIG_FILE_DIR = 'config_file'
DATA_FILE_DIR = 'data_file'


def run_with_optlist(optlist):

    opts = {
        'config-file': None,
        'gui-mode': False,
        'logging-level': 'warning',
        'phase': None,
        'load-data': '',
        'show-plot': False,
        'show-plot-on-cylinder': False,
    }
        
    if len(optlist) == 0:
        print("""usage: python -m loom [OPTION]

    -c CFG_FILE_NAME:
        Read CFG_FILE_NAME to set up the configuration.
        When this option is not set, by default use 
            loom/config_file/default.ini

    -g, --gui-mode:
        Run the graphical user interface.

    -l LEVEL, --logging-level=LEVEL:
        Set logging level to LEVEL. 'warning' is default.

    -p, --phase THETA:
        Generate a spectral network at the phase of THETA.
        Overrides 'phase_range' of the configuration file. 

    --load-data DATA_FILE:
        Load data from a file.
        
    --show-plot:
        Diaplay the spectral network plot.

    --show-plot-on-cylinder:
        Diaplay the spectral network plot on a cylinder,
        assuming that the given Seiberg-Witten differential
        has two punctures at 0 and infinity.
        """)

    else:

        for opt, arg in optlist:
            if (opt == '-c' and len(arg) > 0):
                opts['config-file'] = arg
            elif (opt == '-g' or opt == '--gui-mode'):
                opts['gui-mode'] = True
            elif (opt == '-l' or opt == '--logging-level'):
                opts['logging-level'] = arg
            elif (opt == '-p' or opt == '--phase'):
                opts['phase'] = float(arg)
            elif opt == '--load-data':
                opts['load-data'] = arg
            elif opt == '--show-plot':
                opts['show-plot'] = True
            elif opt == '--show-plot-on-cylinder':
                opts['show-plot-on-cylinder'] = True
        # End of option setting.

        # Set logging.
        if opts['logging-level'] == 'debug':
            logging_level = logging.DEBUG
            logging_format = '%(module)s@%(lineno)d: %(funcName)s: %(message)s'
        elif opts['logging-level'] == 'info':
            logging_level = logging.INFO
            logging_format = '%(process)d: %(message)s'
        else:
            logging_level = logging.WARNING
            logging_format = '%(message)s'

        logging.basicConfig(level=logging_level, format=logging_format, 
                            stream=sys.stdout)

        # Entry point branching
        if opts['gui-mode'] is True:
            return open_gui(opts)
        elif (len(opts['load-data']) == 0):
            logging.warning('option "load-data" requires an argument.')
        elif (len(opts['load-data']) > 0):
            return load_spectral_network(opts)
        else:
            return generate_spectral_network(opts)

# Set options from sys.argv when running on the command line,
# then start running the main code.
def run_with_sys_argv(argv):    
    try:
        optlist, args = getopt.getopt(argv, shortopts, longopts,)
        return run_with_optlist(optlist)

    except getopt.GetoptError:
        print 'Unknown options.'

# Set options from string 'optstr' when running on the interpreter, 
# then start running the main code.
def run(optstr=''):
    return run_with_sys_argv(optstr.split())

# End of definitions

if __name__ == '__main__':
    run_with_sys_argv(sys.argv[1:])
