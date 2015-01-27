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

from config import ConfigData 
from gui import open_gui
from spectral_network import generate_spectral_network

shortopts = 'c:gl:p:w'
longopts = [
    'gui-mode',
    'logging-level',
    'phase',
    'show-plot',
]

def run_with_optlist(optlist):

    opts = {
        'config-file': None,
        'gui-mode': False,
        'logging-level': 'warning',
        'phase': None,
        #'single-network': False,
        'save-data': False,
        'show-plot': False,
    }
        
    if len(optlist) == 0:
        print("""usage: python -m loom [OPTION]

    -c CFG_FILE_NAME:
        read CFG_FILE_NAME to set up the configuration.

    -g, --gui-mode:
        run the graphical user interface.

    -l LEVEL, --logging-level=LEVEL:
        set logging level to LEVEL. 'warning' is default.

    -p, --phase THETA:
        generate a spectral network at the phase of THETA.

    -w:
        save data on a time-stamped file.
        
    --show-plot:
        diaplay the spectral network plot.
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
                #opts['single-network'] = True
            elif opt == '-w':
                opts['save-data'] = True
            elif opt == '--show-plot':
                opts['show-plot'] = True
        # End of option setting.

        config_data = ConfigData(opts)

        if opts['gui-mode'] == True:
            open_gui(opts, config_data)
        else:
            generate_spectral_network(opts, config_data)

# Set options from sys.argv when running on the command line,
# then start running the main code.
def run_with_sys_argv(argv):    
    try:
        optlist, args = getopt.getopt(argv, shortopts, longopts,)
        run_with_optlist(optlist)

    except getopt.GetoptError:
        print 'Unknown options.'

# Set options from string 'optstr' when running on the interpreter, 
# then start running the main code.
def run_with_optstr(optstr):
    run_with_sys_argv(optstr.split())

# End of definitions

run_with_sys_argv(sys.argv[1:])
