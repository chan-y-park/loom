"""
Entry point of 'loom'.

Requires the following libraries:
NumPy 1.8.2
SciPy 0.15.1
SymPy 0.7.6
"""
import sys
import logging
import os
import getopt
import time
import pdb

from config import LoomConfig
from gui import open_gui
from api import (
    set_logging, load_config, load_spectral_network, generate_spectral_network,
    save_spectral_network,
)

shortopts = 'c:g:hl:p:'
longopts = [
    'gui-mode=',
    'help',
    'logging-level',
    'phase=',
    'show-plot',
    'show-plot-on-cylinder',
    'load-data=',
]


def run_with_optlist(optlist):
    if len(optlist) == 0:
        return print_help()

    opts = {
        'config-file': None,
        'gui-mode': False,
        'logging-level': 'info',
        'phase': None,
        'load-data': None,
        'show-plot': False,
        'show-plot-on-cylinder': False,
    }

    for opt, arg in optlist:
        if (opt == '-c' and len(arg) > 0):
            opts['config-file'] = arg
        elif (opt == '-g' or opt == '--gui-mode'):
            opts['gui-mode'] = eval(arg)
        elif (opt == '-h' or opt == '--help'):
            return print_help()
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

    set_logging(opts['logging-level'])

    # Load config & data according to the command line args.
    if opts['load-data'] is not None:
        data_dir = opts['load-data']
        config, spectral_networks = load_spectral_network(data_dir)
    else:
        # Read in the specified config file.
        config_file = opts['config-file']
        # No config file chosen; read the default config file.
        if config_file is None:
            config_file = os.path.join(os.curdir, 'default.ini')
        config = load_config(config_file)
        spectral_networks = []

    # Open GUI.
    if opts['gui-mode'] is True:
        return open_gui(config, spectral_networks)

    # Start in a command line mode and generate spectral networks
    # according to the given arguments.
    if len(spectral_networks) == 0:
        spectral_networks = generate_spectral_network(
            config,
            phase=opts['phase'],
        )
        if (opts['show-plot'] or opts['show-plot-on-cylinder']) is True:
            make_spectral_network_plot(
                config,
                spectral_network_list,
                plot_on_cylinder=opts['show-plot-on-cylinder'],
            )
        if opts['gui-mode'] is False:
            save_spectral_network(config, spectral_networks,
                                  data_dir=None, make_zipped_file=True,)

    return (config, spectral_networks)

def run_with_sys_argv(argv):
    """
    Set options from sys.argv when running on the command line,
    then start running the main code.
    """
    try:
        optlist, args = getopt.getopt(argv, shortopts, longopts,)
        return run_with_optlist(optlist)

    except getopt.GetoptError:
        print 'Unknown options.'

def run(optstr=''):
    """
    Set options from string 'optstr' when running on the interpreter,
    then start running the main code.
    """
    return run_with_sys_argv(optstr.split())

def print_help():
    print("""usage: loom [OPTION]

-c CFG_FILE_NAME:
    Read CFG_FILE_NAME to set up the configuration.
    When this option is not set, by default use
        loom/config_file/default.ini

-g VALUE, --gui-mode=VALUE:
    Run the graphical user interface if VALUE=True,
    otherwise start in a command line mode.

-l LEVEL, --logging-level=LEVEL:
    Set logging level to LEVEL. 'warning' is default.

-p, --phase THETA:
    Generate a spectral network at the phase of THETA.
    Overrides 'phase_range' of the configuration file.

--load-data DATA_FILE:
    Load data from a file. If DATA_FILE is not specified,
    shows a file dialog window to open a directory
    containing data files.

--show-plot:
    Diaplay the spectral network plot.

--show-plot-on-cylinder:
    Diaplay the spectral network plot on a cylinder,
    assuming that the given Seiberg-Witten differential
    has two punctures at 0 and infinity.
    """)


# End of definitions

if __name__ == '__main__':
    run_with_sys_argv(sys.argv[1:])
