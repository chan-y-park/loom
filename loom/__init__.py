import os, platform
import matplotlib
if matplotlib.rcParams['backend'] == 'nbAgg':
    print('Use IPython notebook backend for matplotlib.')
elif platform.system() == 'Linux':
    try:
        os.environ['DISPLAY']
        matplotlib.use('TkAgg')
        print('Use TkAgg backend for matplotlib.')
    except KeyError:
        matplotlib.use('Agg') 
        print('Use Agg backend for matplotlib.')
else:
    print('Use default backend defined in matplotlibrc: '
          '{}'.format(matplotlib.rcParams['backend']))

from __main__ import run
from api import set_logging
from api import generate_spectral_network as generate
from api import load_spectral_network as load
from api import save_spectral_network as save
from api import make_spectral_network_plot as plot
from api import load_config, save_config

print('Setting logging level to "INFO".')
set_logging('info')

__all__ = [
    'config',
    'data_io',
    'gui',
    'intersection',
    'misc',
    'plotting',
    'spectral_network',
    's_wall',
    'run',
    'set_logging',
    'load_config',
    'save_config',
    'load',
    'save',
    'generate',
    'plot',
]
