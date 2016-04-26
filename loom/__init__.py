import os
import time
import platform
#import logging
import sys
import matplotlib
if matplotlib.rcParams['backend'] == 'nbAgg':
    print('Use IPython notebook backend for matplotlib.')
elif matplotlib.rcParams['backend'] == 'WebAgg':
    print('Use WebAgg backend for matplotlib.')
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

# matplotlib configuration
matplotlib.rcParams['figure.figsize'] = [8.0, 8.0]
matplotlib.rcParams['savefig.format'] = 'pdf'
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.pad_inches'] = 0.0

from api import get_loom_dir
from api import set_logging
from api import generate_spectral_network as generate
from api import load_spectral_network as load
from api import save_spectral_network as save
from api import make_spectral_network_plot as plot
from api import load_config, save_config

LOGGING_FILE_PATH = os.path.join(
    get_loom_dir(),
    ('logs/loom_{}-{:02}-{:02} {:02}:{:02}:{:02}.log'
     .format(*time.localtime(time.time())[:6])),
)

set_logging(
    logger_name='loom',
    logging_level_name='INFO',
    logging_stream=sys.stdout,
    logging_file_name=LOGGING_FILE_PATH,
)

__all__ = [
    'config',
    'gui',
    'intersection',
    'misc',
    'plotting',
    'spectral_network',
    's_wall',
    'set_logging',
    'load_config',
    'save_config',
    'load',
    'save',
    'generate',
    'plot',
]
