import time
import sys
import os
import glob
import zipfile, zlib
import logging
import Tkinter as tk
import tkFileDialog
import matplotlib
import pdb

from config import LoomConfig
#from geometry import SWData, get_ramification_points
from trivialization import SWDataWithTrivialization
from spectral_network import SpectralNetwork
from parallel import parallel_get_spectral_network
from plotting import NetworkPlot, NetworkPlotTk

LOGGING_FILE_NAME = 'logs/log.mose.txt'

def set_logging(level):
    if level == 'debug':
        logging_level = logging.DEBUG
        logging_format = '%(module)s@%(lineno)d: %(funcName)s: %(message)s'
    elif level == 'info':
        logging_level = logging.INFO
        logging_format = '%(process)d: %(message)s'
    else:
        logging_level = logging.WARNING
        logging_format = '%(message)s'

    #logging.basicConfig(level=logging_level, format=logging_format, 
    #                    stream=sys.stdout)
    logger = logging.getLogger()
    # Remove other handlers
    for handler in logger.handlers:
        logger.removeHandler(handler)
    logger.setLevel(logging_level)
    formatter = logging.Formatter(logging_format)
    ### create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging_level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    # Log to a file 'log.mose.txt'
    fh = logging.FileHandler(LOGGING_FILE_NAME, 'w')
    fh.setLevel(logging_level)
    fh.setFormatter(formatter)
    logger.addHandler(fh)


def generate_spectral_network(config, phase=None):
    """
    Generate one or more spectral networks according to
    the command-line options and the configuration file
    and return a list of data obtained from SpectralNetwork.get_data()
    """

    phase_range = config['phase_range']
    #sw = SWData(config)
    #ramification_points = get_ramification_points(sw, config['accuracy'])
    sw = SWDataWithTrivialization(config)

    start_time = time.time()
    logging.info('start cpu time: %s', start_time)

    if(phase is not None):
        logging.info('Generate a single spectral network at theta = {}.'
                     .format(phase))
        spectral_network = SpectralNetwork(
            phase=phase, 
        ) 

        spectral_network.grow(sw, config)

        spectral_network_list = [spectral_network]

    elif(phase_range is not None):
        logging.info('Generate multiple spectral networks.')
        logging.info('phase_range = {}.'.format(phase_range))
        spectral_network_list = parallel_get_spectral_network(
            sw, 
            config,
        ) 

    end_time = time.time()
    logging.info('end cpu time: %.8f', end_time)
    logging.info('elapsed cpu time: %.8f', end_time - start_time)

    return spectral_network_list


def load_config(config_file=None):
    if config_file is None:
        root = tk.Tk()
        file_opts = {
            'defaultextension': '.ini',
            'initialdir': os.curdir,
            'initialfile': 'config.ini',
            'parent': root,
            'title': 'Select a configuration file to load.',
        }
        config_file = tkFileDialog.askopenfilename(**file_opts)
        root.destroy()
        if config_file == '':
            return None

    config = LoomConfig()
    config.read(config_file)

    return config
    

def load_spectral_network(data_dir=None):
    if data_dir is None:
        root = tk.Tk()
        dir_opts = {
            'initialdir': os.curdir,
            'mustexist': False,
            'parent': root,
            'title': 'Select a directory that contains data files.',
        }
        data_dir = tkFileDialog.askdirectory(**dir_opts)
        root.destroy()
        if data_dir == '':
            return (None, None)
 
    logging.info('Opening data directory "{}"...'.format(data_dir))

    config = LoomConfig()
    config.read(os.path.join(data_dir, 'config.ini'))

    sw = SWDataWithTrivialization(config)
    spectral_network_list = []

    data_file_list = glob.glob(os.path.join(data_dir, 'data_*.json'))
    data_file_list.sort()
    for data_file in data_file_list:
        logging.info('Loading {}...'.format(data_file))
        spectral_network = SpectralNetwork(
            config=config,
        )
        with open(data_file, 'r') as fp:
            spectral_network.set_from_json_data(fp)
            spectral_network_list.append(spectral_network)

    return (config, spectral_network_list)


def save_config(config, path=None):
    if path is None:
        root = tk.Tk()
        file_opts = {
            'defaultextension': '.ini',
            'initialdir': os.curdir,
            'initialfile': 'config.ini',
            'parent': root,
            'title': 'Save the current configuration to a file.',
        }
        path = tkFileDialog.asksaveasfilename(**file_opts)
        root.destroy()
        if path == '':
            return None

    logging.info('Save configuration to {}.'.format(path))
    with open(path, 'wb') as fp:
        config.parser.write(fp)

    return None


def save_spectral_network(config, spectral_networks, data_dir=None,
                          make_zipped_file=True):
    if data_dir is None:
        # Prepare to save spectral network data to files.
        timestamp = str(int(time.time()))
        data_dir = os.path.join(
            os.curdir,
            'data', 
            timestamp
        )

    logging.info('Make a directory {} to save data.'.format(data_dir))
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # Save configuration to a file.
    config_file_path = os.path.join(data_dir, 'config.ini')
    save_config(config, path=config_file_path)

    # Save spectral network data.
    for i, spectral_network in enumerate(spectral_networks):
        data_file_path = os.path.join(
            data_dir,
            'data_{}.json'.format(
                str(i).zfill(len(str(len(spectral_networks)-1)))
            )
        )
        logging.info('Saving data to {}.'.format(data_file_path))
        with open(data_file_path, 'wb') as fp:
            spectral_network.save_json_data(fp,)

    if make_zipped_file is True:
        file_list = [config_file_path]
        # Make a compressed data file.
        file_list += glob.glob(os.path.join(data_dir, 'data_*.json'))
        zipped_file_path = data_dir + '.zip'
        logging.info('Save compressed data to {}.'.format(zipped_file_path))
        with zipfile.ZipFile(zipped_file_path, 'w',
                             zipfile.ZIP_DEFLATED) as fp:
            for a_file in file_list:
                fp.write(a_file, os.path.relpath(a_file, data_dir))


def make_spectral_network_plot(config, spectral_networks, master=None,
                               show_plot=True, **kwargs):
    spectral_network_plot_title = 'Spectral Network'

    if matplotlib.rcParams['backend'] == 'TkAgg':
        spectral_network_plot = NetworkPlotTk(
            master=master,
            title=spectral_network_plot_title
        )
    else:
        spectral_network_plot = NetworkPlot(
            title=spectral_network_plot_title
        )

    for spectral_network in spectral_networks:
        logging.info('Generating the plot of a spectral network '
                     '@ theta = {}...'.format(spectral_network.phase))
        spectral_network_plot.draw(spectral_network)

    if show_plot is True:
        spectral_network_plot.show()

    #if master is None:
    #    try:
    #        raw_input('Press any key to continue...')
    #    except NameError:
    #        pass
    return spectral_network_plot


