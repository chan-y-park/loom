import time
import sys
import os
import glob
import zipfile, zlib
import logging
import Tkinter as tk
import tkFileDialog
import pdb

from config import LoomConfig
from geometry import SWData, get_ramification_points
from spectral_network import SpectralNetwork
from parallel import parallel_get_spectral_network
from plotting import SpectralNetworkPlot


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

    logging.basicConfig(level=logging_level, format=logging_format, 
                        stream=sys.stdout)


def generate_spectral_network(config, phase=None):
    """
    Generate one or more spectral networks according to
    the command-line options and the configuration file
    and return a list of data obtained from SpectralNetwork.get_data()
    """

    phase_range = config['phase_range']
    sw = SWData(config)
    ramification_points = get_ramification_points(sw, config['accuracy'])

    start_time = time.time()
    logging.info('start cpu time: %s', start_time)

    if(phase is not None):
        logging.info('Generate a single spectral network at theta = {}.'
                     .format(phase))
        spectral_network = SpectralNetwork(
            phase=phase, 
            ramification_points=ramification_points,
            config=config
        ) 

        spectral_network.grow(sw, config)

        spectral_network_list = [spectral_network]

    elif(phase_range is not None):
        logging.info('Generate multiple spectral networks.')
        logging.info('phase_range = {}.'.format(phase_range))
        spectral_network_list = parallel_get_spectral_network(
            sw, ramification_points, config
        ) 

    end_time = time.time()
    logging.info('end cpu time: %.8f', end_time)
    logging.info('elapsed cpu time: %.8f', end_time - start_time)

    return spectral_network_list


def load_config(path=None):
    if path is None:
        root = tk.Tk()
        file_opts = {
            'defaultextension': '.ini',
            'initialdir': os.curdir,
            'initialfile': 'config.ini',
            'parent': root,
            'title': 'Select a configuration file to load.',
        }
        path = tkFileDialog.askopenfilename(**file_opts)
        root.destroy()
        if path == '':
            return None

    config = LoomConfig()
    config.read(path)

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

    sw = SWData(config)
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
                               **kwargs):
    if master is None:
        master=tk.Tk()
        master.withdraw()

    spectral_network_plot = SpectralNetworkPlot(
        master=master,
        config=config,
        **kwargs
    )

    for spectral_network in spectral_networks:
        logging.info('Generating the plot of a spectral network '
                     '@ theta = {}...'.format(spectral_network.phase))
        spectral_network_plot.draw(spectral_network)

    spectral_network_plot.show()

    return master 

