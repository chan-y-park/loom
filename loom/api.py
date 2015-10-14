import time
import os
import glob
import zipfile
import logging
import Tkinter as tk
import tkFileDialog
import matplotlib
import pdb

from multiprocessing import Queue
from logutils_queue import QueueHandler
from config import LoomConfig
from trivialization import SWDataWithTrivialization
from spectral_network import SpectralNetwork
from parallel import parallel_get_spectral_network
from plotting import NetworkPlot, NetworkPlotTk


class SpectralNetworkData:
    """
    A container class of information relevant to
    a set of spectral networks generated from a 
    single Seiberg-Witten data.
    """
    def __init__(self, sw_data, spectral_networks):
        self.sw_data = sw_data
        self.spectral_networks = spectral_networks

def get_logging_formatter(level):
    if level == logging.DEBUG:
        logging_format = '%(module)s@%(lineno)d: %(funcName)s: %(message)s'
    elif level == logging.INFO:
        logging_format = '%(process)d: %(message)s'
    else:
        logging_format = '%(message)s'
    return logging.Formatter(logging_format)

def set_logging(level, queue=None, file_name='logs/log.loom.txt'):
    #print('Setting logging level to "{}".'.format(level))

    logger = logging.getLogger()
    # Remove other handlers.
    logger.handlers = []
    #for handler in handlers:
    #    logger.removeHandler(handler)

    logger.setLevel(level)
    formatter = get_logging_formatter(level)

    # Create a log file.
    fh = logging.FileHandler(file_name, 'w')
    fh.setLevel(level)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    if queue is None:
        # Create a stream handler to stderr.
        logger.addHandler(
            get_logging_handler(level, logging.StreamHandler, None)
        )
    else:
        # Create a queue handler for multiprocessing.
        logger.addHandler(
            get_logging_handler(level, QueueHandler, queue)
        )

def get_logging_handler(level, handler_class, buffer_object):
    h = handler_class(buffer_object)
    h.setLevel(level)
    h.setFormatter(get_logging_formatter(level))
    return h


def load_config(file_path=None):
#    if file_path is None:
#        toplevel = tk.Toplevel()
#        file_opts = {
#            'defaultextension': '.ini',
#            'initialdir': os.curdir,
#            'initialfile': 'config.ini',
#            'parent': toplevel,
#            'title': 'Select a configuration file to load.',
#        }
#        file_path = tkFileDialog.askopenfilename(**file_opts)
#        toplevel.destroy()
#        if file_path == '':
#            return None
    if file_path is None:
        return None
    config = LoomConfig()
    logging.info('Loading configuration from {}...'.format(file_path))
    config.read(file_path)
    logging.info('Finished loading configuration from {}.'.format(file_path))

    return config

def save_config(config, file_path=None):
#    if file_path is None:
#        toplevel = tk.Toplevel()
#        file_opts = {
#            'defaultextension': '.ini',
#            'initialdir': os.curdir,
#            'initialfile': 'config.ini',
#            'parent': toplevel,
#            'title': 'Save the current configuration to a file.',
#        }
#        file_path = tkFileDialog.asksaveasfilename(**file_opts)
#        toplevel.destroy()
#        if file_path == '':
#            return None
#
    if file_path is None:
        return None
    logging.info('Saving configuration to {}.'.format(file_path))
    with open(file_path, 'w') as fp:
        config.parser.write(fp)
    logging.info('Finished saving configuration to {}.'.format(file_path))

    return None

def load_spectral_network(
    data_dir=None,
    logging_queue=None,
    result_queue=None,
):
#    if data_dir is None:
#        return (None, None)
#        toplevel = tk.Toplevel()
#        dir_opts = {
#            'initialdir': os.curdir,
#            'mustexist': False,
#            'parent': toplevel,
#            'title': 'Select a directory that contains data files.',
#        }
#        data_dir = tkFileDialog.askdirectory(**dir_opts)
#        toplevel.destroy()
#        if data_dir == '':
#            return (None, None)
#
    if data_dir is None:
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
        spectral_network = SpectralNetwork()
        with open(data_file, 'r') as fp:
            spectral_network.set_from_json_data(fp)
            spectral_network_list.append(spectral_network)

    data = SpectralNetworkData(sw, spectral_network_list)
    logging.info('Finished loading data from {}.'.format(data_dir))

    rv = (config, data)
    if result_queue is None:
        return rv
    else:
        result_queue.put(rv)
        return None


def save_spectral_network(config, spectral_network_data, data_dir=None,
                          make_zipped_file=True):
    spectral_networks = spectral_network_data.spectral_networks
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
    save_config(config, file_path=config_file_path)

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

    logging.info('Finished saving data to {}.'.format(data_dir))
    return None


def generate_spectral_network(config, phase=None, result_queue=None,):
    """
    Generate one or more spectral networks according to
    the command-line options and the configuration file
    and return a list of data obtained from SpectralNetwork.get_data()
    """

    phase_range = config['phase_range']
    sw = SWDataWithTrivialization(config)

    start_time = time.time()
    logging.info('start cpu time: %s', start_time)

    if(phase is not None):
        logging.info('Generate a single spectral network at theta = {}.'
                     .format(phase))
        spectral_network = SpectralNetwork(
            phase=phase, 
        ) 

        spectral_network.grow(config, sw)

        spectral_networks = [spectral_network]

    elif(phase_range is not None):
        logging.info('Generate multiple spectral networks.')
        logging.info('phase_range = {}.'.format(phase_range))
        spectral_networks = parallel_get_spectral_network(
            sw, 
            config,
        ) 

    end_time = time.time()
    logging.info('end cpu time: %.8f', end_time)
    logging.info('elapsed cpu time: %.8f', end_time - start_time)

    rv = SpectralNetworkData(sw, spectral_networks)
    if result_queue is None:
        return rv
    else:
        result_queue.put(rv)
        return None

def make_spectral_network_plot(spectral_network_data, master=None,
                               show_plot=True, plot_range=None, **kwargs):
    sw_data = spectral_network_data.sw_data
    spectral_networks = spectral_network_data.spectral_networks
    spectral_network_plot_title = 'Spectral Network'

    if matplotlib.rcParams['backend'] == 'TkAgg':
        spectral_network_plot = NetworkPlotTk(
            master=master,
            title=spectral_network_plot_title,
            plot_range=plot_range,
        )
    else:
        spectral_network_plot = NetworkPlot(
            title=spectral_network_plot_title,
            plot_range=plot_range,
        )

    for spectral_network in spectral_networks:
        logging.info('Generating the plot of a spectral network '
                     '@ theta = {}...'.format(spectral_network.phase))
        spectral_network_plot.draw(
            spectral_network, 
            sw_data.branch_points,
            punctures=sw_data.punctures,
            irregular_singularities=sw_data.irregular_singularities,
            g_data=sw_data.g_data,
            **kwargs
        )

    if show_plot is True:
        spectral_network_plot.show()

    return spectral_network_plot
