import time
import os
import glob
import logging
import pdb

from config import LoomConfig
from geometry import SWData, get_ramification_points
from spectral_network import SpectralNetwork
from parallel import parallel_get_spectral_network


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


def load_spectral_network(data_dir):
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

    return spectral_network_list
