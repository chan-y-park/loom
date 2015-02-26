import time
import os
import glob
import zipfile, zlib
import logging
import pdb

from config import LoomConfig
from geometry import SWData, get_ramification_points
from spectral_network import SpectralNetwork
from parallel import parallel_get_spectral_network
from plotting import SpectralNetworkPlot


def generate_spectral_network(config, phase=None, show_plot=False,
                              plot_on_cylinder=False):
    """
    Generate one or more spectral networks according to
    the command-line options and the configuration file
    and return a list of data obtained from SpectralNetwork.get_data()
    """
    file_list = []

    #phase = opts['phase']
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

    # Plot spectral networks.
    if show_plot is True:
        spectral_network_plot = SpectralNetworkPlot(
            config,
            plot_on_cylinder,
            #plot_data_points=True,
            #plot_joints=True,
            #plot_bins=True,
            #plot_segments=True,
        )

        if (len(spectral_network_list) > 0):
            for spectral_network in spectral_network_list:
                spectral_network_plot.draw(spectral_network)

        spectral_network_plot.show()

    # Prepare to save spectral network data to files.
    timestamp = str(int(time.time()))
    data_save_dir = os.path.join(
        config['root_dir'], 
        config['data_dir'], 
        timestamp
    )

    logging.info('Make a directory {} to save data.'.format(data_save_dir))
    os.makedirs(data_save_dir)

    # Save configuration to a file.
    config_file_name = os.path.join(data_save_dir, 'config.ini')
    logging.info('Save configuration to {}.'.format(config_file_name))
    with open(config_file_name, 'wb') as fp:
        config.parser.write(fp)
        file_list.append(config_file_name)

    # Save spectral network data.
    for i, spectral_network in enumerate(spectral_network_list):
        data_file_name = os.path.join(
            data_save_dir,
            'data_{}.json'.format(
                str(i).zfill(len(str(len(spectral_network_list)-1)))
            )
        )
        logging.info('Saving data to {}.'.format(data_file_name))
        with open(data_file_name, 'wb') as fp:
            spectral_network.save_json_data(fp,)

    # Make a compressed data file.
    file_list += glob.glob(os.path.join(data_save_dir, 'data_*.json'))
    zipped_file_name = data_save_dir + '.zip'
    logging.info('Save compressed data to {}.'.format(zipped_file_name))
    with zipfile.ZipFile(zipped_file_name, 'w', zipfile.ZIP_DEFLATED) as fp:
        for a_file in file_list:
            fp.write(a_file, os.path.relpath(a_file, data_save_dir))

    return spectral_network_list


def load_spectral_network(data_dir, plot_on_cylinder=False):
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

    # Make plots from the loaded data
    spectral_network_plot = SpectralNetworkPlot(
        config,
        plot_on_cylinder=plot_on_cylinder,
    )
    for spectral_network in spectral_network_list:
        logging.info('Generating the plot of a spectral network '
                     '@ theta = {}...'.format(spectral_network.phase))
        spectral_network_plot.draw(spectral_network)
    spectral_network_plot.show()

    return spectral_network_list
