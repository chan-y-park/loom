import time
import os
import json
import glob
import zipfile
import logging
import subprocess
import signal
import matplotlib
import pdb

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


def set_logging(
    logger_name='loom',
    logging_level=logging.INFO, 
    logging_queue=None, 
    logging_stream=None,
    logging_file_name=None,
    remove_handlers=True,
    use_rotating_file_handler=False,
):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging_level)
    formatter = get_logging_formatter(logging_level)

    if remove_handlers is True:
        # Remove other handlers.
        logger.handlers = []

    if logging_file_name is not None:
        # Create a log file.
        if use_rotating_file_handler is True:
            fh = logging.handlers.RotatingFileHandler(
                logging_file_name,
                mode='w',
                maxBytes=10*1024*1024,
                backupCount=10,
            )
        else:
            fh = logging.FileHandler(logging_file_name, 'w')
        fh.setLevel(logging_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    if logging_stream is not None:
        # Create a stream handler to stderr.
        logger.addHandler(
            get_logging_handler(
                logging_level, logging.StreamHandler, logging_stream
            )
        )
    if logging_queue is not None:
        # Create a queue handler for multiprocessing.
        logger.addHandler(
            get_logging_handler(
                logging_level, QueueHandler, logging_queue
            )
        )


def get_logging_handler(level, handler_class, buffer_object):
    h = handler_class(buffer_object)
    h.setLevel(level)
    h.setFormatter(get_logging_formatter(level))
    return h


def load_config(file_path=None, logger_name='loom',):
    logger = logging.getLogger(logger_name)
    if file_path is None:
        return None
    config = LoomConfig(logger_name=logger_name)
    logger.info('Loading configuration from {}...'.format(file_path))
    with open(file_path, 'r') as fp:
        config.read(fp)
    logger.info('Finished loading configuration from {}.'.format(file_path))

    return config


def save_config(config, file_path=None, logger_name='loom',):
    logger = logging.getLogger(logger_name)
    if file_path is None:
        return None
    logger.info('Saving configuration to {}.'.format(file_path))
    with open(file_path, 'w') as fp:
        config.parser.write(fp)
    logger.info('Finished saving configuration to {}.'.format(file_path))

    return None


def load_spectral_network(
    data_dir=None,
    logging_queue=None,
    result_queue=None,
    logger_name='loom',
):
    logger = logging.getLogger(logger_name)
    

    if data_dir is None:
        return (None, None)
    logger.info('Opening data directory "{}"...'.format(data_dir))

    config = LoomConfig(logger_name=logger_name)
    with open(os.path.join(data_dir, 'config.ini')) as fp:
        config.read(fp)

    # Check the version of the saved data.
    current_version = get_current_branch_version()
    version_file_path = os.path.join(data_dir, 'version')
    try:
        with open(version_file_path, 'r') as fp:
            data_version = fp.read()
    except IOError:
        data_version = None
    if current_version != data_version:
        logger.debug('The version of the data is different '
                     'from the current version of loom.')

    sw_data_file_path = os.path.join(data_dir, 'sw_data.json')
    with open(sw_data_file_path, 'r') as fp:
        json_data = json.load(fp)
        sw_data = SWDataWithTrivialization(config, logger_name=logger_name,
                                           json_data=json_data,)

    spectral_networks = []

    data_file_list = glob.glob(os.path.join(data_dir, 'data_*.json'))
    data_file_list.sort()
    for data_file in data_file_list:
        logger.info('Loading {}...'.format(data_file))
        spectral_network = SpectralNetwork(logger_name=logger_name)
        with open(data_file, 'r') as fp:
            json_data = json.load(fp)
            spectral_network.set_from_json_data(json_data, sw_data)
            spectral_networks.append(spectral_network)

    data = SpectralNetworkData(sw_data, spectral_networks)
    logger.info('Finished loading data from {}.'.format(data_dir))

    rv = (config, data)
    if result_queue is None:
        return rv
    else:
        result_queue.put(rv)
        return None


def save_spectral_network(
    config, 
    spectral_network_data,
    data_dir=None,
    make_zipped_file=True,
    logger_name='loom',
):
    logger = logging.getLogger(logger_name)

    sw_data = spectral_network_data.sw_data
    spectral_networks = spectral_network_data.spectral_networks

    if data_dir is None:
        # Prepare to save spectral network data to files.
        timestamp = str(int(time.time()))
        data_dir = os.path.join(
            get_loom_dir(),
            'data',
            timestamp
        )

    logger.info('Make a directory {} to save data.'.format(data_dir))
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # Save version data.
    version_file_path = os.path.join(data_dir, 'version')
    version = get_current_branch_version() 
    with open(version_file_path, 'w') as fp:
        fp.write(version)

    # Save configuration to a file.
    config_file_path = os.path.join(data_dir, 'config.ini')
    save_config(config, file_path=config_file_path, logger_name=logger_name)

    # Save geometric & trivialization data.
    sw_data_file_path = os.path.join(data_dir, 'sw_data.json')
    logger.info('Saving data to {}.'.format(sw_data_file_path))
    with open(sw_data_file_path, 'wb') as fp:
        json_data = sw_data.get_json_data()
        json.dump(json_data, fp,)

    # Save spectral network data.
    for i, spectral_network in enumerate(spectral_networks):
        data_file_path = os.path.join(
            data_dir,
            'data_{}.json'.format(
                str(i).zfill(len(str(len(spectral_networks) - 1)))
            )
        )
        logger.info('Saving data to {}.'.format(data_file_path))
        with open(data_file_path, 'wb') as fp:
            json_data = spectral_network.get_json_data()
            json.dump(json_data, fp,)

    if make_zipped_file is True:
        # Make a compressed data file.
        file_list = [config_file_path, sw_data_file_path]
        file_list += glob.glob(os.path.join(data_dir, 'data_*.json'))
        zipped_file_name = os.path.basename(os.path.normpath(data_dir)) 
        zipped_file_path = data_dir + '{}.zip'.format(zipped_file_name)
        logger.info('Save compressed data to {}.'.format(zipped_file_path))
        with zipfile.ZipFile(zipped_file_path, 'w',
                             zipfile.ZIP_DEFLATED) as fp:
            for a_file in file_list:
                fp.write(a_file, os.path.relpath(a_file, data_dir))

    logger.info('Finished saving data to {}.'.format(data_dir))
    return None


def stop_signal_handler(signum, frame):
    if signum == signal.SIGTERM:
        raise SystemExit('SIGTERM catched by generate_spectral_network.')


def generate_spectral_network(
    config,
    phase=None,
    n_processes=0,
    result_queue=None,
    logging_queue=None,
    logger_name='loom',
):
    """
    Generate one or more spectral networks according to
    the command-line options and the configuration file
    and return a list of data obtained from SpectralNetwork.get_data()
    """

    logger = logging.getLogger(logger_name)

    if phase is None:
        phase = config['phase']
    else:
        config['phase'] = phase
    try:
        start_time = time.time()
        start_date_time = (
            '{:02}-{:02}-{:02} {:02}:{:02}:{:02}'
            .format(*time.localtime(start_time)[:6])
        )
        logger.info('Started @ {}'.format(start_date_time))

        sw = SWDataWithTrivialization(config, logger_name=logger_name)

        if(isinstance(phase, float)):
            logger.info('Generate a single spectral network at theta = {}.'
                        .format(phase))
            spectral_network = SpectralNetwork(
                phase=phase, 
                logger_name=logger_name,
            ) 

            spectral_network.grow(config, sw)

            spectral_networks = [spectral_network]

        elif(isinstance(phase, list)):
            logger.info('Generate multiple spectral networks.')
            logger.info('phases = {}.'.format(phase))
            spectral_networks = parallel_get_spectral_network(
                sw, 
                config,
                n_processes,
                logger_name=logger_name,
            ) 

        end_time = time.time()
        end_date_time = (
            '{:02}-{:02}-{:02} {:02}:{:02}:{:02}'
            .format(*time.localtime(end_time)[:6])
        )
        logger.info('Finished @ {}'.format(end_date_time))
        logger.info('elapsed cpu time: %.3f', end_time - start_time)

    except (KeyboardInterrupt, SystemExit) as e:
        logger.warning('loom.api caught {} while generating spectral networks.'
                       .format(type(e)))
        sw = None
        spectral_networks = None

    spectral_network_data = SpectralNetworkData(sw, spectral_networks)
    if logging_queue is not None:
        # Put a mark that generating spectral networks is done.
        try:
            logging_queue.put_nowait(None)
        except:
            logger.warn("Failed in putting a finish mark "
                        "in the logging queue.")

    if result_queue is None:
        return spectral_network_data
    else:
        rv = (config, spectral_network_data)
        result_queue.put(rv)
        return None


def make_spectral_network_plot(
    spectral_network_data,
    master=None,
    show_plot=True,
    plot_range=None,
    logger_name='loom',
    **kwargs
):
    logger = logging.getLogger(logger_name)
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
        logger.info('Generating the plot of a spectral network '
                    '@ theta = {}...'.format(spectral_network.phase))
        plot_legend = spectral_network_plot.draw(
            spectral_network, 
            sw_data.branch_points,
            punctures=(sw_data.regular_punctures + 
                       sw_data.irregular_punctures),
            irregular_singularities=sw_data.irregular_singularities,
            g_data=sw_data.g_data,
            **kwargs
        )
        
        logger.info(plot_legend)

    if show_plot is True:
        spectral_network_plot.show()

    return spectral_network_plot


def get_loom_dir():
    loom_dir = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        '..',
    )
    return os.path.realpath(loom_dir)


def get_current_branch_version():
    version = subprocess.check_output(
        ['git', '-C', get_loom_dir(), 'rev-parse', 'HEAD']
    )

    return version
