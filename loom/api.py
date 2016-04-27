import time
import os
import json
import glob
import zipfile
import logging
import subprocess
import matplotlib
# import traceback
import pdb

from logging.handlers import RotatingFileHandler
from logutils_queue import QueueHandler
from config import LoomConfig
from trivialization import SWDataWithTrivialization
from spectral_network import SpectralNetwork
from parallel import parallel_get_spectral_network
from plotting import NetworkPlot, NetworkPlotTk
from misc import get_phases_from_dict 
from misc import get_phase_dict 


class SpectralNetworkData:
    """
    A container class of information relevant to
    a set of spectral networks generated from a
    single Seiberg-Witten data.
    """
    def __init__(
        self,
        sw_data=None,
        spectral_networks=None,
        config=None,
        config_file_path=None,
        data_dir=None,
        logger_name='loom',
    ):
        self.logger_name = logger_name
        self.config = config
        self.sw_data = sw_data
        self.spectral_networks = spectral_networks
        self.data_attributes = ['sw_data', 'spectral_networks',]

        if config_file_path is not None:
            self.config = load_config(
                file_path=config_file_path,
                logger_name=self.logger_name,
            )
        elif data_dir is not None:
            load_spectral_network(
                data_dir=data_dir,
                spectral_network_data=self,
                logger_name=self.logger_name,
            )

    def save(self, data_dir=None):
        save_spectral_network(
            self.config,
            self,
            data_dir=data_dir,
            logger_name=self.logger_name,
        )

    def generate(
        self, phases=None, n_processes=0, extend=False,
        result_queue=None, logging_queue=None, cache_dir=None,
    ):
        logger = logging.getLogger(self.logger_name)
        if cache_dir is not None and os.path.exists(cache_dir) is False:
            os.makedirs(cache_dir)
       
        accuracy = self.config['accuracy']
        data_file_prefix = None

        try:
            if extend is False:
                data_file_prefix = 'data'
                start_time = time.time()
                logger.info(
                    'Started @ {}'.format(get_date_time_str(start_time))
                )
                self.sw_data = SWDataWithTrivialization(
                    self.config, logger_name=self.logger_name
                )
                if cache_dir is not None:
                    sw_data_file_path = os.path.join(cache_dir, 'sw_data.json')
                    self.sw_data.save(sw_data_file_path)

                if phases is None:
                    phase = get_phases_from_dict(
                        self.config['phase'], accuracy
                    )

            else:
                data_file_prefix = 'data_'
                # Extend spectral networks with additional phases.
                if phases is None:       
                    # No additional phase to add.
                    logger.warning(
                        'No additional phase given, '
                        'stop extending the spectral networks.'
                    )
                    return None
#                else:
#                    additional_phase_dict = get_phase_dict(phase)
#                    prev_phases = [sn.phase for sn in self.spectral_networks]
#                    prev_num_of_spectral_networks = len(prev_phases)
#                    add_config_phase(self.config, additional_phase_dict)
#                    new_phases = get_phases_from_dict(
#                        self.config['phase'], accuracy
#                    )
#                    phases = []
#                    for a_phase in new_phases:
#                        delta = [
#                            abs(a_phase - prev_phase) 
#                            for prev_phase in prev_phases
#                        ]
#                        if min(delta) > accuracy:
#                            phases.append(a_phase)

            if len(phases) == 0:
                logger.warning('No phase to generate.')
                spectral_networks = []

            elif(len(phases) == 1):
                phase = phases[0]
                logger.info('Generate a single spectral network at theta = {}.'
                            .format(phase))
                spectral_network = SpectralNetwork(
                    phase=phase,
                    logger_name=self.logger_name,
                )

                if cache_dir is not None:
                    cache_file_path = os.path.join(
                        cache_dir,
                        '{}_0.json'.format(data_file_prefix),
                    )
                else:
                    cache_file_path = None
                spectral_network.grow(
                    config=self.config, sw_data=self.sw_data,
                    cache_file_path=cache_file_path,
                )

                spectral_networks = [spectral_network]

            else:
#                phases = get_phases_from_dict(self.config['phase'], accuracy)
#                if extend is True:
#                    phases = [
#                        a_phase for a_phase in phases
#                        if (
#                            min(
#                                [abs(sn.phase - a_phase)
#                                 for sn in self.spectral_networks]
#                            ) > accuracy
#                        )
#                    ]
                logger.info('Generate multiple spectral networks.')
                logger.info('Number of phases: {}'.format(len(phases)))
                seed_spectral_networks = [
                    SpectralNetwork(
                        phase=a_phase,
                        logger_name=self.logger_name,
                    )
                    for a_phase in phases
                ]
                
                spectral_networks = parallel_get_spectral_network(
                    config=self.config,
                    sw_data=self.sw_data,
                    spectral_networks=seed_spectral_networks,
                    n_processes=n_processes,
                    logger_name=self.logger_name,
                    cache_dir=cache_dir,
                    data_file_prefix=data_file_prefix,
                )

        except (KeyboardInterrupt, SystemExit) as e:
            logger.warning(
                'SpectralNetworkData.generate() '
                'caught {} while generating spectral networks.'
                .format(type(e))
            )

        # TODO: handle and print all the other exceptions
        # to web UI.
        if extend is True:
            self.spectral_networks += spectral_networks
            self.spectral_networks.sort(key=lambda sn: sn.phase)
        else:
            self.spectral_networks = spectral_networks
            end_time = time.time()
            logger.info('Finished @ {}'.format(get_date_time_str(end_time)))
            logger.info('elapsed cpu time: %.3f', end_time - start_time)

        if logging_queue is not None:
            # Put a mark that generating spectral networks is done.
            try:
                logging_queue.put_nowait(None)
            except:
                logger.warn(
                    'Failed in putting a finish mark in the logging queue.'
                )

        if result_queue is not None:
            result_queue.put(self)

        if cache_dir is not None and extend is False:
            version_file_path = os.path.join(cache_dir, 'version')
            save_version(version_file_path)
            # NOTE: The following should be placed
            # at the last stage of spectral network generation. 
            config_file_path = os.path.join(cache_dir, 'config.ini')
            self.config.save(config_file_path)

    def extend(
        self,
        additional_n_steps=0,
        new_mass_limit=None,
        additional_iterations=0,
        additional_phases=None,
        n_processes=0,
        result_queue=None,
        logging_queue=None,
        cache_dir=None,
    ):
        logger = logging.getLogger(self.logger_name)
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)

        if additional_n_steps > 0:
            # Extend existing S-walls.
            self.config['num_of_steps'] += additional_n_steps

        if additional_iterations > 0:
            # Grow a spectral network for additional iterations.
            self.config['num_of_iterations'] += additional_iterations

#        if (
#            new_mass_limit is not None and 
#            new_mass_limit != self.config['mass_limit']
#        ):
#            if additional_n_steps == 0:
#                logger.warning(
#                    'Changing the mass limit without increasing '
#                    'number of steps has no effect.'
#                    'The mass limit will be kept to {}'
#                    .format(self.config['mass_limit'])
#                )
#            else:
#                self.config['mass_limit'] = new_mass_limit
            if new_mass_limit is not None:
                self.config['mass_limit'] = new_mass_limit
        
        start_time = time.time()
        logger.info('Started @ {}'.format(get_date_time_str(start_time)))

        # First rotate back the z-plane to the location
        # where spectral networks are generated.
        self.rotate_back()

        if (
            additional_n_steps > 0 or 
            additional_iterations > 0 or
            new_mass_limit is not None
        ):
            try:
                logger.info('Extending spectral networks...')
                if len(self.spectral_networks) == 1:
                    sn = self.spectral_networks[0]
                    sn.grow(
                        self.config,
                        self.sw_data,
                        additional_iterations=additional_iterations,
                        additional_n_steps=additional_n_steps,
                        new_mass_limit=new_mass_limit,
                    )
                    if cache_dir is not None:
                        sn_data_file_path = os.path.join(
                            cache_dir, 'data_0.json',
                        )
                        sn.save(sn_data_file_path)

                else:
                    self.spectral_networks = parallel_get_spectral_network(
                        config=self.config,
                        sw_data=self.sw_data,
                        spectral_networks=self.spectral_networks,
                        n_processes=n_processes,
                        additional_iterations=additional_iterations,
                        additional_n_steps=additional_n_steps,
                        new_mass_limit=new_mass_limit,
                        logger_name=self.logger_name,
                        cache_dir=cache_dir,
                    )
            except (KeyboardInterrupt, SystemExit) as e:
                logger.warning(
                    'SpectralNetworkData.extend() '
                    'caught {} while generating spectral networks.'
                    .format(type(e))
                )
        else:
            # No extension of exhisting spectral networks.
            if cache_dir is not None:
                # Save current spectral networks to cache.
                for i, sn in enumerate(self.spectral_networks):
                    cache_file_path = os.path.join(
                        cache_dir,
                        'data_{}.json'.format(i)
                    )
                    logger.info(
                        'Saving cache data to {}.'.format(cache_file_path)
                    )
                    sn.save(cache_file_path)

        if additional_phases is not None:
            additional_phase_dict = get_phase_dict(additional_phases)
            prev_phases = [sn.phase for sn in self.spectral_networks]
            add_config_phase(self.config, additional_phase_dict)
            new_phases = get_phases_from_dict(
                self.config['phase'], self.config['accuracy'],
            )
            phases = []
            for a_phase in new_phases:
                delta = [
                    abs(a_phase - prev_phase) 
                    for prev_phase in prev_phases
                ]
                if min(delta) > self.config['accuracy']:
                    phases.append(a_phase)

            logger.info(
                'Adding spectral networks with phase = {}'
                .format(additional_phases)
            )
            self.generate(
                phases=phases,
                n_processes=n_processes,
                extend=True,
                cache_dir=cache_dir,
            )

        end_time = time.time()
        logger.info('Finished @ {}'.format(get_date_time_str(end_time)))
        logger.info('elapsed cpu time: %.3f', end_time - start_time)

        if logging_queue is not None:
            # Put a mark that generating spectral networks is done.
            try:
                logging_queue.put_nowait(None)
            except:
                logger.warn(
                    'Failed in putting a finish mark in the logging queue.'
                )

        if result_queue is not None:
            result_queue.put(self)

        if cache_dir is not None:
            version_file_path = os.path.join(cache_dir, 'version')
            save_version(version_file_path)
            sw_data_file_path = os.path.join(cache_dir, 'sw_data.json')
            self.sw_data.save(sw_data_file_path)
            # NOTE: The following should be placed
            # at the last stage of spectral network generation. 
            config_file_path = os.path.join(cache_dir, 'config.ini')
            self.config.save(config_file_path)

    def plot(self, plot_range=None):
        # TODO: Implement if needed.
        return None

    def set_z_rotation(self, z_rotation):
        self.sw_data.set_z_rotation(z_rotation)
        for sn in self.spectral_networks:
            sn.set_z_rotation(z_rotation)

    def reset_z_rotation(self):
        z_r = self.sw_data.z_plane_rotation
        self.set_z_rotation(z_r)

    def rotate_back(self):
        bc_r = self.sw_data.branch_cut_rotation
        self.set_z_rotation(1/bc_r)


def get_logging_formatter(level):
    if level == logging.DEBUG:
        logging_format = (
            '%(process)d: %(module)s@%(lineno)d: %(funcName)s: %(message)s'
        )
    elif level == logging.INFO:
        logging_format = '%(process)d: %(message)s'
    elif level == logging.WARNING:
        logging_format = (
            '%(process)d: === WARNING ===\n'
            '%(process)d: %(funcName)s: %(message)s'
        )
    elif level == logging.ERROR:
        logging_format = (
            '%(process)d: ### ERROR ###\n'
            '%(process)d: %(funcName)s: %(message)s'
        )
    else:
        logging_format = '%(process)d: %(message)s'
    return logging.Formatter(logging_format)


class LoomLoggingFormatter(logging.Formatter):
    debug_format = (
        '%(process)d: %(module)s@%(lineno)d: '
        '%(funcName)s: %(message)s'
    )
    info_format = '%(process)d: %(message)s'
    warning_format = '%(process)d: === WARNING === %(funcName)s: %(message)s'
    error_format = '%(process)d: ### ERROR ### %(funcName)s: %(message)s'
    default_format = '%(message)s'

    def __init__(self, fmt=default_format):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):
        saved_format = self._fmt

        if record.levelno == logging.DEBUG:
            self._fmt = LoomLoggingFormatter.debug_format
        elif record.levelno == logging.INFO:
            self._fmt = LoomLoggingFormatter.info_format
        elif record.levelno == logging.WARNING:
            self._fmt = LoomLoggingFormatter.warning_format
        elif record.levelno == logging.ERROR:
            self._fmt = LoomLoggingFormatter.error_format

        result = logging.Formatter.format(self, record)

        self._fmt = saved_format

        return result


def set_logging(
    logger_name='loom',
    logging_level=None,
    logging_level_name='INFO',
    logging_queue=None,
    logging_stream=None,
    logging_file_name=None,
    remove_handlers=True,
    use_rotating_file_handler=False,
):
    #logging.basicConfig()
    logger = logging.getLogger(logger_name)
    if logging_level is None:
        logging_level = logging.getLevelName(logging_level_name)
    logger.setLevel(logging_level)
    #formatter = get_logging_formatter(logging_level)
    formatter = LoomLoggingFormatter()

    #stdout_h = logging.StreamHandler(sys.stdout)
    #stdout_h.setFormatter(formatter)
    #logging.root.addHandler(stdout_h)
    #logging.root.setLevel(logging_level)

    if remove_handlers is True:
        # Remove other handlers.
        logger.handlers = []

    if logging_file_name is not None:
        # Create a log file.
        if use_rotating_file_handler is True:
            fh = RotatingFileHandler(
                logging_file_name,
                mode='w',
                maxBytes=(10 * 1024 * 1024),
                backupCount=10,
            )
        else:
            fh = logging.FileHandler(logging_file_name, 'w')
        fh.setLevel(logging_level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    if logging_stream is not None:
        # Create a stream handler to logging_stream.
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
    f = LoomLoggingFormatter()
    h.setFormatter(f)
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
#    with open(file_path, 'w') as fp:
#        config.parser.write(fp)
    config.save(file_path)
    logger.info('Finished saving configuration to {}.'.format(file_path))

    return None


def load_spectral_network(
    data_dir=None,
    spectral_network_data=None,
    #logging_queue=None,
    #result_queue=None,
    logger_name='loom',
):
    if data_dir is None:
        return (None, None)

    logger = logging.getLogger(logger_name)
    logger.info('Opening data directory "{}"...'.format(data_dir))

    config = LoomConfig(logger_name=logger_name)
    with open(os.path.join(data_dir, 'config.ini')) as fp:
        config.read(fp)

    # Check the version of the saved data.
    try:
        current_version = get_current_branch_version()
        version_file_path = os.path.join(data_dir, 'version')
        with open(version_file_path, 'r') as fp:
            data_version = fp.read()
    except:
        data_version = None

    if data_version is not None and current_version != data_version:
        logger.debug('The version of the data is different '
                     'from the current version of loom.')

    sw_data_file_path = os.path.join(data_dir, 'sw_data.json')
    with open(sw_data_file_path, 'r') as fp:
        json_data = json.load(fp)
        sw_data = SWDataWithTrivialization(config, logger_name=logger_name,
                                           json_data=json_data,)

    spectral_networks = []

    data_file_list = glob.glob(os.path.join(data_dir, 'data_*.json'))
#    data_file_list.sort()
    for data_file in data_file_list:
        logger.info('Loading {}...'.format(data_file))
        spectral_network = SpectralNetwork(logger_name=logger_name)
#        with open(data_file, 'r') as fp:
#            json_data = json.load(fp)
#            spectral_network.set_from_json_data(json_data, sw_data)
        spectral_network.load(data_file, sw_data)
        spectral_networks.append(spectral_network)
    spectral_networks.sort(key=lambda sn: sn.phase)

    logger.info('Finished loading data from {}.'.format(data_dir))

    if spectral_network_data is not None:
        spectral_network_data.config = config
        spectral_network_data.sw_data = sw_data
        spectral_network_data.spectral_networks = spectral_networks
        return None

    data = SpectralNetworkData(
        sw_data=sw_data,
        spectral_networks=spectral_networks,
    )
    return (config, data)
    #rv = (config, data)
    #if result_queue is None:
    #    return rv
    #else:
    #    result_queue.put(rv)
    #    return None


def save_spectral_network(
    config,
    spectral_network_data,
    data_dir=None,
    make_zipped_file=False,
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
    save_version(version_file_path)

    # Save configuration to a file.
    config_file_path = os.path.join(data_dir, 'config.ini')
    save_config(config, file_path=config_file_path, logger_name=logger_name)

    # Save geometric & trivialization data.
    sw_data_file_path = os.path.join(data_dir, 'sw_data.json')
    logger.info('Saving data to {}.'.format(sw_data_file_path))
#    with open(sw_data_file_path, 'wb') as fp:
#        json_data = sw_data.get_json_data()
#        json.dump(json_data, fp,)
    sw_data.save(sw_data_file_path)

    # Save spectral network data.
    for i, spectral_network in enumerate(spectral_networks):
        data_file_path = os.path.join(
            data_dir,
            'data_{}.json'.format(
                str(i).zfill(len(str(len(spectral_networks) - 1)))
            )
        )
        logger.info('Saving data to {}.'.format(data_file_path))
#        with open(data_file_path, 'wb') as fp:
#            json_data = spectral_network.get_json_data()
#            json.dump(json_data, fp,)
        spectral_network.save(data_file_path)

    if make_zipped_file is True:
        # Make a compressed data file.
        file_list = [config_file_path, sw_data_file_path]
        file_list += glob.glob(os.path.join(data_dir, 'data_*.json'))
        # TODO: the following routine for getting zipped_file_name has a bug.
        zipped_file_name = os.path.basename(os.path.normpath(data_dir))
        zipped_file_path = data_dir + '{}.zip'.format(zipped_file_name)
        logger.info('Save compressed data to {}.'.format(zipped_file_path))
        with zipfile.ZipFile(zipped_file_path, 'w',
                             zipfile.ZIP_DEFLATED) as fp:
            for a_file in file_list:
                fp.write(a_file, os.path.relpath(a_file, data_dir))

    logger.info('Finished saving data to {}.'.format(data_dir))
    return None


# XXX: generate_spectral_network() will be deprecated.
# Use SpectralNetworkData.generate().
def generate_spectral_network(
    config,
    phase=None,
    n_processes=0,
    #result_queue=None,
    #logging_queue=None,
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

    start_time = time.time()
    logger.info('Started @ {}'.format(get_date_time_str(start_time)))

    try:
        sw_data = SWDataWithTrivialization(config, logger_name=logger_name)

        if(isinstance(phase, float)):
            logger.info('Generate a single spectral network at theta = {}.'
                        .format(phase))
            spectral_network = SpectralNetwork(
                phase=phase,
                logger_name=logger_name,
            )

            spectral_network.grow(config, sw_data)

            spectral_networks = [spectral_network]

        elif(isinstance(phase, list)):
            logger.info('Generate multiple spectral networks.')
            logger.info('phases = {}.'.format(phase))
            spectral_networks = parallel_get_spectral_network(
                config=config,
                sw_data=sw_data,
                n_processes=n_processes,
                logger_name=logger_name,
            )

    except (KeyboardInterrupt, SystemExit) as e:
        logger.warning(
            'loom.api.generate_spectral_network() '
            'caught {} while generating spectral networks.'
            .format(type(e))
        )
        sw_data = None
        spectral_networks = None

    end_time = time.time()
    logger.info('Finished @ {}'.format(get_date_time_str(end_time)))
    logger.info('elapsed cpu time: %.3f', end_time - start_time)

    spectral_network_data = SpectralNetworkData(
        sw_data=sw_data, spectral_networks=spectral_networks,
    )

    #if logging_queue is not None:
    #    # Put a mark that generating spectral networks is done.
    #    try:
    #        logging_queue.put_nowait(None)
    #    except:
    #        logger.warn(
    #            'Failed in putting a finish mark in the logging queue.'
    #        )

    #if result_queue is not None:
    #    rv = (config, spectral_network_data)
    #    result_queue.put(rv)
    #    return None
    #else:
    #    return spectral_network_data

    return spectral_network_data


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

    # Rotate the z-plane into the location defined by the curve.
    spectral_network_data.reset_z_rotation()

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
            branch_cut_rotation=sw_data.branch_cut_rotation,
            **kwargs
        )
        logger.info(plot_legend)

    # Set the z-plane rotation back.
    # TODO: Decide whether to save a rotated data or a raw data.
    #spectral_network_data.set_z_rotation(1/z_r)

    if show_plot is True:
        spectral_network_plot.show()

    return spectral_network_plot


def get_date_time_str(a_time):
    a_date_time_str = (
        '{:02}-{:02}-{:02} {:02}:{:02}:{:02}'
        .format(*time.localtime(a_time)[:6])
    )
    return a_date_time_str


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


def save_version(file_path):
    try:
        version = get_current_branch_version()
    except:
        version = "None"
    with open(file_path, 'w') as fp:
        fp.write(version)


def add_config_phase(config, phase_dict):
    """
    Add to config['phase'] a phase dict.
    """
    config_phase = config['phase']
    config_phase['single'] += phase_dict['single']
    config_phase['range'] += phase_dict['range']

    config['phase'] = config_phase

