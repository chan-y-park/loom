import time
import os
import json
import glob
import sympy
import zipfile
import logging
import subprocess
import matplotlib
import mpldatacursor
# import traceback
# import pdb 

from logging.handlers import RotatingFileHandler
from matplotlib import pyplot
from sympy import oo

from logutils_queue import QueueHandler
from config import LoomConfig
from trivialization import SWDataWithTrivialization
from spectral_network import SpectralNetwork
from parallel import parallel_get_spectral_network
# TODO: plotting.py will be deprecated; use plot_ui.py
from plotting import NetworkPlot, NetworkPlotTk
from misc import get_phases_from_dict
from misc import get_phase_dict
from misc import parse_sym_dict_str


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
        self.data_attributes = ['sw_data', 'spectral_networks']

        if config_file_path is not None:
            self.config = LoomConfig(
                file_path=config_file_path,
                logger_name=self.logger_name,
            )
        elif data_dir is not None:
            self.load(data_dir=data_dir)

    def load(
        self, data_dir=None, logger_name=None,
        result_queue=None, logging_queue=None,
    ):
        if data_dir is None:
            return None

        if logger_name is None:
            logger_name = self.logger_name
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

        spectral_networks = []

        if config['parameter_sequence'] is None:
            # This is a standard family of networks
            # with different phases.
            sw_data_file_path = os.path.join(data_dir, 'sw_data.json')
            with open(sw_data_file_path, 'r') as fp:
                json_data = json.load(fp)
                sw_data = SWDataWithTrivialization(
                    config, logger_name=logger_name,
                    json_data=json_data,
                )

            data_file_list = glob.glob(os.path.join(data_dir, 'data_*.json'))
            for data_file in data_file_list:
                logger.info('Loading {}...'.format(data_file))
                spectral_network = SpectralNetwork(logger_name=logger_name)
                spectral_network.load(data_file, sw_data)
                spectral_networks.append(spectral_network)
            spectral_networks.sort(key=lambda sn: sn.phase)

            logger.info('Finished loading data from {}.'.format(data_dir))

        else:
            # This is a family of networks with different parameters.
            var_str, val_str = parse_sym_dict_str(
                config['parameter_sequence'],
                multi_parameter=True
            )
            n_steps = sympy.sympify(val_str)[2]
            sw_data = []
            for i in range(n_steps):
                sw_data_file_path = os.path.join(
                    data_dir, 'sw_data_'+str(i)+'.json'
                )
                with open(sw_data_file_path, 'r') as fp:
                    json_data = json.load(fp)
                    sw_data_i = SWDataWithTrivialization(
                        config, logger_name=logger_name,
                        json_data=json_data,
                    )
                sw_data.append(sw_data_i)

                # Now load the spectral network for the i-th value of the 
                # parameter
                data_file_i = os.path.join(
                    data_dir, 'data_' + str(i) + '.json'
                )
                logger.info('Loading {}...'.format(data_file_i))
                spectral_network = SpectralNetwork(logger_name=logger_name)
                spectral_network.load(data_file_i, sw_data_i)
                spectral_networks.append(spectral_network)

                logger.info('Finished loading data from {}.'.format(data_dir))

        self.config = config
        self.sw_data = sw_data
        self.spectral_networks = spectral_networks

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

    def save(
        self, data_dir=None, make_zipped_file=False,
        logger_name=None,
    ):
        if logger_name is None:
            logger_name = self.logger_name
        logger = logging.getLogger(logger_name)

        sw_data = self.sw_data
        spectral_networks = self.spectral_networks

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
        self.config.save(
            file_path=config_file_path,
            logger_name=logger_name
        )

        # Save geometric & trivialization data.
        if self.config.data['parameter_sequence'] is None:
            sw_data_file_path = os.path.join(data_dir, 'sw_data.json')
            logger.info('Saving data to {}.'.format(sw_data_file_path))
            sw_data.save(sw_data_file_path)
        else:
            for i, sw_data_i in enumerate(sw_data):
                sw_data_file_path = os.path.join(
                    data_dir, 'sw_data_'+str(i)+'.json'
                )
                logger.info('Saving data to {}.'.format(sw_data_file_path))
                sw_data_i.save(sw_data_file_path)

        # Save spectral network data.
        if self.config.data['parameter_sequence'] is None:
            for i, spectral_network in enumerate(spectral_networks):
                data_file_path = os.path.join(
                    data_dir,
                    'data_{}.json'.format(
                        str(i).zfill(len(str(len(spectral_networks) - 1)))
                    )
                )
                logger.info('Saving data to {}.'.format(data_file_path))
                spectral_network.save(data_file_path)
        else:
            for i, spectral_network in enumerate(spectral_networks):
                data_file_path = os.path.join(
                    data_dir,
                    'data_{}.json'.format(
                        str(i).zfill(len(str(len(spectral_networks) - 1)))
                    )
                )
                logger.info('Saving data to {}.'.format(data_file_path))
                spectral_network.save(data_file_path)


        if make_zipped_file is True:
            # Make a compressed data file.
            file_list = [config_file_path, sw_data_file_path]
            file_list += glob.glob(os.path.join(data_dir, 'data_*.json'))
            # TODO: the following routine for getting zipped_file_name
            # has a bug.
            zipped_file_name = os.path.basename(os.path.normpath(data_dir))
            zipped_file_path = data_dir + '{}.zip'.format(zipped_file_name)
            logger.info('Save compressed data to {}.'.format(zipped_file_path))
            with zipfile.ZipFile(zipped_file_path, 'w',
                                 zipfile.ZIP_DEFLATED) as fp:
                for a_file in file_list:
                    fp.write(a_file, os.path.relpath(a_file, data_dir))

        logger.info('Finished saving data to {}.'.format(data_dir))
        return None

    def generate(
        self, phases=None, n_processes=0, extend=False,
        result_queue=None, logging_queue=None, cache_dir=None,
    ):
        logger = logging.getLogger(self.logger_name)
        if cache_dir is not None and os.path.exists(cache_dir) is False:
            os.makedirs(cache_dir)

        accuracy = self.config['accuracy']
        data_file_prefix = None

        if self.config.data['parameter_sequence'] is None:
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
                        sw_data_file_path = os.path.join(
                            cache_dir, 'sw_data.json'
                        )
                        self.sw_data.save(sw_data_file_path)

                    if phases is None:
                        self.config['phase'] = get_phase_dict(
                            self.config['phase']
                        )
                        phases = get_phases_from_dict(
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

                if type(phases) == dict:
                    ph_from_dict = get_phases_from_dict(phases, accuracy)
                    phases = ph_from_dict

                if len(phases) == 0:
                    logger.warning('No phase to generate.')
                    spectral_networks = []

                elif(len(phases) == 1):
                    phase = phases[0]
                    logger.info(
                        'Generate a single spectral network at theta = {}.'
                        .format(phase)
                    )
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

        else:
            phase = 1.570795
            logger.info(
                'Study multiple parameters at fixed phase {}'
                .format(phase)
            )

            original_parameters = self.config['differential_parameters']
            spectral_networks = []
            sw_data_sequence = []
            parameter_sequence = {}
            var_str, val_str = parse_sym_dict_str(
                self.config['parameter_sequence'],
                multi_parameter=True
            )

            if (
                self.config['ramification_point_finding_method']
                =='from_branch_points'
            ):
                branch_points_sequence = sympy.sympify(
                self.config['branch_points_sequence']
            )
                branch_points_sequence_str = map(str, branch_points_sequence)

            if extend is True:
                raise NotImplementedError

            val = sympy.sympify(val_str)
            val_i = val[0]
            val_f = val[1]
            n_steps = val[2]
            step = (val_f - val_i) / (n_steps - 1)
            val_list = [val_i + j * step for j in range(n_steps)]

            logger.info(
                'will use these parameter values {}'.format(val_list)
            )

            for i, val_i in enumerate(val_list):
                logger.info(
                    'Producing Spectral Network for {} = {}'
                    .format(var_str, val_i)
                )
                # First insert the parameter value from the 
                # parameter sequence
                diff_params = parse_sym_dict_str(
                    original_parameters
                )
                diff_param_str = '{ '
                for pair in diff_params:
                    diff_param_str = (
                        diff_param_str + pair[0]+' : '+pair[1] +', '
                    )
                diff_param_str = (
                    diff_param_str + var_str +' : ' + str(val_i) + ' }'
                )
                self.config['differential_parameters'] = diff_param_str
                # Then, if branch points are given, also insert the 
                # branch point
                if (
                    self.config['ramification_point_finding_method']
                    =='from_branch_points'
                ):
                    self.config['branch_points'] = (
                        branch_points_sequence_str[i]
                    )

                data_file_prefix = 'multi_data_'+str(i)+'_'
                start_time = time.time()
                logger.info(
                    'Started @ {}'.format(get_date_time_str(start_time))
                )
                self.sw_data = SWDataWithTrivialization(
                    self.config, logger_name=self.logger_name
                )
                sw_data_sequence.append(self.sw_data)

                if cache_dir is not None:
                    sw_data_file_path = os.path.join(
                        cache_dir, 'sw_data' + str(i) + '.json'
                    )
                    self.sw_data.save(sw_data_file_path)

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

                spectral_networks.append(spectral_network)
                
            self.spectral_networks = spectral_networks
            self.sw_data = sw_data_sequence
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
        if cache_dir is not None and os.path.exists(cache_dir) is not True:
            os.makedirs(cache_dir)

        if additional_n_steps > 0:
            # Extend existing S-walls.
            self.config['num_of_steps'] += additional_n_steps

        if additional_iterations > 0:
            # Grow a spectral network for additional iterations.
            self.config['num_of_iterations'] += additional_iterations

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
            self.config['phase'] = get_phase_dict(self.config['phase'])
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

    def find_two_way_streets(self, search_radius=None,):
        soliton_tree_data = []
        for sn in self.spectral_networks:
            soliton_tree_data.append(
                sn.find_two_way_streets(
                    config=self.config,
                    sw_data=self.sw_data,
                    search_radius=search_radius,
                )
            )
        return soliton_tree_data

    def set_z_rotation(self, z_rotation):
        self.sw_data.set_z_rotation(z_rotation)
        for sn in self.spectral_networks:
            sn.set_z_rotation(z_rotation)

    def reset_z_rotation(self):
        # TODO: handle z-rotation for multi-parameter case
        # note that each parameter choice may have its own
        # different z_r
        if type(self.sw_data) is list:
            pass
        else:
            z_r = self.sw_data.z_plane_rotation
            self.set_z_rotation(z_r)

    def rotate_back(self):
        # TODO: handle z-rotation for multi-parameter case
        # note that each parameter choice may have its own
        # different z_r
        if type(self.sw_data) is list:
            pass
        else:
            bc_r = self.sw_data.branch_cut_rotation
            self.set_z_rotation(1/bc_r)

    def plot_s_walls(
        self,
        walls,
        plot_range=None,
        plot_data_points=False,
        file_path=None,
    ):
        """
        Plot walls on a complex plane for debugging purpose.
        """
        pyplot.figure()
        pyplot.axes().set_aspect('equal')

        sw_data = self.sw_data

        for pp in (sw_data.regular_punctures + sw_data.irregular_punctures):
            if pp.z == oo:
                continue
            pyplot.plot(
                pp.z.real, pp.z.imag, 'o',
                markeredgewidth=2, markersize=8,
                color='k', markerfacecolor='none',
                label=pp.label,
            )

        for bp in sw_data.branch_points:
            pyplot.plot(
                bp.z.real, bp.z.imag, 'x',
                markeredgewidth=2, markersize=8,
                color='k', label=bp.label,
            )

        for wall in walls:
            xs = wall.z.real
            ys = wall.z.imag
            pyplot.plot(xs, ys, '-', label=wall.label)

            if(plot_data_points is True):
                pyplot.plot(xs, ys, 'o', color='k', markersize=4)

        if plot_range is None:
            pyplot.autoscale(enable=True, axis='both', tight=True)
            pyplot.margins(0.1)
        else:
            [[x_min, x_max], [y_min, y_max]] = plot_range
            pyplot.xlim(x_min, x_max)
            pyplot.ylim(y_min, y_max)

        if file_path is not None:
            pyplot.savefig(file_path)
        else:
            mpldatacursor.datacursor(
                formatter='{label}'.format,
                hover=True,
            )

            pyplot.show()


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
    logger = logging.getLogger(logger_name)
    if logging_level is None:
        logging_level = logging.getLevelName(logging_level_name)
    logger.setLevel(logging_level)
    formatter = LoomLoggingFormatter()

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


# XXX: save_config() will be deprecated.
# Use LoomConfig.__init__().
def load_config(file_path=None, logger_name='loom',):
    config = LoomConfig(file_path=file_path, logger_name=logger_name)

    return config


# XXX: save_config() will be deprecated.
# Use LoomConfig.save().
def save_config(config, file_path=None, logger_name='loom',):
    config.save(file_path=file_path, logger_name=logger_name)

    return None


# XXX: load_spectral_network() will be deprecated.
# Use SpectralNetworkData.load().
def load_spectral_network(
    data_dir=None,
    result_queue=None,
    logger_name='loom',
):
    data = SpectralNetworkData(
        data_dir=data_dir,
        logger_name=logger_name,
    )
    if result_queue is not None:
        result_queue.put(data)
    else:
        return (data.config, data)


# XXX: save_spectral_network() will be deprecated.
# Use SpectralNetworkData.save().
def save_spectral_network(
    config,
    spectral_network_data,
    data_dir=None,
    make_zipped_file=False,
    logger_name='loom',
):
    spectral_network_data.save(
        data_dir=data_dir,
        make_zipped_file=make_zipped_file,
        logger_name=logger_name,
    )


# XXX: generate_spectral_network() will be deprecated.
# Use SpectralNetworkData.generate().
# PL: not clear what the new workflow will be after these deprecations,
# please give some directions
def generate_spectral_network(
    config,
    n_processes=0,
    result_queue=None,
    logger_name='loom',
):
    """
    Generate one or more spectral networks according to
    the command-line options and the configuration file
    and return a list of data obtained from SpectralNetwork.get_data()
    """

    spectral_network_data = SpectralNetworkData(
        config=config,
        logger_name=logger_name,
    )

    spectral_network_data.generate(
        n_processes=n_processes,
        result_queue=result_queue,
    )

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
    if type(spectral_network_data.sw_data) != list:
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
        # TODO: When using plot_ui.py, remove the following
        # and uncomment the next following lines.
        plot_legend = spectral_network_plot.draw(
            spectral_network,
            sw_data.branch_points,
            punctures=(sw_data.regular_punctures +
                       sw_data.irregular_punctures),
            irregular_singularities=sw_data.irregular_singularities,
            g_data=sw_data.g_data,
            branch_cut_rotation=sw_data.branch_cut_rotation,
            logger_name=logger_name,
            **kwargs
        )
        # XXX: Use the following with plot_ui.py
#        spectral_network_plot.draw(
#            sw_data=sw_data,
#            spectral_network=spectral_network,
#            logger_name=logger_name,
#            **kwargs
#        )
#        plot_legend = spectral_network_plot.get_legend()
        logger.info(plot_legend)

        if show_plot is True:
            spectral_network_plot.show()
        
        return spectral_network_plot

    else:
        # This case is assumed to correspond to having a parameter sequence
        plot_sequence = []
        for i, sw_data in enumerate(spectral_network_data.sw_data):
            spectral_networks = [spectral_network_data.spectral_networks[i]]
            spectral_network_plot_title = 'Spectral Network '+str(i)

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
            plot_sequence.append(spectral_network_plot)

            # Rotate the z-plane into the location defined by the curve.
            z_r = sw_data.z_plane_rotation
            sw_data.set_z_rotation(z_r)
            for sn in spectral_networks:
                sn.set_z_rotation(z_r)

            for spectral_network in spectral_networks:
                plot_legend = spectral_network_plot.draw(
                    spectral_network,
                    sw_data.branch_points,
                    punctures=(sw_data.regular_punctures +
                               sw_data.irregular_punctures),
                    irregular_singularities=sw_data.irregular_singularities,
                    g_data=sw_data.g_data,
                    branch_cut_rotation=sw_data.branch_cut_rotation,
                    logger_name=logger_name,
                    **kwargs
                )
                logger.info(plot_legend)

        return plot_sequence


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
