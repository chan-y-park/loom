import ConfigParser
import os
import sys
import logging

from math import pi

CONFIG_FILE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               'config_file')


class ConfigData:
    def __init__(self, opts):
        #self.phase = opts['phase']
        #self.single_network = opts['single-network']
        if(opts['phase'] is not None):
            self.phase_range = None

        self.set_logger(opts['logging-level'])
        self.read_config_from_file(opts['config-file'])
        
    def set_logger(self, level):
        if len(level) == 0:
            level = 'warning'

        if level == 'debug':
            logging_level = logging.DEBUG
            logging_format = '%(module)s@%(lineno)d: %(funcName)s: %(message)s'
        elif level == 'info':
            logging_level = logging.INFO
            logging_format = '%(process)d: %(message)s'
        elif level == 'warning':
            logging_level = logging.WARNING
            logging_format = '%(message)s'

        logging.basicConfig(level=logging_level, format=logging_format, 
                            stream=sys.stdout)
        #self.logger = logging.getLogger('loom_log')

    def read_config_from_file(self, config_file):
        if config_file is None:
            config_file = 'default.ini'
        config_parser = ConfigParser.SafeConfigParser()
        config_file_full_path = os.path.join(CONFIG_FILE_DIR, config_file)
        logging.debug('config file full path: %s', config_file_full_path)
        config_parser.read(config_file_full_path)
        # directories
        self.root_dir = config_parser.get('directory', 'root_dir')
        self.config_dir = config_parser.get('directory', 'config_dir')
        self.data_dir = config_parser.get('directory', 'data_dir')

        self.sw_curve_eq_string = config_parser.get('analytic expressions',
                                                    'sw_curve')
        self.sw_diff_v_string = config_parser.get('analytic expressions',
                                                  'sw_diff_v')
        # construct a dict of curve parameters
        self.sw_parameters = {}
        for param in config_parser.options('Seiberg-Witten parameters'):
            self.sw_parameters[param] = config_parser.getfloat(
                'Seiberg-Witten parameters', param
            )
        
        # numerical parameters
        self.phase_range = eval(config_parser.get(
            'numerical parameters',
            'phase_range'
        ))
        self.z_range_limits = eval(config_parser.get(
            'numerical parameters',
            'z_range_limits'
        ))
        self.num_of_steps = config_parser.getint(
            'numerical parameters',
            'num_of_steps'
        )
        self.num_of_iterations = config_parser.getint(
            'numerical parameters',
            'num_of_iterations'
        )
        self.size_of_small_step = eval(config_parser.get(
            'numerical parameters',
            'size_of_small_step'
        ))
        self.size_of_large_step = eval(config_parser.get(
            'numerical parameters',
            'size_of_large_step'
        ))
        self.size_of_neighborhood = eval(config_parser.get(
            'numerical parameters',
            'size_of_neighborhood'
        ))
        self.size_of_bin = eval(config_parser.get(
            'numerical parameters',
            'size_of_bin'
        ))
        self.accuracy = eval(config_parser.get(
            'numerical parameters',
            'accuracy'
        ))
        self.n_processes = config_parser.getint(
            'numerical parameters',
            'n_processes'
        )
