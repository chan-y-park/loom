import ConfigParser
import os
import sys
import logging

CONFIG_FILE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                              'config_file')

class ConfigData:
    def __init__(self, opts):
        self.set_logger(opts['logging-level'])
        self.phase = opts['phase']
        self.single_network = opts['single-network']
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
        if len(config_file) == 0:
            config_file = 'default.ini'
        config_parser = ConfigParser.SafeConfigParser()
        config_file_full_path = os.path.join(CONFIG_FILE_DIR, config_file)
        logging.debug('config file full path: %s', config_file_full_path)
        config_parser.read(config_file_full_path)

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

        self.accuracy = config_parser.getfloat('numerical parameters',
                                          'accuracy')
