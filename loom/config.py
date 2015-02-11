import ConfigParser
import os
import sys
import logging
import pdb

from math import pi


class ConfigData:
    def __init__(self, opts):
        if (opts['phase'] is not None):
            self.phase_range = None
        self.SL2C_params = None

    def read_config_from_file(self, config_file):
        logging.info('config file: %s', config_file)
        config_parser = ConfigParser.SafeConfigParser()
        config_parser.read(config_file)

        # directories
        self.root_dir = config_parser.get('directory', 'root_dir')
        self.config_dir = config_parser.get('directory', 'config_dir')
        self.data_dir = config_parser.get('directory', 'data_dir')

        # analytic expressions
        self.sw_curve_eq_string = config_parser.get('analytic expressions',
                                                    'sw_curve')
        self.sw_diff_v_string = config_parser.get('analytic expressions',
                                                  'sw_diff_v')
        # construct a dict of curve parameters
        self.sw_parameters = {}
        for param in config_parser.options('Seiberg-Witten parameters'):
            self.sw_parameters[param] = eval(config_parser.get(
                'Seiberg-Witten parameters', param
            ))
        
        # numerical parameters
        self.SL2C_params = eval(config_parser.get(
            'numerical parameters',
            'SL2C_params'
        ))
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

        self.config_parser = config_parser
