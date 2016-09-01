import ConfigParser
import logging
# import pdb


class LoomConfig:
    """
    A container class of the configuration data.
    Saves the configuration data as a Python dictionary.
    """
    def __init__(self, file_path=None, logger_name='loom'):
        logger = logging.getLogger(logger_name)
        self.logger_name = logger_name
        # The data dict has option values as Python objects,
        # some of them being strings.
        self.data = {}
        self.parser = ConfigParser.SafeConfigParser()
#        # attribute options has the following structure:
#        # options['section']['option'] = 'label'
#        self.options = {
#            'Seiberg-Witten data': {
#                'description': 'Description',
#                'root_system': 'Root system',
#                'representation': 'Representation',
#                'casimir_differentials': 'Casimir differentials',
#                'differential_parameters': 'Parameters of differentials',
#                'parameter_sequence': 'A sequence of values for a parameter',
#                'regular_punctures': 'Regular punctures',
#                'irregular_punctures': 'Irregular punctures',
#                'ramification_points': 'Ramification points',
#                'branch_points': 'Branch points',
#                'branch_points_sequence': 'Branch points sequence',
#                'mt_params': 'Mobius transformation',
#                'ramification_point_finding_method':
#                    'Ramification point finding method',
##                'integration_method': 'Integration Method',
#            },
#            'numerical parameters': {
#                'accuracy': 'Accuracy',
#                'plot_range': 'Plot range',
#                'num_of_steps': 'Number of steps',
#                'num_of_iterations': 'Number of iterations',
#                'size_of_small_step': 'Size of a small step',
#                'size_of_large_step': 'Size of a large step',
#                'size_of_bp_neighborhood':
#                    'Size of a branch point neighborhood',
#                'size_of_puncture_cutoff': 'Size of a puncture cutoff',
#                'mass_limit': 'Mass limit',
#                'phase': 'Phase (single value or range)',
#            },
#            'settings': {
#                'trivialize': 'Trivialize the cover',
#                'use_scipy_ode': 'Use SciPy ODE solver',
#            },
#        }
#
#        # {deprecated option: option,}
#        self.deprecated_options = {
#            'phase_range': 'phase',
#            'punctures': 'irregular_punctures',
#            'size_of_neighborhood': 'size_of_bp_neighborhood',
#            'size_of_bin': None,
#            'size_of_ramification_pt_cutoff': None,
#            'n_processes': None,
#        }

        # attribute options has the following structure:
        # options['section']['option'] = (default value) 
        self.options = {
            'Seiberg-Witten data': {
                'description': None,
                'root_system': None,
                'representation': None,
                'casimir_differentials': None,
                'differential_parameters': None,
                'parameter_sequence': None,
                'regular_punctures': None,
                'irregular_punctures': None,
                'ramification_points': None,
                'branch_points': None,
                'branch_points_sequence': None,
                'mt_params': None,
                'ramification_point_finding_method': None,
#                'integration_method': None,
            },
            'numerical parameters': {
                'accuracy': None,
                'plot_range': None,
                'num_of_steps': None,
                'num_of_iterations': None,
                'size_of_small_step': None,
                'size_of_large_step': None,
                'size_of_bp_neighborhood': None,
                'size_of_puncture_cutoff': None,
                'mass_limit': None,
                'phase': None,
            },
            'settings': {
                'trivialize': True,
                'use_scipy_ode': True,
            },
        }

        # {deprecated option: option,}
        self.deprecated_options = {
            'phase_range': 'phase',
            'punctures': 'irregular_punctures',
            'size_of_neighborhood': 'size_of_bp_neighborhood',
            'size_of_bin': None,
            'size_of_ramification_pt_cutoff': None,
            'n_processes': None,
            'integration_method': 'use_scipy_ode',
        }

        if file_path is not None:
            logger.info('Loading configuration from {}...'.format(file_path))
            with open(file_path, 'r') as fp:
                self.read(fp)
            logger.info('Finished loading configuration from {}.'
                        .format(file_path))

        # Default settings.
        if self.parser.has_section('settings') is False:
            self.parser.add_section('settings')
        for option, default in self.options['settings'].iteritems():
            try:
                if self[option] is not None:
                    continue
            except KeyError:
                pass
            self[option] = default


    def __setitem__(self, option, value):
        try:
            # Update the data dict.
            self.data[option] = value
            # Update the parser
            for section in self.options:
                if option in self.options[section]:
                    self.parser.set(section, option, str(value))
        except KeyError:
            # Not an available the option, raise an error.
            raise KeyError('Unknown option \'{}\'.'.format(option))

    def __getitem__(self, option):
        try:
            return self.data[option]
        except KeyError:
            # Not an available the option, raise an error.
            raise KeyError('Unknown option \'{}\'.'.format(option))

    def __delitem__(self, option):
        try:
            # Remove the option from the data dict.
            del self.data[option]
            # Remove the option from the parser.
            for section in self.options:
                if option in self.options[section]:
                    self.parser.remove_option(section, option)
        except KeyError:
            raise KeyError('Unknown option \'{}\'.'.format(option))

    def get_label(self, option):
        for section in self.options.keys():
            if option in self.options[section]:
                return self.options[section][option]

        raise ValueError('Unknown option \'{}\'.'.format(option))

    def keys(self):
        return self.data.keys()

    def iteritems(self):
        return self.data.iteritems()

    def read(self, config_file):
        """
        Read an .ini file and load the configuration data.
        """
        logger = logging.getLogger(self.logger_name)

        self.parser.readfp(config_file)

        for section in self.parser.sections():
            try:
                not_configured_options = self.options[section].keys()
            except KeyError:
                self.parser.remove_section(section)
                raise ValueError('Unknown section \'{}\'. in the config file.'
                                 .format(section))

            for option in self.parser.options(section):
                parser_value = self.parser.get(section, option)

                # Check deprecated options
                if option in self.deprecated_options:
                    old_option = option
                    option = self.deprecated_options[old_option]
                    logger.warning(
                        'Option \'{}\' is deprecated.'.format(old_option)
                    )
                    if option is None:
                        # There is no newer version of this option,
                        # discard the deprecated option.
                        continue
                    elif old_option is 'integration_method':
                        # XXX: temporary
                        section = 'settings'
                        # option = 'use_scipy_ode'
                        if parser_value == 'ode_int':
                            parser_value = 'True'
                        elif parser_value == 'manual':
                            parser_value = 'False'
                        else:
                            raise RuntimeError(
                                'Unknown option: {} = {}'
                                .format(old_option, parser_value)
                            )

                    logger.warning('Use \'{}\' instead.'.format(option))
                    self.parser.remove_option(section, old_option)
                    self.parser.set(section, option, parser_value)

                if (parser_value == 'None'):
                    value = eval(parser_value)
                elif (section == 'Seiberg-Witten data'):
                    value = parser_value
                elif (section == 'numerical parameters'
                      or section == 'settings'):
                    value = eval(parser_value)
                else:
                    raise ValueError(
                        'Option \'{}\' in an unknown section \'{}\'.'
                        .format(option, section)
                    )

                try:
                    not_configured_options.remove(option)
                except ValueError:
                    raise ValueError('Unknown option \'{}\'.'.format(option))

                self.data[option] = value

            for option in not_configured_options:
                self.parser.set(section, option, 'None')
                self[option] = None

#            # default integration is ode_int
#            if self.data['integration_method'] is None:
#                self.data['integration_method'] = 'ode_int'

    def save(self, file_path=None, logger_name=None,):
        logger = logging.getLogger(logger_name)
        if file_path is None:
            return None
        logger.info('Saving configuration to {}.'.format(file_path))
        with open(file_path, 'w') as fp:
            self.parser.write(fp)
        logger.info('Finished saving configuration to {}.'.format(file_path))
