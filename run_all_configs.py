#!/usr/bin/env python

import os
import sys
import glob
import logging
import loom

logger_name = 'loom_test'

config_file_paths = glob.glob(
    os.path.join(
        loom.api.get_loom_dir(),
        'config/*.ini',
    )
)

loom.set_logging(
    logger_name=logger_name,
    logging_level=logging.INFO,
    logging_stream=sys.stdout,
    logging_file_name=os.path.join(
        loom.api.get_loom_dir(),
        'logs/loom_test.log',
    ),
)

logger = logging.getLogger(logger_name)

logger.info('Files to test:')
for file_path in config_file_paths:
    logger.info('{}'.format(file_path))

for file_path in config_file_paths:
    logger.info('Testing {}...'.format(file_path))
    c = loom.load_config(file_path, logger_name=logger_name)
    d = loom.generate(c, phase=1.0, logger_name=logger_name)
