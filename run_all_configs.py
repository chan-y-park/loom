#!/usr/bin/env python
"""
Run all the configs in the /config folder (but not recursively).
Generate a spectral network at a single phase
and produce a Bokeh plot for each config in /test_results folder
"""

import os
import sys
import glob
import logging
import loom
import pdb

from bokeh.io import save
from loom.bokeh_plot import get_spectral_network_bokeh_plot

passed_configs = []
# passed_configs = [
#     "/home/chan/loom/config/coset_A_3.ini",
#     "/home/chan/loom/config/D_4_AD_from_A_1.ini",
#     "/home/chan/loom/config/AD_SO_6.ini",
#     "/home/chan/loom/config/pure_SO_4.ini",
#     "/home/chan/loom/config/pure_SO_8.ini",
#     "/home/chan/loom/config/D_3_AD_from_SU_4.ini",
#     "/home/chan/loom/config/coset_A_3_2nd_rep.ini",
#     "/home/chan/loom/config/AD_from_SU_4_N_f_2.ini",
#     "/home/chan/loom/config/coset_A_4_2nd_rep.ini",
#     "/home/chan/loom/config/default.ini",
#     "/home/chan/loom/config/SU_4_N_f_1.ini",
#     "/home/chan/loom/config/AD_from_pure_SU_4.ini",
#     "/home/chan/loom/config/D_4_AD_from_SO_8.ini",
#     "/home/chan/loom/config/coset_D_4.ini",
#     "/home/chan/loom/config/coset_D_3.ini",
#     "/home/chan/loom/config/pure_SO_6.ini",
#     "/home/chan/loom/config/AD_SO_4.ini",
#     "/home/chan/loom/config/seeding_type_II.ini",
#     "/home/chan/loom/config/triangle_4.ini",
#     "/home/chan/loom/config/D_4_AD_from_A_2.ini",
#     "/home/chan/loom/config/pure_SU_2.ini",
#     "/home/chan/loom/config/seeding_type_I.ini",
#     "/home/chan/loom/config/seeding_type_III_A.ini",
#     "/home/chan/loom/config/ex_issue_2.ini",
# ]

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
    if file_path in passed_configs:
        logger.info('{} already passed the test, skip it.'.format(file_path))
        continue
    logger.info('Testing {}...'.format(file_path))
    c = loom.load_config(file_path, logger_name=logger_name)
    d = loom.generate(c, phase=1.0, logger_name=logger_name)
    p = get_spectral_network_bokeh_plot(d, plot_range=c['plot_range'],
                                        notebook=True)
    file_dir, file_name = os.path.split(file_path)
    name, ext = os.path.splitext(file_name)
    save(obj=p, filename='test_results/{}.html'.format(name), title=name,)
