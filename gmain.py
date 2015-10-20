#!/usr/bin/env python

import os
import logging
#import sys
#from loom.__main__ import run_with_sys_argv

from loom.gui import open_gui

#argv = sys.argv[1:]
#argv += ['--gui-mode', 'True']

#run_with_sys_argv(argv)
config_file_path = os.path.join(os.curdir, 'default.ini')
open_gui(config_file_path=config_file_path, logging_level=logging.INFO)

