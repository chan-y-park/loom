#!/usr/bin/env python

import os
import logging
#import sys

from loom.api import get_loom_dir
from loom.gui import open_gui

#argv = sys.argv[1:]

config_file_path = os.path.join(get_loom_dir(), 'config/default.ini')
open_gui(config_file_path=config_file_path, logging_level=logging.INFO)

