#!/usr/bin/env python

import scipy
import os
import logging
import sys
import pdb
import flask
import time

LOOM_DIR = os.path.dirname(os.path.realpath(__file__))
#os.chdir(LOOM_DIR)
#sys.path.insert(0, LOOM_DIR)
sys.stdout = sys.stderr
print 'loom WSGI working directory: {}'.format(os.getcwd())

from loom.web_ui import get_application

#argv = sys.argv[1:]
#argv += ['--gui-mode', 'True']

#run_with_sys_argv(argv)
DEFAULT_CONFIG_FILE = os.path.join(LOOM_DIR, 'default.ini')
application = get_application(DEFAULT_CONFIG_FILE, logging_level=logging.INFO)

if __name__ == '__main__':
    host = '0.0.0.0'
    port = 8888
    try:
        application.run(
            host=host,
            port=port,
            debug=True,
            use_reloader=True,
            threaded=True,
        )
    except (KeyboardInterrupt, SystemExit):
        raise
