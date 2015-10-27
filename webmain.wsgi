#!/usr/bin/env python

import scipy
import os
import logging
import sys
import pdb
import flask
import time

sys.stdout = sys.stderr
print 'loom WSGI working directory: {}'.format(os.getcwd())

from loom.api import get_loom_dir
from loom.web_ui import get_application

#argv = sys.argv[1:]

default_config_file = os.path.join(get_loom_dir(), 'config/default.ini')
application = get_application(default_config_file, logging_level=logging.INFO)

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
