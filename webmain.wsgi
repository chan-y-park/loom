#!/usr/bin/env python

import os
import signal
import logging
import sys

import scipy
import flask
import time

import pdb

sys.stdout = sys.stderr
print 'loom WSGI working directory: {}'.format(os.getcwd())

from loom.api import get_loom_dir
from loom.web_ui import get_application

#argv = sys.argv[1:]

default_config_file = os.path.join(get_loom_dir(), 'config/default.ini')
application = get_application(default_config_file, logging_level=logging.INFO)

if __name__ == '__main__':
    pid_file = os.path.join(get_loom_dir(), 'logs/webmain_pid')
    try:
        with open(pid_file, 'r') as fp:
            old_pid = int(fp.read())
            os.kill(old_pid, signal.SIGKILL)
    except (IOError, OSError):
        pass
    with open(pid_file, 'w') as fp:
        pid = os.getpid()
        fp.write(str(pid))
    host = '0.0.0.0'
    port = 8888
    try:
        application.run(
            host=host,
            port=port,
            debug=True,
            #use_reloader=True,
            use_reloader=False,
            threaded=True,
        )
    except (KeyboardInterrupt, SystemExit):
        raise
