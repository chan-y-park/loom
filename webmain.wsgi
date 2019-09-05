#!/usr/bin/env python3

import os
import signal
import logging
import sys
import getopt

import scipy
import flask
import time

# import pdb

argv = sys.argv[1:]
print('__name__={}'.format(__name__))
print('loom WSGI working directory: {}'.format(os.getcwd()))
if '_mod_wsgi_' in __name__:
    f = open(os.devnull, 'w')
    sys.stdout = f
    print('stdout redirecting failed if this message is shown.')

from loom.api import get_loom_dir
from loom.web_ui import get_application

default_config_file = os.path.join(get_loom_dir(), 'config/default.ini')
application = get_application(
    default_config_file,
    logging_level=logging.INFO
)

if __name__ == '__main__':
    host = '0.0.0.0'
    port = 8888
    optlist, args = getopt.getopt(argv, 'p:')
    for opt, arg in optlist:
        if opt == '-p':
            port = int(arg)

    pid_file = os.path.join(get_loom_dir(), 'logs/webmain_pid')
#    try:
#        with open(pid_file, 'r') as fp:
#            old_pid = int(fp.read())
#            os.kill(old_pid, signal.SIGKILL)
#    except (IOError, OSError):
#        pass
    with open(pid_file, 'w') as fp:
        pid = os.getpid()
        fp.write(str(pid))
    try:
        application.run(
            host=host,
            port=port,
            debug=True,
            use_reloader=False,
            threaded=True,
        )
    except (KeyboardInterrupt, SystemExit):
        raise
