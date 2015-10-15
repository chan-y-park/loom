#!/usr/bin/env python

import os
import logging
#import sys
import pdb

from loom.web_ui import get_application

#argv = sys.argv[1:]
#argv += ['--gui-mode', 'True']

#run_with_sys_argv(argv)
config_file = os.path.join(os.curdir, 'default.ini')
application = get_application(config_file, logging_level=logging.INFO)

if __name__ == '__main__':
    host = '0.0.0.0'
    port = 8888
    print 'Listening on {}:{}.'.format(host, port)
    application.run(
        host=host,
        port=port,
        debug=True,
        use_reloader=False,
    )
