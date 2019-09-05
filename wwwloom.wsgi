#!/usr/bin/env python3

import os
import sys
import logging

# Adapt PATH and activate virtualenv
BASEDIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, BASEDIR)
activate_this = BASEDIR + '/.venv/bin/activate_this.py'
with open(activate_this) as file_:
    exec(file_.read(), dict(__file__=activate_this))

if '_mod_wsgi_' in __name__:
    f = open(os.devnull, 'w')
    sys.stdout = f
    print('stdout redirecting failed if this message is shown.')

sys.path.insert(0, BASEDIR + '/loom')

from loom.api import get_loom_dir
from loom.web_ui import get_application

default_config_file = os.path.join(get_loom_dir(), 'config/default.ini')
application = get_application(
    default_config_file,
    logging_level=logging.INFO
)
