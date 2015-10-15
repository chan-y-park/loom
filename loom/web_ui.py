import flask
import multiprocessing

from api import (
    set_logging,
)

# Flask configuration
DEBUG = True
SECRET_KEY = 'web_loom_key'

class WebLoomApplication(flask.Flask):
    def __init__(self, config_file, logging_level):
        super(WebLoomApplication, self).__init__('web_loom')
        self.logging_queue = multiprocessing.Queue()
        set_logging(logging_level, queue=self.logging_queue,)

def index():
    return 'Hello World from loom!'

def get_application(config_file, logging_level):
    application = WebLoomApplication(config_file, logging_level)
    application.config.from_object(__name__)
    application.add_url_rule('/', 'index', index)
    return application
