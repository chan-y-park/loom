import multiprocessing
import time
import pdb
import flask

from api import (
    set_logging,
)

# Flask configuration
DEBUG = True
SECRET_KEY = 'web_loom_key'

class WebLoomApplication(flask.Flask):
    def __init__(self, config_file, logging_level):
        super(WebLoomApplication, self).__init__('web_loom')
        self.app_id = str(int(time.time()))     # timestamp
        self.logging_queue = multiprocessing.Queue()
        set_logging(logging_level, queue=self.logging_queue,)

    def test(self):
        print flask.current_app.app_id

def index():
    flask.current_app.test()
    return 'index'

def progress():
    flask.current_app.test()
    return 'progress'

def get_application(config_file, logging_level):
    application = WebLoomApplication(config_file, logging_level)
    application.config.from_object(__name__)
    application.add_url_rule('/', 'index', index)
    application.add_url_rule('/progress', 'progress', progress)
    return application
