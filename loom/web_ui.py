import multiprocessing
import time
import pdb
import flask
import sys
import logging
import uuid

from StringIO import StringIO
from Queue import Empty as QueueEmpty
#sys.stdout = sys.stderr

from api import (
    get_logging_handler,
    set_logging,
    load_config,
    generate_spectral_network,
)

# Flask configuration
DEBUG = True
SECRET_KEY = 'web_loom_key'
PARENT_LOGGER_NAME = 'web_loom_logger'

class LoomDB(object):
    def __init__(self):
        self.logging_queues = {}
        self.result_queues = {}
        self.loom_processes = {}

    def start_loom_process(
        self, process_uuid, logging_level, loom_config, phase=None,
    ):
        logging_queue = multiprocessing.Queue()
        self.logging_queues[process_uuid] = logging_queue
        logger_name = get_logger_name(process_uuid)
        set_logging(
            logger_name=logger_name,
            logging_level=logging_level,
            logging_queue=logging_queue,
        )

        result_queue = multiprocessing.Queue()
        self.result_queues[process_uuid] = result_queue

        loom_process = multiprocessing.Process(
            target = generate_spectral_network,
            args=(
                loom_config,
            ),
            kwargs=dict(
                phase=phase,
                result_queue=result_queue,
                logger_name=logger_name,
            ),
        )
        self.loom_processes[process_uuid] = loom_process
        loom_process.start()

        return None

    def get_log_message(
        self, process_uuid, logging_stream, logging_stream_handler
    ):
        record = self.logging_queues[process_uuid].get_nowait()
        if record is not None:
            logging_stream_handler.handle(record)
            logs = logging_stream.getvalue()
            logging_stream.truncate(0)
            return logs

    def yield_log_message(self, process_uuid, logging_level,):
        logging_stream = StringIO()
        logging_stream_handler = get_logging_handler(
            logging_level,
            logging.StreamHandler,
            logging_stream,
        )

        while self.result_queues[process_uuid].empty() is True:
            try:
                logs = self.get_log_message(process_uuid, logging_stream,
                                            logging_stream_handler,)
                yield '<br>{}</br>\n'.format(logs)
            except QueueEmpty:
                pass 
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                import sys, traceback
                print >> sys.stderr, 'logging_listener_process:'
                traceback.print_exc(file=sys.stderr)

        # Get the remaining logs, if any.
        while True:
            try:
                logs = self.get_log_message(process_uuid, logging_stream,
                                            logging_stream_handler,)
                yield '<br>{}</br>\n'.format(logs)
            except QueueEmpty:
                break 
            except (KeyboardInterrupt, SystemExit):
                raise
        # Draw 'Plot' button.
        yield '<input type=submit class="btn btn-default" value=Plot>'
                
    def finish_loom_process(self, process_uuid):
        logger_name = get_logger_name(process_uuid)
        logger = logging.getLogger(logger_name)
        result_queue = self.result_queues[process_uuid]
        loom_process = self.loom_processes[process_uuid]
        if result_queue.empty() is True:
            if loom_process.is_alive():
                logger.warning('Job {} still running.'.format(name))
                return None
            else:
                logger.warning(
                    'Generating spectral networks failed: '
                    'pid = {}, exitcode = {}.'
                    .format(generate_process.pid,
                            generate_process.exitcode,)
                )

        spectral_network_data = result_queue.get()
        loom_process.join()
        sw_data = spectral_network_data.sw_data
        spectral_networks = spectral_network_data.spectral_networks
        logger.info('Finished generating spectral network data.')
        # Remove the result queue
        # Remove the logging queue.
        # Remove the logging queue handler.
        # Call the plotting function.
        return None

class WebLoomApplication(flask.Flask):
    def __init__(self, config_file, logging_level):
        super(WebLoomApplication, self).__init__('web_loom')
        set_logging(
            logger_name=PARENT_LOGGER_NAME,
            logging_level=logging_level,
            logging_file_name='logs/web_loom.txt',
        )
        self.loom_db = LoomDB()

def index():
    return flask.render_template('index.html')

def config():
    loom_config = None
    if flask.request.method == 'GET':
        return flask.render_template('config.html')
    elif flask.request.method == 'POST':
        phase = eval(flask.request.form['phase'])
        process_uuid = str(uuid.uuid4())
        logger_name = get_logger_name(process_uuid)
        loom_config = get_loom_config(flask.request.form, logger_name) 
        #app = flask.current_app._get_current_object()
        app = flask.current_app
        app.loom_db.start_loom_process(
            process_uuid, logging.INFO, loom_config, phase,
        )
        return flask.redirect(
            flask.url_for('progress', process_uuid=process_uuid) 
        )

def progress_stream_template(template_name, **context):
    # http://flask.pocoo.org/docs/patterns/streaming/#streaming-from-templates
    app = flask.current_app
    app.update_template_context(context)
    t = app.jinja_env.get_template(template_name)
    rv = t.stream(context)
    return rv

def progress(process_uuid):
    if flask.request.method == 'GET':
        app = flask.current_app
        return flask.Response(
            progress_stream_template(
                'progress.html',
                process_uuid=process_uuid,
                progress=app.loom_db.yield_log_message(
                    process_uuid, logging.INFO,
                )
            )
        )
    elif flask.request.method == 'POST':
        return flask.redirect(
            flask.url_for(
                'plot',
                #process_uuid=flask.request.form['process_uuid'],
                process_uuid=process_uuid,
            )
        )

def plot(process_uuid):
    # Finish loom_process
    app = flask.current_app
    app.loom_db.finish_loom_process(process_uuid)
    # Make a Bokeh plot
    return flask.render_template(
        'plot.html',
        # Bokeh results
    )

def get_application(config_file, logging_level):
    application = WebLoomApplication(config_file, logging_level)
    application.config.from_object(__name__)
    application.add_url_rule(
        '/', 'index', index, methods=['GET'],
    )
    application.add_url_rule(
        '/config', 'config', config, methods=['GET', 'POST'],
    )
    application.add_url_rule(
        '/progress/<process_uuid>', 'progress', progress,
        methods=['GET', 'POST'],
    )
    application.add_url_rule(
        '/plot/<process_uuid>', 'plot', plot,
        methods=['GET'],
    )
    return application

def get_logger_name(uuid):
    return PARENT_LOGGER_NAME + '.' + uuid

def get_loom_config(request_dict, logger_name):
    loom_config = load_config('default.ini', logger_name=logger_name)
    return loom_config
