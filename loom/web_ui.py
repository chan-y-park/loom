import os
import multiprocessing
import time
import pdb
import flask
import sys
import logging
import uuid
import numpy
import bokeh

from StringIO import StringIO
from Queue import Empty as QueueEmpty
#sys.stdout = sys.stderr


from api import (
    get_loom_dir,
    get_logging_handler,
    set_logging,
    load_config,
    generate_spectral_network,
)

from bokeh_plot import get_spectral_network_bokeh_plot

# Flask configuration
DEBUG = True
SECRET_KEY = 'web_loom_key'
PARENT_LOGGER_NAME = 'loom'
WEB_APP_NAME = 'web_loom'
LOGGING_FILE_PATH = os.path.join(
    get_loom_dir(),
    'logs/web_loom.log',
)


# TODO: kill an orphaned process gracefully.
class LoomDB(object):
    """
    The internal DB to manage loom processes.
    """
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
                logging_queue=logging_queue,
                logger_name=logger_name,
            ),
        )
        self.loom_processes[process_uuid] = loom_process
        loom_process.start()

        return None

    def get_log_message(
        self, process_uuid, logging_stream, logging_stream_handler
    ):
        record = self.logging_queues[process_uuid].get(True, 3)
        if record is not None:
            logging_stream_handler.handle(record)
            logs = logging_stream.getvalue()
            logging_stream.truncate(0)
            return logs
        else:
            raise QueueEmpty

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
                yield 'data: {}\n\n'.format(logs)
            except QueueEmpty:
                pass 
            except (KeyboardInterrupt, SystemExit):
                raise
            except:
                import traceback
                print >> sys.stderr, 'logging_listener_process:'
                traceback.print_exc(file=sys.stderr)

        # Get the remaining logs, if any.
        while True:
            try:
                logs = self.get_log_message(process_uuid, logging_stream,
                                            logging_stream_handler,)
                yield 'data: {}\n\n'.format(logs)
            except QueueEmpty:
                break 
            except (KeyboardInterrupt, SystemExit):
                raise
        yield 'event: finish\ndata: \n\n'
                
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
                    .format(loom_process.pid,
                            loom_process.exitcode,)
                )

        spectral_network_data = result_queue.get()
        sw_data = spectral_network_data.sw_data
        spectral_networks = spectral_network_data.spectral_networks
        logger.info('Finished generating spectral network data.')

        loom_process.join()
        del self.loom_processes[process_uuid]

        # Remove the result queue
        del self.result_queues[process_uuid]
        # Remove the logging queue.
        del self.logging_queues[process_uuid]
        # Remove the logging queue handler.
        logger.handlers = []
        # Remove the logger.
        del logging.Logger.manager.loggerDict[logger_name]
        # Call the plotting function.
        bokeh_layout = get_spectral_network_bokeh_plot(spectral_network_data)
        script, div = bokeh.embed.components(bokeh_layout)
        legend = get_plot_legend(sw_data)
        return (script, div, legend,) 


class WebLoomApplication(flask.Flask):
    """
    A wrapper of Flask application containing an instance of LoomDB.
    """
    def __init__(self, config_file, logging_level):
        super(WebLoomApplication, self).__init__(WEB_APP_NAME)
        set_logging(
            logger_name=WEB_APP_NAME,
            logging_level=logging_level,
            logging_file_name=LOGGING_FILE_PATH,
        )
        self.loom_db = LoomDB()

###
# View functions
###
def index():
    return flask.render_template('index.html')

def config():
    loom_config = None
    if flask.request.method == 'GET':
        return flask.render_template(
            'config.html',
            event_source_url = None,
        )
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
        return flask.render_template(
            'config.html',
            event_source_url = flask.url_for(
                'logging_stream', process_uuid=process_uuid,
            ),
            text_area_content = "Start loom, uuid = {}".format(process_uuid),
            plot_url = flask.url_for(
                'plot', process_uuid=process_uuid,
            ),
        )

def logging_stream(process_uuid):
    if flask.request.headers.get('accept') == 'text/event-stream':
        app = flask.current_app
        return flask.Response(
            app.loom_db.yield_log_message(process_uuid, logging.INFO),
            mimetype='text/event-stream',
        )
        
def plot(process_uuid):
    # Finish loom_process
    app = flask.current_app
    script, div, legend = app.loom_db.finish_loom_process(process_uuid)
    # Make a Bokeh plot
    return flask.render_template(
        'plot.html',
        plot_script=script,
        plot_div=div,
        plot_legend=legend,
    )

def admin():
    # TODO: password-protect this page.
    app = flask.current_app
    loom_db = app.loom_db
    pdb.set_trace()

###
# Entry point
###

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
        '/plot/<process_uuid>', 'plot', plot,
        methods=['POST'],
    )
    application.add_url_rule(
        '/logging_stream/<process_uuid>', 'logging_stream', logging_stream,
        methods=['GET'],
    )
    application.add_url_rule(
        '/admin', 'admin', admin,
        methods=['GET'],
    )
    return application

###
# Misc. web UIs
###

def get_logger_name(uuid):
    return WEB_APP_NAME + '.' + uuid

def get_loom_config(request_dict, logger_name):
    loom_config = load_config('default.ini', logger_name=logger_name)
    return loom_config

def get_plot_legend(sw_data):
    legend = ''
    g_data = sw_data.g_data
    roots = g_data.roots
    weights = g_data.weights
    #root_dictionary = make_root_dictionary(g_data)
    #weight_dictionary = make_weight_dictionary(g_data)
    #root_labels = root_dictionary.keys()
    #roots = root_dictionary.values()
    weight_pairs=[
        [str('(mu_'+str(p[0])+', mu_'+str(p[1])+')') 
         for p in g_data.ordered_weight_pairs(rt)]
        for rt in roots
    ]
    #weight_labels = weight_dictionary.keys()
    #weights = weight_dictionary.values()

    legend += ('\t--- The Root System ---\n')
    for i in range(len(roots)):
        legend += (
            'alpha_' + str(i) + ' : {}\n'.format(list(roots[i])) +
            'ordered weight pairs : {}\n'.format(weight_pairs[i])
        )

    legend += ('\t--- The Weight System ---\n')
    for i in range(len(weights)):
        legend += (
            'nu_' + str(i) + ' : {}\n'.format(list(weights[i]))
        )

    legend += ('\t--- The Branch Points ---\n')
    for bp in sw_data.branch_points:
        root_labels = []
        for pr in bp.positive_roots:
            for i, r in enumerate(roots):
                if numpy.array_equal(pr, r):
                    root_labels.append('alpha_{}'.format(i)) 
        legend += (
            bp.label + 
            '\tposition : {}\n'.format(bp.z) +
            '\t\troot type : {}\n'.format(root_labels) +
            '\t\tmonodromy matrix : \n{}\n'.format(bp.monodromy)
        )

    legend += ('\t--- The Irregular Singularities ---\n')
    for irs in sw_data.irregular_singularities:
        legend += (
            irs.label + 
            '\tposition : {}\n'.format(irs.z) + 
            '\tmonodomry matrix : \n{}\n'.format(irs.monodromy)
        )

    return legend


