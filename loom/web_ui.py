import os
import subprocess
import signal
import multiprocessing
import threading
import time
# import pdb
import flask
import sys
import logging
import uuid
import json
import zipfile

from cStringIO import StringIO
from io import BytesIO
from Queue import Empty as QueueEmpty

from api import (
    get_loom_dir,
    get_logging_handler,
    set_logging,
    load_config,
    load_spectral_network,
    save_spectral_network,
    generate_spectral_network,
    get_current_branch_version,
)
from config import LoomConfig
from bokeh_plot import get_spectral_network_bokeh_plot
from plotting import get_legend
from misc import get_data_size_of
from cmath import pi

# Flask configuration
DEBUG = True
SECRET_KEY = 'web_loom_key'
PARENT_LOGGER_NAME = 'loom'
WEB_APP_NAME = 'web_loom'
STAT_LOGGER_NAME = 'stat_of_' + WEB_APP_NAME
DEFAULT_NUM_PROCESSES = 4
DB_CLEANUP_CYCLE_SECS = 60
LOOM_PROCESS_JOIN_TIMEOUT_SECS = 3

# Array of config options. 
# Entries that will be placed in the same row
# are in the same row of this array.
config_options = [
    ['description'],
    ['casimir_differentials'],
    ['differential_parameters'],
    ['regular_punctures'],
    ['irregular_punctures'],
    ['plot_range'],
    ['num_of_steps'],
    ['num_of_iterations'],
    ['mass_limit'],
    ['phase'],
]

advanced_config_options = [
    ['mt_params'], 
    ['ramification_point_finding_method'],
    ['size_of_small_step'],
    ['size_of_large_step'],
    ['size_of_bp_neighborhood'],
    ['size_of_puncture_cutoff'],
    ['accuracy'],
]

# TODO: kill an orphaned process gracefully.
class LoomDB(object):
    """
    The internal DB to manage loom processes.
    """
    def __init__(self):
        self.logging_queues = {}
        self.result_queues = {}
        self.loom_processes = {}
        self.is_alive = {}

        signal.signal(signal.SIGINT, self.loom_db_stop_signal_handler)
        signal.signal(signal.SIGTERM, self.loom_db_stop_signal_handler)

        self.db_manager = threading.Thread(
            target=self.db_manager,
        )
        self.db_manager_stop = threading.Event()
        self.db_manager_stop.clear()
        self.db_manager.start()

    def loom_db_stop_signal_handler(self, signum, frame):
        logger_name = get_logger_name()
        logger = logging.getLogger(logger_name)

        if signum == signal.SIGINT:
            msg = 'LoomDB caught SIGINT; raises KeyboardInterrrupt.'
            e = KeyboardInterrupt
        elif signum == signal.SIGTERM:
            msg = 'LoomDB caught SIGTERM; raises SystemExit.'
            e = SystemExit

        logger.warning(msg)
        self.db_manager_stop.set()
        raise e

    def db_manager(self):
        """
        A thread that will manage the DB
        and clean up data of previous clients.
        """
        logger_name = get_logger_name()
        logger = logging.getLogger(logger_name)

        # Wait for DB_CLEANUP_CYCLE_SECS and do clean-ups.
        # When self.db_manager_stop event is set,
        # break the while-loop and finish all the processes.
        while not self.db_manager_stop.wait(DB_CLEANUP_CYCLE_SECS):
            to_delete = []
            try:
                for process_uuid, alive in self.is_alive.iteritems():
                    if alive is True:
                        # Reset the heartbeat counter so that
                        # it can be cleaned up later.
                        self.is_alive[process_uuid] = False
                    else:
                        to_delete.append(process_uuid)

                        try:
                            # Flush the result queue.
                            result_queue = self.result_queues[process_uuid]
                            while result_queue.empty() is False:
                                result_queue.get_nowait()
                            # Remove the result queue
                            del self.result_queues[process_uuid]
                        except KeyError:
                            logger.warning(
                                "Removing result queue {} failed: "
                                "no such a queue exists."
                                .format(process_uuid)
                            ) 
                            pass

                for process_uuid in to_delete:
                    try:
                        del self.is_alive[process_uuid]
                    except KeyError:
                        logger.warning(
                            "Deleting is_alive[{}] failed."
                            .format(process_uuid)
                        ) 
                        pass

            except (KeyboardInterrupt, SystemExit) as e:
                msg = ('LoomDB manager thread caught {}; '
                       'finishes the manager.'.format(type(e)))
                logger.warning(msg)
                raise e(msg)
        
        # Received a stop event; finish all the processes.
        process_uuids = self.loom_processes.keys()
        for process_uuid in process_uuids:
            logger.info('Finishing process {}...'.format(process_uuid))
            self.finish_loom_process(process_uuid, join_timeout=0)
            try:
                # Flush the result queue.
                result_queue = self.result_queues[process_uuid]
                while result_queue.empty() is False:
                    result_queue.get_nowait()
                # Remove the result queue
                del self.result_queues[process_uuid]
            except KeyError:
                logger.warning(
                    "Removing result queue {} failed: "
                    "no such a queue exists."
                    .format(process_uuid)
                ) 
                pass
        logger.info('LoomDB manager thread is finished.')

    def start_loom_process(
        self, process_uuid, logging_level, loom_config, n_processes=None,
    ):
        if n_processes is None:
            n_processes = DEFAULT_NUM_PROCESSES

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
            target=generate_spectral_network,
            args=(loom_config,),
            kwargs=dict(
                n_processes=n_processes,
                result_queue=result_queue,
                logging_queue=logging_queue,
                logger_name=logger_name,
            ),
        )
        self.loom_processes[process_uuid] = loom_process
        self.is_alive[process_uuid] = True
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
        stop_finish = False

        try:
            result_queue = self.result_queues[process_uuid]
        except KeyError:
            yield 'event: key_error\ndata: \n\n'
            raise StopIteration

        while result_queue.empty() is True:
            try:
                logs = self.get_log_message(process_uuid, logging_stream,
                                            logging_stream_handler,)
                yield 'data: {}\n\n'.format(logs)
            except QueueEmpty:
                pass 
            except (KeyboardInterrupt, SystemExit):
                raise
            except KeyError:
                yield 'event: key_error\ndata: \n\n'
                raise StopIteration
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
            except KeyError:
                yield 'event: finish\ndata: \n\n'
                raise StopIteration

        # Recevied all the logs, finish the SSE stream.
        yield 'event: finish\ndata: \n\n'
        raise StopIteration

    def get_result(self, process_uuid):
        logger_name = get_logger_name(process_uuid)
        logger = logging.getLogger(logger_name)

        rv = None

        result_queue = self.result_queues[process_uuid]
        loom_process = self.loom_processes[process_uuid]

        if result_queue.empty() is True:
            if loom_process.is_alive():
                logger.warning('Process {} still alive.'.format(process_uuid))
                return rv 
            else:
                logger.warning(
                    'Generating spectral networks failed: '
                    'pid = {}, exitcode = {}.'
                    .format(loom_process.pid,
                            loom_process.exitcode,)
                )
                return rv 
        else:
            # Result queue has the data returned from the loom_process.
            rv = result_queue.get()
            logger.info('Finished generating spectral network data.')

            # Start a thread to record a stat log.
            user_ip = flask.request.remote_addr
            # XXX: Use the following behind a Proxy server.
            # user_ip = flask.request.environ.get('HTTP_X_REAL_IP',
            #                                     request.remote_addr)
            stat_thread = threading.Thread(
                target=record_stat,
                args=(rv, STAT_LOGGER_NAME, user_ip, process_uuid),
            )
            stat_thread.start()

            self.finish_loom_process(process_uuid)
            return rv 

    def finish_loom_process(self, process_uuid,
                            join_timeout=LOOM_PROCESS_JOIN_TIMEOUT_SECS):
        logger_name = get_logger_name(process_uuid)
        logger = logging.getLogger(logger_name)
        web_loom_logger = logging.getLogger(get_logger_name())

        try:
            # Terminate the loom_process.
            loom_process = self.loom_processes[process_uuid]

            loom_process.join(join_timeout)
            if loom_process.is_alive():
                loom_process.terminate()

            del self.loom_processes[process_uuid]
        except KeyError:
            web_loom_logger.warning(
                "Terminating loom_process {} failed: no such a process exists."
                .format(logger_name)
            ) 
            pass

        try:
            # Remove the logging queue handler.
            logger.handlers = []
            # Remove the logger.
            del logging.Logger.manager.loggerDict[logger_name]
        except KeyError:
            web_loom_logger.warning(
                "Removing logger {} failed: no such a logger exists."
                .format(logger_name)
            ) 
            pass

        try:
            # Flush the logging queue.
            logging_queue = self.logging_queues[process_uuid]
            while logging_queue.empty() is True:
                try:
                    logging_queue.get_nowait()
                except QueueEmpty:
                    break
            # Remove the logging queue.
            del self.logging_queues[process_uuid]
        except KeyError:
            web_loom_logger.warning(
                "Removing logging queue {} failed: no such a queue exists."
                .format(process_uuid)
            ) 
            pass


class WebLoomApplication(flask.Flask):
    """
    A wrapper of Flask application containing an instance of LoomDB.
    """
    def __init__(self, config_file, logging_level):
        super(WebLoomApplication, self).__init__(WEB_APP_NAME)
        # Set a logger for loom that has a rotating file handler.
        set_logging(
            logger_name=WEB_APP_NAME,
            logging_level=logging_level,
            logging_file_name=get_logging_file_path(WEB_APP_NAME),
            use_rotating_file_handler=True,
        )
        # Set a logger for recording the stat
        # with a non-rotating file handler.
        set_logging(
            logger_name=STAT_LOGGER_NAME,
            logging_file_name=get_logging_file_path(STAT_LOGGER_NAME),
        )
        self.loom_db = LoomDB()


###
# View functions
###


def index():
    # Make a list of contributors from git logs.
    ps = subprocess.Popen(
        ['git', 'log', '--date-order', '--reverse', '--format="%aN"'], 
        stdout=subprocess.PIPE,
    )
    loom_contributors_str = subprocess.check_output(
        ['awk', '!x[$0]++'],
        stdin=ps.stdout
    ).strip().split("\n")
    ps.wait()

    loom_contributors = []
    for name_str in loom_contributors_str:
        name = name_str.strip('"')
        if name == 'plonghi':
            loom_contributors.append('Pietro Longhi')
        elif name == 'chan':
            if 'Chan Y. Park' not in loom_contributors:
                loom_contributors.append('Chan Y. Park')
        else:
            loom_contributors.append(name)

    return flask.render_template(
        'index.html',
        loom_contributors=loom_contributors,     
    )


def config(n_processes=None):
    # XXX: n_processes is a string.
    loom_config = None
    event_source_url = None
    text_area_content = '' 
    process_uuid = None

    if flask.request.method == 'GET':
        # Load the default configuration.
        loom_config = get_loom_config()

    elif flask.request.method == 'POST':
        try:
            uploaded_config_file = flask.request.files['config_file']
        except KeyError:
            uploaded_config_file = None

        if uploaded_config_file is not None:
            if uploaded_config_file.filename == '':
                # Load button clicked without a selected file.
                # Load the default configuration.
                loom_config = get_loom_config()
            else:
                loom_config = LoomConfig(logger_name=get_logger_name())
                loom_config.read(uploaded_config_file)
        else:
            if n_processes is None:
                n_processes = flask.request.form['n_processes']
            process_uuid = str(uuid.uuid4())
            logger_name = get_logger_name(process_uuid)
            loom_config = get_loom_config(flask.request.form, logger_name) 

            app = flask.current_app
            app.loom_db.start_loom_process(
                process_uuid, logging.INFO, loom_config,
                n_processes=eval(n_processes),
            )
            event_source_url = flask.url_for(
                'logging_stream', process_uuid=process_uuid,
            )
            text_area_content = (
                "Start loom, uuid = {}".format(process_uuid)
            )

    return flask.render_template(
        'config.html',
        config_options=config_options,
        advanced_config_options=advanced_config_options,
        loom_config=loom_config,
        n_processes=n_processes,
        process_uuid=process_uuid,
        event_source_url=event_source_url,
        text_area_content=text_area_content,
    )


def logging_stream(process_uuid):
    if flask.request.headers.get('accept') == 'text/event-stream':
        app = flask.current_app
        return flask.Response(
            app.loom_db.yield_log_message(process_uuid, logging.INFO),
            mimetype='text/event-stream',
        )


def save_config():
    loom_config = get_loom_config(flask.request.form)
    loom_config_fp = BytesIO()
    loom_config.parser.write(loom_config_fp)
    loom_config_fp.seek(0)
    rv = flask.send_file(loom_config_fp, mimetype='text/plain',
                         as_attachment=True,
                         attachment_filename='config.ini',
                         add_etags=True,)
    return rv
        

def plot():
    loom_db = flask.current_app.loom_db
    rotate_back = False

    if flask.request.method == 'POST':
        try:
            rotate_back = bool(flask.request.form['rotate_back'])
        except KeyError:
            pass
        process_uuid = flask.request.form['process_uuid']
        progress_log = flask.request.form['progress_log']

        if rotate_back is True:
            rv = loom_db.result_queues[process_uuid].get()
        else:
            # Finish loom_process
            rv = loom_db.get_result(process_uuid)

    elif flask.request.method == 'GET':
        data_dir = flask.request.args['data']
        process_uuid = str(uuid.uuid4()) 
        progress_log = None
        full_data_dir = os.path.join(
            get_loom_dir(), 'data', data_dir
        )
        rv = load_spectral_network(
            full_data_dir, logger_name=get_logger_name()
        )
        loom_db.result_queues[process_uuid] = multiprocessing.Queue()

    loom_config, spectral_network_data = rv

    if rotate_back is True:
        #reset_z_rotation = False
        bc_r = spectral_network_data.sw_data.branch_cut_rotation
        print 'branch_cut_rotation = {}'.format(bc_r)
        spectral_network_data.set_z_rotation(1/bc_r)
        for p in spectral_network_data.sw_data.regular_punctures:
            print '{}: z = {}'.format(p.label, p.z)
    else:
        spectral_network_data.reset_z_rotation()
        #reset_z_rotation = True
 
    # Put data back into the queue for future use.
    loom_db.result_queues[process_uuid].put(rv)

    return render_plot_template(
        loom_config, spectral_network_data, process_uuid=process_uuid,
        progress_log=progress_log, #reset_z_rotation=reset_z_rotation,
    )


def save_data_to_server():
    if flask.request.method == 'POST':
        process_uuid = flask.request.form['process_uuid']
        data_name = flask.request.form['data_name']

        logger_name = get_logger_name(process_uuid)
        logger = logging.getLogger(logger_name)

        loom_db = flask.current_app.loom_db
        rv = loom_db.result_queues[process_uuid].get()
        loom_config, spectral_network_data = rv

        data_dir = os.path.join(get_loom_dir(), 'data', data_name,)
        # XXX: check if there is an existing dir with the same name.
        if os.path.exists(data_dir):
            msg = (
                'Data with name "{}" already exists, '
                'chose a different name.'.format(data_name)
            )
        else:
            save_spectral_network(
                loom_config,
                spectral_network_data,
                data_dir=data_dir,
                logger_name=logger_name,
            )
            msg = 'Data successfully saved as "{}".'.format(data_name)

    return flask.render_template(
        'save_message.html',
        msg=msg,
    )


def download_data(process_uuid):
    loom_db = flask.current_app.loom_db
    rv = loom_db.result_queues[process_uuid].get()
    loom_config, spectral_network_data = rv
    sw_data = spectral_network_data.sw_data
    spectral_networks = spectral_network_data.spectral_networks
    data = {}
    
    data['version'] = get_current_branch_version()

    # fp = StringIO()
    fp = BytesIO()
    loom_config.parser.write(fp)
    fp.seek(0)
    data['config.ini'] = fp.read()

    data['sw_data.json'] = json.dumps(sw_data.get_json_data())

    for i, spectral_network in enumerate(spectral_networks):
        file_name_idx = str(i).zfill(len(str(len(spectral_networks) - 1)))
        data['data_{}.json'.format(file_name_idx)] = json.dumps(
            spectral_network.get_json_data()
        )

    data_zip_fp = BytesIO()
    with zipfile.ZipFile(data_zip_fp, 'w') as zfp:
        for file_name, data_str in data.iteritems():
            zip_info = zipfile.ZipInfo(file_name)
            zip_info.date_time = time.localtime(time.time())[:6]
            zip_info.compress_type = zipfile.ZIP_DEFLATED
            zip_info.external_attr = 0777 << 16L
            zfp.writestr(zip_info, data_str)
    data_zip_fp.seek(0)

    loom_db.result_queues[process_uuid].put(rv)

    return flask.send_file(
        data_zip_fp,
        attachment_filename='loom_data_{}.zip'.format(process_uuid),
        as_attachment=True,
    )


def download_plot(process_uuid):
    loom_db = flask.current_app.loom_db

    rv = loom_db.result_queues[process_uuid].get()
    loom_config, spectral_network_data = rv
    loom_db.result_queues[process_uuid].put(rv)

    plot_html_zip_fp = BytesIO()
    with zipfile.ZipFile(plot_html_zip_fp, 'w') as zfp:
        zip_info = zipfile.ZipInfo('loom_plot_{}.html'.format(process_uuid))
        zip_info.date_time = time.localtime(time.time())[:6]
        zip_info.compress_type = zipfile.ZIP_DEFLATED
        zip_info.external_attr = 0777 << 16L
        zfp.writestr(
            zip_info,
            render_plot_template(loom_config, spectral_network_data,
                                 download=True,),
        )
    plot_html_zip_fp.seek(0)

    return flask.send_file(
        plot_html_zip_fp,
        attachment_filename='loom_plot_{}.html.zip'.format(process_uuid),
        as_attachment=True,
    )


def download_E6_E7_data():
    return flask.send_file(
        'data/E6_E7_data.tar.gz',
        as_attachment=True,
    )


def keep_alive(process_uuid):
    """
    Receive heartbeats from clients.
    """
    logger = logging.getLogger(get_logger_name())
    app = flask.current_app
    try:
        is_alive = app.loom_db.is_alive
        if is_alive[process_uuid] is False:
            is_alive[process_uuid] = True
        logger.debug(
            'is_alive[{}] = {}'
            .format(process_uuid, app.loom_db.is_alive[process_uuid])
        )
    except KeyError:
        pass
    return ('', 204)
        

def admin():
    # TODO: password-protect this page.
    # app = flask.current_app
    # loom_db = app.loom_db
    # pdb.set_trace()
    return ('', 204)


def shutdown():
    app = flask.current_app
    app.loom_db.db_manager_stop.set()

    f = flask.request.environ.get('werkzeug.server.shutdown')
    if f is None:
        raise RuntimeError('Not running with the Werkzeug Server')
    f()
    return 'Server shutting down...'


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
        '/config/<n_processes>', 'config', config, methods=['GET', 'POST'],
    )
    application.add_url_rule(
        '/config', 'config', config, methods=['GET', 'POST'],
    )
    application.add_url_rule(
        '/save_config', 'save_config', save_config, methods=['POST'],
    )
    application.add_url_rule(
        '/plot', 'plot', plot,
        methods=['GET', 'POST'],
    )
    application.add_url_rule(
        '/save_data_to_server', 'save_data_to_server',
        save_data_to_server,
        methods=['POST'],
    )
    application.add_url_rule(
        '/download_data/<process_uuid>', 'download_data', download_data,
        methods=['POST'],
    )
    application.add_url_rule(
        '/download_plot/<process_uuid>', 'download_plot', download_plot,
        methods=['POST'],
    )
    application.add_url_rule(
        '/logging_stream/<process_uuid>', 'logging_stream', logging_stream,
        methods=['GET'],
    )
    application.add_url_rule(
        '/keep_alive/<process_uuid>', 'keep_alive', keep_alive,
        methods=['GET'],
    )
    application.add_url_rule(
        '/E6_E7_data', 'download_E6_E7_data', download_E6_E7_data,
        methods=['GET'],
    )
    application.add_url_rule(
        '/admin', 'admin', admin,
        methods=['GET'],
    )
    application.add_url_rule(
        '/shutdown', 'shutdown', shutdown,
        methods=['GET'],
    )
    return application


###
# Misc. web UIs
###


def get_logger_name(uuid=None):
    logger_name = WEB_APP_NAME
    if uuid is not None:
        logger_name += '.' + uuid
    return logger_name


def get_loom_config(request_dict=None, logger_name=get_logger_name()):
    logger = logging.getLogger(logger_name)

    default_config_file = os.path.join(
        get_loom_dir(),
        'config/default.ini',
    )
    loom_config = load_config(default_config_file, logger_name=logger_name)

    if request_dict is not None:
        # Update config with form data.
        root_system = request_dict['type'] + request_dict['rank']
        for section in loom_config.parser.sections():
            for option in loom_config.parser.options(section):
                try:
                    if option == 'root_system':
                        value = root_system
                    else:
                        value = request_dict[option]
                    if (section == 'numerical parameters'):
                        loom_config[option] = eval(value)
                    else:
                        loom_config[option] = value
                    loom_config.parser.set(section, option, value)
                except KeyError:
                    logger.warning(
                        'No entry for option = {}, skip it.'
                        .format(option)
                    )
                    pass

    return loom_config


def render_plot_template(
    loom_config, spectral_network_data, process_uuid=None,
    progress_log=None, download=False, #rotate_back=False,
):
    download_data_url = download_plot_url = None
    sw_data = spectral_network_data.sw_data

    # Make a Bokeh plot
    bokeh_plot_script, div = get_spectral_network_bokeh_plot(
        spectral_network_data,
        plot_range=loom_config['plot_range'],
        logger_name=get_logger_name(),
        #reset_z_rotation=reset_z_rotation,
    )

    initial_phase = '{:.3f}'.format(
        spectral_network_data.spectral_networks[0].phase / pi
    )

    legend = get_legend(
        g_data=sw_data.g_data,
        regular_punctures=sw_data.regular_punctures,
        branch_points=sw_data.branch_points,
        irregular_singularities=sw_data.irregular_singularities,
    )

    if download is False:
        download_data_url = flask.url_for(
            'download_data',
            process_uuid=process_uuid,
        )
        download_plot_url = flask.url_for(
            'download_plot',
            process_uuid=process_uuid,
        )

    with open('static/bokeh_callbacks.js', 'r') as fp:
        bokeh_custom_script = fp.read()

    return flask.render_template(
        'plot.html',
        process_uuid=process_uuid,
        bokeh_plot_script=bokeh_plot_script,
        div=div,
        progress_log=progress_log,
        plot_legend=legend,
        bokeh_custom_script=bokeh_custom_script,
        download_data_url=download_data_url,
        download_plot_url=download_plot_url,
        loom_config=loom_config,
        config_options=config_options,
        advanced_config_options=advanced_config_options,
        initial_phase=initial_phase,
    )


def get_logging_file_path(logger_name):
    logging_file_path = os.path.join(
        get_loom_dir(),
        ('logs/{}_{}-{:02}-{:02} {:02}:{:02}:{:02}.log'
         .format(logger_name, *time.localtime(time.time())[:6])),
    )
    return logging_file_path


def record_stat(rv, stat_logger_name, ip, uuid):
    config, spectral_network_data = rv
    # Make a zipped config file and record the stat.
    config_file_name = '{}.ini'.format(uuid)
    zipfile_path = os.path.join(
        get_loom_dir(), 'logs', '{}.zip'.format(config_file_name),
    )
    config_fp = BytesIO()
    config.parser.write(config_fp)
    config_fp.seek(0)
    with zipfile.ZipFile(zipfile_path, 'w') as zfp:
        zip_info = zipfile.ZipInfo(config_file_name)
        zip_info.date_time = time.localtime(time.time())[:6]
        zip_info.compress_type = zipfile.ZIP_DEFLATED
        #zip_info.external_attr = 0777 << 16L
        zip_info.external_attr = 040755 << 16L
        zfp.writestr(zip_info, config_fp.read())

    data_size = get_data_size_of(spectral_network_data)

    stat_logger = logging.getLogger(stat_logger_name)
    stat_logger.info('{}, {}, {}'.format(ip, uuid, data_size))
