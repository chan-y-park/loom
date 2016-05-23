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
import zipfile
import glob
import shutil

from cStringIO import StringIO
from io import BytesIO
from Queue import Empty as QueueEmpty

from api import (
    SpectralNetworkData,
    get_loom_dir,
    get_logging_handler,
    set_logging,
    load_config,
)
from config import LoomConfig
from bokeh_plot import get_spectral_network_bokeh_plot
from plotting import get_legend
from misc import get_data_size_of
from misc import get_phases_from_dict
from misc import get_phase_dict
from cmath import pi

# Flask configuration
DEBUG = True
SECRET_KEY = 'web_loom_key'
# Global variables
PARENT_LOGGER_NAME = 'loom'
WEB_APP_NAME = 'web_loom'
STAT_LOGGER_NAME = 'stat_of_' + WEB_APP_NAME
DEFAULT_NUM_PROCESSES = 4
DB_CLEANUP_CYCLE_SECS = 3600
LOOM_PROCESS_JOIN_TIMEOUT_SECS = 3
LOOM_CACHE_DIR = 'cache'

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
    ['ramification_points'],
    ['branch_points'],
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
    def __init__(self, logging_level=None):
        self.logging_queues = {}
        self.result_queues = {}
        self.loom_processes = {}
        self.logging_level = logging_level

        signal.signal(signal.SIGINT, self.loom_db_stop_signal_handler)

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
            pass

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
        self,
        process_uuid=None,
        loom_config=None,
        n_processes=None,
        full_data_dir=None,
        task=None,
        data_name=None,
        saved_data=None,
        additional_n_steps=0, new_mass_limit=None,
        additional_iterations=0, additional_phases=None,
    ):
        logging_level = self.logging_level
        if task == 'generate' or task == 'extend':
            # Prepare a cache directory for new data.
            cache_dir = get_cache_dir(process_uuid)
            if os.path.exists(cache_dir) is False:
                os.makedirs(cache_dir)
            logging_file_name=os.path.join(cache_dir, 'log'),
        elif task == 'load':
            # Do not create a logging file
            # when loading a saved data.
            logging_file_name=None
        else:
            raise RuntimeError('Unknown task for loom: {}'.format(task))

        if n_processes is None:
            n_processes = DEFAULT_NUM_PROCESSES

        logging_queue = multiprocessing.Queue()
        self.logging_queues[process_uuid] = logging_queue
        logger_name = get_logger_name(process_uuid)
        set_logging(
            logger_name=logger_name,
            logging_level=logging_level,
            logging_queue=logging_queue,
            logging_file_name=logging_file_name,
        )

        result_queue = multiprocessing.Queue()
        self.result_queues[process_uuid] = result_queue

        if task == 'load':
            spectral_network_data = SpectralNetworkData(
                logger_name=logger_name,
            )
            loom_process = multiprocessing.Process(
                target=spectral_network_data.load,
                kwargs=dict(
                    data_dir=get_full_data_dir(data_name, saved_data),
                    result_queue=result_queue,
                    logging_queue=logging_queue,
                ),
            )

#        elif (
#            additional_n_steps == 0 and
#            additional_iterations == 0 and
#            new_mass_limit is None and
#            additional_phases is None
#        ):
        elif task == 'generate':
            spectral_network_data = SpectralNetworkData(
                config=loom_config,
                logger_name=logger_name,
            )

            loom_config['phase'] = get_phase_dict(loom_config['phase'])

            phases = get_phases_from_dict(
                loom_config['phase'], loom_config['accuracy'],
            )

            loom_process = multiprocessing.Process(
                target=spectral_network_data.generate,
                kwargs=dict(
                    phases=phases,
                    n_processes=n_processes,
                    result_queue=result_queue,
                    logging_queue=logging_queue,
                    cache_dir=cache_dir,
                ),
            )

        elif task == 'extend':
            full_data_dir = get_full_data_dir(data_name, saved_data)
            spectral_network_data = SpectralNetworkData(
                data_dir=full_data_dir,
                logger_name=logger_name,
            )
            loom_process = multiprocessing.Process(
                target=spectral_network_data.extend,
                kwargs=dict(
                    additional_n_steps=additional_n_steps,
                    new_mass_limit=new_mass_limit,
                    additional_iterations=additional_iterations,
                    additional_phases=additional_phases,
                    n_processes=n_processes,
                    result_queue=result_queue,
                    logging_queue=logging_queue,
                    cache_dir=cache_dir,
                )
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

        result_queue = self.result_queues[process_uuid]
        loom_process = self.loom_processes[process_uuid]

        if result_queue.empty() is True:
            if loom_process.is_alive():
                logger.warning('Process {} still alive.'.format(process_uuid))
                return None
            else:
                logger.warning(
                    'Generating spectral networks failed: '
                    'pid = {}, exitcode = {}.'
                    .format(loom_process.pid,
                            loom_process.exitcode,)
                )
                return None
        else:
            # Result queue has the data returned from the loom_process.
            spectral_network_data = result_queue.get()
            logger.info(
                'Process {} finished generating spectral network data.'
                .format(process_uuid)
            )

            # Start a thread to record a stat log.
            user_ip = flask.request.remote_addr
            # XXX: Use the following behind a Proxy server.
            # user_ip = flask.request.environ.get('HTTP_X_REAL_IP',
            #                                     request.remote_addr)
            # XXX: Use the file sizes instead by checking the directory
            # using the process_uuid. Do it inside the thread.
            data_size = get_data_size_of(spectral_network_data)
            stat_thread = threading.Thread(
                target=record_stat,
                args=(STAT_LOGGER_NAME, user_ip, process_uuid, data_size),
            )
            stat_thread.start()

            self.finish_loom_process(process_uuid)
            return spectral_network_data

    def finish_loom_process(
        self, process_uuid,
        join_timeout=LOOM_PROCESS_JOIN_TIMEOUT_SECS
    ):
        logger_name = get_logger_name(process_uuid)
        logger = logging.getLogger(logger_name)
        web_loom_logger = logging.getLogger(get_logger_name())

        try:
            # Terminate the loom_process.
            loom_process = self.loom_processes[process_uuid]

            loom_process.join(join_timeout)
            if loom_process.is_alive():
                web_loom_logger.warning(
                    'Process {} did not join successfully after timeout.'
                    .format(process_uuid)
                )
                loom_process.terminate()
            else:
                web_loom_logger.info(
                    'Process {} joined successfully.'
                    .format(process_uuid)
                )

            del self.loom_processes[process_uuid]
        except KeyError:
            web_loom_logger.warning(
                "Terminating loom_process {} failed: no such a process exists."
                .format(process_uuid)
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

        try:
            # Flush the result queue.
            result_queue = self.result_queues[process_uuid]
            while result_queue.empty() is True:
                try:
                    result_queue.get_nowait()
                except QueueEmpty:
                    break
            # Remove the result queue.
            del self.result_queues[process_uuid]
        except KeyError:
            web_loom_logger.warning(
                "Removing result queue {} failed: no such a queue exists."
                .format(process_uuid)
            )
            pass


class WebLoomApplication(flask.Flask):
    """
    A wrapper of Flask application containing an instance of LoomDB.
    """
    def __init__(self, config_file, logging_level):
        super(WebLoomApplication, self).__init__(WEB_APP_NAME)
        # Set a logger for loom web frontend.
        set_logging(
            logger_name=WEB_APP_NAME,
            logging_level=logging_level,
            logging_file_name=get_logging_file_path(WEB_APP_NAME),
        )
        # Set a logger for recording the stat
        # with a non-rotating file handler.
        set_logging(
            logger_name=STAT_LOGGER_NAME,
            logging_level=logging_level,
            logging_file_name=get_logging_file_path(STAT_LOGGER_NAME),
        )
        self.loom_db = LoomDB(logging_level)


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
    full_data_dir = None

    if flask.request.method == 'GET':
        # Load the default configuration.
        loom_config = get_loom_config()
        if n_processes is None:
            n_processes_val = DEFAULT_NUM_PROCESSES
        else:
            n_processes_val = eval(n_processes)

    elif flask.request.method == 'POST':
        if n_processes is None:
            n_processes_val = eval(flask.request.form['n_processes'])
        else:
            n_processes_val = eval(n_processes)

        try:
            uploaded_config_file = flask.request.files['config_file']
        except KeyError:
            uploaded_config_file = None

        try:
            data_name = flask.request.form['data_name']
        except KeyError:
            data_name = None

        if uploaded_config_file is not None:
            # Load configuration from the uploaded file.
            if uploaded_config_file.filename == '':
                # Load button clicked without a selected file.
                # Load the default configuration.
                loom_config = get_loom_config()
            else:
                loom_config = LoomConfig(logger_name=get_logger_name())
                loom_config.read(uploaded_config_file)
        elif data_name is not None:
            config_file_path = os.path.join(
                get_loom_dir(),
                'data',
                data_name,
                'config.ini',
            )
            loom_config = LoomConfig(
                file_path=config_file_path,
                logger_name=get_logger_name(),
            )
#        else:
#            # Generate/Extend spectral networks.
#
#            additional_params = {
#                'additional_n_steps': 0,
#                'new_mass_limit': None,
#                'additional_iterations': 0,
#                'additional_phases': None
#            }
#
#            for key in additional_params.keys():
#                try:
#                    additional_params[key] = eval(
#                        flask.request.form[key]
#                    )
#                except (KeyError, SyntaxError):
#                    pass
#
#            logger_name = get_logger_name(process_uuid)
#
#            if (
#                additional_params['additional_n_steps'] == 0 and
#                additional_params['additional_iterations'] == 0 and
#                additional_params['new_mass_limit'] is None and
#                additional_params['additional_phases'] is None
#            ):
#                loom_config = get_loom_config(flask.request.form, logger_name)
#            else:
#                loom_config = None
#                try:
#                    process_uuid = flask.request.form['process_uuid']
#                    saved_data = eval(flask.request.form['saved_data'])
#                    full_data_dir = get_full_data_dir(process_uuid, saved_data)
#                except KeyError:
#                    raise RuntimeError(
#                        'Need data to extend spectral networks.'
#                    )
#
#            process_uuid = str(uuid.uuid4())
#
#            app = flask.current_app
#            app.loom_db.start_loom_process(
#                process_uuid=process_uuid,
#                loom_config=loom_config,
#                n_processes=n_processes_val,
#                full_data_dir=full_data_dir,
#                **additional_params
#            )
#            event_source_url = flask.url_for(
#                'logging_stream', process_uuid=process_uuid,
#            )
#            text_area_content = (
#                "Start loom, uuid = {}".format(process_uuid)
#            )

    return flask.render_template(
        'config.html',
        config_options=config_options,
        advanced_config_options=advanced_config_options,
        loom_config=loom_config,
        n_processes=n_processes_val,
        process_uuid=process_uuid,
        event_source_url=event_source_url,
        text_area_content=text_area_content,
    )


def load(n_processes=None):
    # XXX: n_processes is a string.
    loom_config = None
    event_source_url = None
    text_area_content = ''
    process_uuid = None
    #full_data_dir = None

    if n_processes is None:
        n_processes_val = DEFAULT_NUM_PROCESSES
    else:
        n_processes_val = eval(n_processes)

    if flask.request.method == 'GET':
        full_data_directories = glob.glob(
            os.path.join(get_loom_dir(), 'data', "*",)
        )
        full_data_directories.sort()
        data_names = [
            os.path.split(full_data_dir)[1]
            for full_data_dir in full_data_directories
        ]
    elif flask.request.method == 'POST':
        print('POST')

    return flask.render_template(
        'load.html',
        data_names=data_names,
        n_processes=n_processes_val,
    )


def progress():
    event_source_url = None
    text_area_content = ''
    process_uuid = None
    saved_data = None
    n_processes_val = eval(flask.request.form['n_processes'])

    try:
        data_name = flask.request.form['data_name']
    except KeyError:
        data_name = None

    if data_name is not None:
        # Load a saved data and plot it.
        saved_data = True

        try:
            n_processes = flask.request.args['n_processes']
        except KeyError:
            n_processes = None

        process_uuid = str(uuid.uuid4())
    else:
        additional_params = {
            'additional_n_steps': 0,
            'new_mass_limit': None,
            'additional_iterations': 0,
            'additional_phases': None
        }

        for key in additional_params.keys():
            try:
                additional_params[key] = eval(
                    flask.request.form[key]
                )
            except (KeyError, SyntaxError):
                pass

        if (
            additional_params['additional_n_steps'] == 0 and
            additional_params['additional_iterations'] == 0 and
            additional_params['new_mass_limit'] is None and
            additional_params['additional_phases'] is None
        ):
            # Generate a new spectral network.
            saved_data = False
            process_uuid = str(uuid.uuid4())
            logger_name = get_logger_name(process_uuid)
            loom_config = get_loom_config(flask.request.form, logger_name)
        else:
            # Extend a spectral network.
            loom_config = None
            try:
                process_uuid = flask.request.form['process_uuid']
                saved_data = eval(flask.request.form['saved_data'])
                full_data_dir = get_full_data_dir(process_uuid, saved_data)
            except KeyError:
                raise RuntimeError(
                    'Need data to extend spectral networks.'
                )
            process_uuid = str(uuid.uuid4())

    app = flask.current_app
    app.loom_db.start_loom_process(
        process_uuid=process_uuid,
        loom_config=loom_config,
        n_processes=n_processes_val,
        full_data_dir=full_data_dir,
        **additional_params
    )

    event_source_url = flask.url_for(
        'logging_stream', process_uuid=process_uuid,
    )
    text_area_content = (
        "Start loom, uuid = {}".format(process_uuid)
    )

    return flask.render_template(
        'progress.html',
        n_processes=n_processes_val,
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
    rotate_back = None
    plot_two_way_streets = None
    saved_data = None
    search_radius = None

    if flask.request.method == 'POST':
        try:
            rotate_back = eval(flask.request.form['rotate_back'])
        except KeyError:
            rotate_back = False

        try:
            plot_two_way_streets = eval(
                flask.request.form['plot_two_way_streets']
            )
            search_radius_str = flask.request.form['search_radius']
            if search_radius_str != '':
                search_radius = eval(search_radius_str)
        except KeyError:
            plot_two_way_streets = False
            search_radius = None

        process_uuid = flask.request.form['process_uuid']
        progress_log = flask.request.form['progress_log']
        n_processes = flask.request.form['n_processes']
        saved_data = eval(flask.request.form['saved_data'])

#        if rotate_back is True or plot_two_way_streets is True:
#            full_data_dir = get_full_data_dir(process_uuid, saved_data)
#
#            spectral_network_data = SpectralNetworkData(
#                data_dir=full_data_dir,
#            )
#        else:
#            # Finish loom_process
#            spectral_network_data = loom_db.get_result(process_uuid)
        spectral_network_data = loom_db.get_result(process_uuid)

    elif flask.request.method == 'GET':
        # Load saved data.
        saved_data = True
        data_dir = flask.request.args['data']
        try:
            n_processes = flask.request.args['n_processes']
        except KeyError:
            n_processes = None
        try:
            plot_two_way_streets = eval(
                flask.request.args['plot_two_way_streets']
            )
        except KeyError:
            plot_two_way_streets = False
        try:
            search_radius_str = flask.request.args['search_radius']
            if search_radius_str != '':
                search_radius = eval(search_radius_str)
        except KeyError:
            search_radius = None

        process_uuid = data_dir
        progress_log = None
        full_data_dir = get_full_data_dir(process_uuid, saved_data)
        spectral_network_data = SpectralNetworkData(
            data_dir=full_data_dir,
        )

    if rotate_back is True:
        spectral_network_data.rotate_back()
    else:
        spectral_network_data.reset_z_rotation()

    return render_plot_template(
        spectral_network_data,
        process_uuid=process_uuid,
        progress_log=progress_log,
        n_processes=n_processes,
        saved_data=saved_data,
        plot_two_way_streets=plot_two_way_streets,
        search_radius=search_radius
    )


def save_data_to_server():
    if flask.request.method == 'POST':
        process_uuid = flask.request.form['process_uuid']
        saved_data = eval(flask.request.form['saved_data'])
        data_name = flask.request.form['data_name']
    else:
        raise RuntimeError

    data_dir_to_save = os.path.join(get_loom_dir(), 'data', data_name,)
    if os.path.exists(data_dir_to_save):
        msg = (
            'Data with name "{}" already exists, '
            'chose a different name.'.format(data_name)
        )
    else:
        os.makedirs(data_dir_to_save)
        full_data_dir = get_full_data_dir(process_uuid, saved_data)
        files_to_copy = get_data_file_path_list(full_data_dir)
        for src in files_to_copy:
            shutil.copy(src, data_dir_to_save)
        msg = 'Data successfully saved as "{}".'.format(data_name)

    return flask.render_template(
        'save_message.html',
        msg=msg,
    )


def download_data():
    if flask.request.method == 'POST':
        process_uuid = flask.request.form['process_uuid']
        saved_data = eval(flask.request.form['saved_data'])
    else:
        raise RuntimeError

    full_data_dir = get_full_data_dir(process_uuid, saved_data)
    files_to_zip = get_data_file_path_list(full_data_dir)

    data_zip_fp = BytesIO()
    with zipfile.ZipFile(data_zip_fp, 'w') as zfp:
        for file_path in files_to_zip:
            arcname = os.path.basename(file_path)
            zfp.write(file_path, arcname, compress_type=zipfile.ZIP_DEFLATED)
    data_zip_fp.seek(0)

    return flask.send_file(
        data_zip_fp,
        attachment_filename='loom_data_{}.zip'.format(process_uuid),
        as_attachment=True,
    )


def download_plot():
    if flask.request.method == 'POST':
        process_uuid = flask.request.form['process_uuid']
        saved_data = eval(flask.request.form['saved_data'])
        plot_two_way_streets = eval(
            flask.request.form['plot_two_way_streets']
        )
        search_radius_str = flask.request.form['search_radius']
        if search_radius_str != '':
            search_radius = eval(search_radius_str)
        else:
            search_radius = None
    else:
        raise RuntimeError

    full_data_dir = get_full_data_dir(process_uuid, saved_data)
    spectral_network_data = SpectralNetworkData(data_dir=full_data_dir)
    spectral_network_data.reset_z_rotation()

    plot_html_zip_fp = BytesIO()
    with zipfile.ZipFile(plot_html_zip_fp, 'w') as zfp:
        zip_info = zipfile.ZipInfo('loom_plot_{}.html'.format(process_uuid))
        zip_info.date_time = time.localtime(time.time())[:6]
        zip_info.compress_type = zipfile.ZIP_DEFLATED
        zip_info.external_attr = 040664 << 16L
        zfp.writestr(
            zip_info,
            render_plot_template(
                spectral_network_data,
                process_uuid=process_uuid,
                saved_data=saved_data,
                download=True,
                plot_two_way_streets=plot_two_way_streets,
                search_radius=search_radius
            ),
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
        '/load/<n_processes>', 'load', load, methods=['GET', 'POST'],
    )
    application.add_url_rule(
        '/load', 'load', load, methods=['GET', 'POST'],
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
        '/download_data', 'download_data', download_data,
        methods=['POST'],
    )
    application.add_url_rule(
        '/download_plot', 'download_plot', download_plot,
        methods=['POST'],
    )
    application.add_url_rule(
        '/logging_stream/<process_uuid>', 'logging_stream', logging_stream,
        methods=['GET'],
    )
    application.add_url_rule(
        '/E6_E7_data', 'download_E6_E7_data', download_E6_E7_data,
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
#    loom_config = load_config(default_config_file, logger_name=logger_name)
    loom_config = LoomConfig(
        file_path=default_config_file,
        logger_name=logger_name,
    )

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
                    if (
                        section == 'numerical parameters'
                        or value == 'None'
                    ):
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
    spectral_network_data, process_uuid=None,
    progress_log=None, n_processes=None,
    download=False, saved_data=False,
    plot_two_way_streets=False, search_radius=None,
):
    loom_config = spectral_network_data.config
    sw_data = spectral_network_data.sw_data
    soliton_tree_data = None

    if plot_two_way_streets is True:
        soliton_tree_data = spectral_network_data.find_two_way_streets(
            search_radius=search_radius,
        )

    # Make a Bokeh plot
    bokeh_plot_script, div = get_spectral_network_bokeh_plot(
        spectral_network_data,
        plot_range=loom_config['plot_range'],
        plot_two_way_streets=plot_two_way_streets,
        soliton_tree_data=soliton_tree_data,
        logger_name=get_logger_name(),
    )

#    if plot_two_way_streets is True:
#        for sn in spectral_network_data.spectral_networks:
#            if len(sn.streets) == 0:
#                continue
#            else:
#                initial_phase = '{:.3f}'.format(sn.phase / pi)
#                break
#    else:
#        initial_phase = '{:.3f}'.format(
#            spectral_network_data.spectral_networks[0].phase / pi
#        )
    initial_phase = '{:.3f}'.format(
        spectral_network_data.spectral_networks[0].phase / pi
    )

    legend = get_legend(
        g_data=sw_data.g_data,
        regular_punctures=sw_data.regular_punctures,
        branch_points=sw_data.branch_points,
        irregular_singularities=sw_data.irregular_singularities,
    )

    with open('static/bokeh_callbacks.js', 'r') as fp:
        bokeh_custom_script = fp.read()

    if len(spectral_network_data.spectral_networks) > 1:
        show_sn_slider = True
    else:
        show_sn_slider = False

    return flask.render_template(
        'plot.html',
        process_uuid=process_uuid,
        bokeh_plot_script=bokeh_plot_script,
        div=div,
        progress_log=progress_log,
        plot_legend=legend,
        bokeh_custom_script=bokeh_custom_script,
        download=str(download),
        loom_config=loom_config,
        config_options=config_options,
        advanced_config_options=advanced_config_options,
        initial_phase=initial_phase,
        n_processes=n_processes,
        saved_data=saved_data,
        default_search_radius=loom_config['size_of_bp_neighborhood'],
        plot_two_way_streets=str(plot_two_way_streets),
        search_radius=search_radius,
        show_sn_slider=str(show_sn_slider),
    )


def get_logging_file_path(logger_name):
    logging_file_path = os.path.join(
        get_loom_dir(),
        ('logs/{}_{}-{:02}-{:02} {:02}:{:02}:{:02}.log'
         .format(logger_name, *time.localtime(time.time())[:6])),
    )
    return logging_file_path


def get_cache_dir(process_uuid):
    cache_dir = os.path.join(
        get_loom_dir(),
        LOOM_CACHE_DIR,
        process_uuid,
    )
    return cache_dir


def get_full_data_dir(process_uuid, saved_data):
    if saved_data is True:
        data_dir = process_uuid
        full_data_dir = os.path.join(
            get_loom_dir(), 'data', data_dir
        )
    else:
        full_data_dir = get_cache_dir(process_uuid)

    return full_data_dir


def get_data_file_path_list(data_dir):
    data_file_path_list = glob.glob(os.path.join(data_dir, 'data_*.json'))
    data_file_path_list.sort()

    data_file_path_list += [
        os.path.join(data_dir, 'version'),
        os.path.join(data_dir, 'config.ini'),
        os.path.join(data_dir, 'sw_data.json'),
    ]

    return data_file_path_list


def record_stat(stat_logger_name, ip, uuid, data_size):
    stat_logger = logging.getLogger(stat_logger_name)
    stat_logger.info('{}, {}, {}'.format(ip, uuid, data_size))
