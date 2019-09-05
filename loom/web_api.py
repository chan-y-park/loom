import multiprocessing
#import signal
import threading
import logging
import time
import os
import sys
import glob
# import pdb

import flask

from io import StringIO

from queue import Empty as QueueEmpty

from api import (
    SpectralNetworkData,
    get_loom_dir,
    get_logging_handler,
    set_logging,
)
from config import LoomConfig
from misc import (
    get_phases_from_dict,
    get_phase_dict,
    get_data_size_of,
)

# Global variables
PARENT_LOGGER_NAME = 'loom'
WEB_APP_NAME = 'web_loom'
STAT_LOGGER_NAME = 'stat_of_' + WEB_APP_NAME
DEFAULT_NUM_PROCESSES = 4
#LOOM_PROCESS_JOIN_TIMEOUT_SECS = 3
LOOM_CACHE_DIR = 'cache'
DB_CLEANUP_CYCLE_SECS = 3600
RESULT_QUEUE_LIFETIME_SECS = (24 * 3600)
# Set RESULT_QUEUES_MAXSIZE > 0 to set the maximum number of queues.
#RESULT_QUEUES_MAXSIZE = 0
RESULT_QUEUES_MAXSIZE = 10


class LoomDBQueue(object):
    def __init__(self):
        self.queue = multiprocessing.Queue()
        self.timestamp = None

    def put(self, obj, block=True, timeout=None):
        self.queue.put(obj, block, timeout)
        self.timestamp = time.time()

    def put_nowait(self, obj):
        self.queue.put(obj, block=False)

    def get(self):
        return self.queue.get()

    def get_nowait(self):
        return self.queue.get_nowait()

    def empty(self):
        return self.queue.empty()


# TODO: kill an orphaned process gracefully.
class LoomDB(object):
    """
    The internal DB to manage loom processes.
    """
    def __init__(self, logging_level=None):
        self.logging_queues = {}
        self.result_queues = {}
        self.loom_processes = {}
        self.timestamps = {}
        self.result_queues_maxsize = RESULT_QUEUES_MAXSIZE
        self.logging_level = logging_level

#        signal.signal(signal.SIGINT, self.loom_db_stop_signal_handler)

        self.db_manager = threading.Thread(
            target=self.db_manager,
        )
        self.db_manager_stop = threading.Event()
        self.db_manager_stop.clear()
        self.db_manager.start()

#    def loom_db_stop_signal_handler(self, signum, frame):
#        logger_name = get_logger_name()
#        logger = logging.getLogger(logger_name)
#
#        if signum == signal.SIGINT:
#            msg = 'LoomDB caught SIGINT; raises KeyboardInterrrupt.'
#            e = KeyboardInterrupt
#
#        logger.warning(msg)
#        self.db_manager_stop.set()
#        raise e

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
        result_queues = self.result_queues
        while not self.db_manager_stop.wait(DB_CLEANUP_CYCLE_SECS):
            to_delete = []
            for process_uuid, result_queue in result_queues.items():
                timestamp = result_queue.timestamp
                if timestamp is None:
                    continue
                age = time.time() - timestamp
                if age > RESULT_QUEUE_LIFETIME_SECS:
                    logger.info(
                        'Deleting an old queue {}...'
                        .format(process_uuid)
                    )
                    to_delete.append(process_uuid)
            for process_uuid in to_delete:
                self.delete_queue(result_queues, process_uuid)

        # Received a stop event; finish all the processes.
        process_uuids = self.loom_processes.keys()
        for process_uuid in process_uuids:
            logger.info('Finishing process {}...'.format(process_uuid))
            self.finish_loom_process(process_uuid, join_timeout=0)
            self.delete_queue(self.result_queues, process_uuid)

        logger.info('LoomDB manager thread is finished.')

    def delete_queue(self, queues, process_uuid):
        logger = logging.getLogger(get_logger_name())
        try:
            # Flush the queue.
            queue = queues[process_uuid]
            while queue.empty() is False:
                try:
                    queue.get_nowait()
                except QueueEmpty:
                    break
            # Remove the result queue.
            del queues[process_uuid]
        except KeyError:
            logger.warning(
                "Removing queue {} failed: no such a queue exists."
                .format(process_uuid)
            )

    def get_result_queue(self, process_uuid, create=True):
        logger_name = get_logger_name()
        logger = logging.getLogger(logger_name)
        result_queues = self.result_queues

        # NOTE: To keep the size of result_queues precisely,
        # need to wrap the following with a lock.
        # Here we don't use a lock because result_queues_maxsize
        # is just a guideline.
        try:
            # Try recycling the previous result queue.
            result_queue = result_queues[process_uuid]
        except KeyError:
            if create is False:
                return None

            # If there are more queues than the maximum size,
            # remove the oldest queue.
            maxsize = self.result_queues_maxsize
            if (maxsize > 0 and len(result_queues) >= maxsize):
                oldest_process_uuid = min(
                    result_queues.keys(), key=result_queues.get
                )
                logger.info(
                    'Deleting the oldest result_queue {}...'
                    .format(oldest_process_uuid)
                )
                self.delete_queue(self.result_queues, oldest_process_uuid)

            # Create a new result queue.
            result_queue = LoomDBQueue()
            result_queues[process_uuid] = result_queue

        return result_queue

    def start_loom_process(
        self,
        loom_config=None,
        spectral_network_data=None,
        full_data_dir=None,
        process_uuid=None,
        n_processes=None,
        task=None,
        saved_data=None,
        data_name=None,
        rotate_back=None,
        plot_two_way_streets=None,
        search_radius=None,
        additional_n_steps=0, new_mass_limit=None,
        additional_iterations=0, additional_phases=None,
    ):
        logging_level = self.logging_level
        if task == 'generate' or task == 'extend':
            # Prepare a cache directory for new data.
            cache_dir = get_cache_dir(process_uuid)
            if os.path.exists(cache_dir) is False:
                os.makedirs(cache_dir)
            logging_file_name = os.path.join(cache_dir, 'log')
        elif (
            task == 'load' or
            task == 'rotate_back' or
            task == 'plot_two_way_streets'
        ):
            # Do not create a logging file
            logging_file_name = None
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
        logger = logging.getLogger(logger_name)

        result_queue = self.get_result_queue(process_uuid, create=True)

        if task == 'rotate_back':
            logger.info('Loading spectral networks to rotate back...')
        elif task == 'plot_two_way_streets':
            logger.info('Loading spectral networks to find two-way streets...')

        if (
            task == 'load' or
            task == 'rotate_back' or
            task == 'plot_two_way_streets'
        ):
            spectral_network_data = SpectralNetworkData(
                logger_name=logger_name,
            )
            loom_process = multiprocessing.Process(
                target=spectral_network_data.load,
                kwargs=dict(
                    data_dir=full_data_dir,
                    result_queue=result_queue,
                    logging_queue=logging_queue,
                ),
            )

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
            if spectral_network_data is None:
                spectral_network_data = SpectralNetworkData(
                    data_dir=full_data_dir,
                    logger_name=logger_name,
                )
            else:
                spectral_network_data.logger_name = logger_name
                for sn in spectral_network_data.spectral_networks:
                    sn.logger_name = logger_name
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
            return

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
                return
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
                return

        # Recevied all the logs, finish the SSE stream.
        yield 'event: finish\ndata: \n\n'
        return

    def get_result(self, process_uuid):
        logger_name = get_logger_name(process_uuid)
        logger = logging.getLogger(logger_name)

        #result_queue = self.result_queues[process_uuid]
        result_queue = self.get_result_queue(process_uuid, create=False)
        if result_queue is None:
            logger.error(
                'No result queue for Process {}.'
                .format(process_uuid)
            )
            return None
        loom_process = self.loom_processes[process_uuid]

        if result_queue.empty() is True:
            if loom_process.is_alive():
                logger.warning('Process {} still alive.'.format(process_uuid))
                return None
            else:
                logger.error(
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
        #join_timeout=LOOM_PROCESS_JOIN_TIMEOUT_SECS
        join_timeout=0,
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

        self.delete_queue(self.logging_queues, process_uuid)


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


def get_full_data_dir(
    process_uuid=None,
    data_name=None,
    saved_data=None,
):
    if saved_data is True:
        if data_name is None:
            raise RuntimeError
        full_data_dir = os.path.join(
            get_loom_dir(), 'data', data_name
        )
    else:
        if process_uuid is None:
            raise RuntimeError
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
