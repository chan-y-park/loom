import time
import os
import glob
import zipfile
import logging
import Tkinter as tk
import tkFileDialog
import matplotlib
import pdb

from multiprocessing import Queue

from config import LoomConfig
from trivialization import SWDataWithTrivialization
from spectral_network import SpectralNetwork
from parallel import parallel_get_spectral_network
from plotting import NetworkPlot, NetworkPlotTk
#from logutils.queue import QueueHandler#, QueueListener

#
# QueueHandler
# Copyright (C) 2010-2013 Vinay Sajip. See LICENSE.txt for details.
#
"""
This module contains classes which help you work with queues. A typical
application is when you want to log from performance-critical threads, but
where the handlers you want to use are slow (for example,
:class:`~logging.handlers.SMTPHandler`). In that case, you can create a queue,
pass it to a :class:`QueueHandler` instance and use that instance with your
loggers. Elsewhere, you can instantiate a :class:`QueueListener` with the same
queue and some slow handlers, and call :meth:`~QueueListener.start` on it.
This will start monitoring the queue on a separate thread and call all the
configured handlers *on that thread*, so that your logging thread is not held
up by the slow handlers.

Note that as well as in-process queues, you can use these classes with queues
from the :mod:`multiprocessing` module.

**N.B.** This is part of the standard library since Python 3.2, so the
version here is for use with earlier Python versions.
"""
class QueueHandler(logging.Handler):
    """
    This handler sends events to a queue. Typically, it would be used together
    with a multiprocessing Queue to centralise logging to file in one process
    (in a multi-process application), so as to avoid file write contention
    between processes.
    
    :param queue: The queue to send `LogRecords` to.
    """

    def __init__(self, queue):
        """
        Initialise an instance, using the passed queue.
        """
        logging.Handler.__init__(self)
        self.queue = queue

    def enqueue(self, record):
        """
        Enqueue a record.

        The base implementation uses :meth:`~queue.Queue.put_nowait`. You may
        want to override this method if you want to use blocking, timeouts or
        custom queue implementations.
        
        :param record: The record to enqueue.
        """
        #self.queue.put_nowait(record)
        self.queue.put(record)

    def prepare(self, record):
        """
        Prepares a record for queuing. The object returned by this method is
        enqueued.

        The base implementation formats the record to merge the message
        and arguments, and removes unpickleable items from the record
        in-place.

        You might want to override this method if you want to convert
        the record to a dict or JSON string, or send a modified copy
        of the record while leaving the original intact.
        
        :param record: The record to prepare.
        """
        # The format operation gets traceback text into record.exc_text
        # (if there's exception data), and also puts the message into
        # record.message. We can then use this to replace the original
        # msg + args, as these might be unpickleable. We also zap the
        # exc_info attribute, as it's no longer needed and, if not None,
        # will typically not be pickleable.
        self.format(record)
        record.msg = record.message
        record.args = None
        record.exc_info = None
        return record

    def emit(self, record):
        """
        Emit a record.

        Writes the LogRecord to the queue, preparing it for pickling first.
        
        :param record: The record to emit.
        """
        try:
            self.enqueue(self.prepare(record))
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)


class SpectralNetworkData:
    """
    A container class of information relevant to
    a set of spectral networks generated from a 
    single Seiberg-Witten data.
    """
    def __init__(self, sw_data, spectral_networks):
        self.sw_data = sw_data
        self.spectral_networks = spectral_networks

def get_logging_formatter(level):
    if level == logging.DEBUG:
        logging_format = '%(module)s@%(lineno)d: %(funcName)s: %(message)s'
    elif level == logging.INFO:
        logging_format = '%(process)d: %(message)s'
    else:
        logging_format = '%(message)s'
    return logging.Formatter(logging_format)

def set_logging(level, queue=None, file_name='logs/log.loom.txt'):
    #print('Setting logging level to "{}".'.format(level))

    logger = logging.getLogger()
    logger.setLevel(level)
    formatter = get_logging_formatter(level)
    # Remove other handlers.
    for handler in logger.handlers:
        logger.removeHandler(handler)

    # Create a log file.
    fh = logging.FileHandler(file_name, 'w')
    fh.setLevel(level)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    if queue is None:
        # Create a stream handler to stderr.
        logger.addHandler(
            get_logging_handler(level, logging.StreamHandler, None)
        )
    else:
        # Create a queue handler for multiprocessing.
        logger.addHandler(
            get_logging_handler(level, QueueHandler, queue)
        )

def get_logging_handler(level, handler_class, buffer_object):
    h = handler_class(buffer_object)
    h.setLevel(level)
    h.setFormatter(get_logging_formatter(level))
    return h

def generate_spectral_network(config, phase=None):
    """
    Generate one or more spectral networks according to
    the command-line options and the configuration file
    and return a list of data obtained from SpectralNetwork.get_data()
    """

    phase_range = config['phase_range']
    sw = SWDataWithTrivialization(config)

    start_time = time.time()
    logging.info('start cpu time: %s', start_time)

    if(phase is not None):
        logging.info('Generate a single spectral network at theta = {}.'
                     .format(phase))
        spectral_network = SpectralNetwork(
            phase=phase, 
        ) 

        spectral_network.grow(config, sw)

        spectral_networks = [spectral_network]

    elif(phase_range is not None):
        logging.info('Generate multiple spectral networks.')
        logging.info('phase_range = {}.'.format(phase_range))
        spectral_networks = parallel_get_spectral_network(
            sw, 
            config,
        ) 

    end_time = time.time()
    logging.info('end cpu time: %.8f', end_time)
    logging.info('elapsed cpu time: %.8f', end_time - start_time)

    return SpectralNetworkData(sw, spectral_networks)


def load_config(config_file=None):
    if config_file is None:
        toplevel = tk.Toplevel()
        file_opts = {
            'defaultextension': '.ini',
            'initialdir': os.curdir,
            'initialfile': 'config.ini',
            'parent': toplevel,
            'title': 'Select a configuration file to load.',
        }
        config_file = tkFileDialog.askopenfilename(**file_opts)
        toplevel.destroy()
        if config_file == '':
            return None

    config = LoomConfig()
    config.read(config_file)

    return config


def load_spectral_network(
    data_dir=None,
    logging_level=None,
    logging_queue=None,
    result_queue=None,
):
    if data_dir is None:
        return (None, None)
        toplevel = tk.Toplevel()
        dir_opts = {
            'initialdir': os.curdir,
            'mustexist': False,
            'parent': toplevel,
            'title': 'Select a directory that contains data files.',
        }
        data_dir = tkFileDialog.askdirectory(**dir_opts)
        toplevel.destroy()
        if data_dir == '':
            return (None, None)

    logging.info('Opening data directory "{}"...'.format(data_dir))

    config = LoomConfig()
    config.read(os.path.join(data_dir, 'config.ini'))

    sw = SWDataWithTrivialization(config)
    spectral_network_list = []

    data_file_list = glob.glob(os.path.join(data_dir, 'data_*.json'))
    data_file_list.sort()
    for data_file in data_file_list:
        logging.info('Loading {}...'.format(data_file))
        spectral_network = SpectralNetwork()
        with open(data_file, 'r') as fp:
            spectral_network.set_from_json_data(fp)
            spectral_network_list.append(spectral_network)

    data = SpectralNetworkData(sw, spectral_network_list)
    logging.info('Finished loading data from {}.'.format(data_dir))

    if logging_queue is not None:
        logging_queue.put_nowait(None) 
        
    rv = (config, data)
    if result_queue is None:
        return rv
    else:
        result_queue.put(rv)

def save_config(config, path=None):
    if path is None:
        root = tk.Tk()
        file_opts = {
            'defaultextension': '.ini',
            'initialdir': os.curdir,
            'initialfile': 'config.ini',
            'parent': root,
            'title': 'Save the current configuration to a file.',
        }
        path = tkFileDialog.asksaveasfilename(**file_opts)
        root.destroy()
        if path == '':
            return None

    logging.info('Save configuration to {}.'.format(path))
    with open(path, 'wb') as fp:
        config.parser.write(fp)

    return None


def save_spectral_network(config, spectral_network_data, data_dir=None,
                          make_zipped_file=True):
    spectral_networks = spectral_network_data.spectral_networks
    if data_dir is None:
        ### Prepare to save spectral network data to files.
        timestamp = str(int(time.time()))
        data_dir = os.path.join(
            os.curdir,
            'data',
            timestamp
        )

    logging.info('Make a directory {} to save data.'.format(data_dir))
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    ### Save configuration to a file.
    config_file_path = os.path.join(data_dir, 'config.ini')
    save_config(config, path=config_file_path)

    ### Save spectral network data.
    for i, spectral_network in enumerate(spectral_networks):
        data_file_path = os.path.join(
            data_dir,
            'data_{}.json'.format(
                str(i).zfill(len(str(len(spectral_networks)-1)))
            )
        )
        logging.info('Saving data to {}.'.format(data_file_path))
        with open(data_file_path, 'wb') as fp:
            spectral_network.save_json_data(fp,)

    if make_zipped_file is True:
        file_list = [config_file_path]
        ### Make a compressed data file.
        file_list += glob.glob(os.path.join(data_dir, 'data_*.json'))
        zipped_file_path = data_dir + '.zip'
        logging.info('Save compressed data to {}.'.format(zipped_file_path))
        with zipfile.ZipFile(zipped_file_path, 'w',
                             zipfile.ZIP_DEFLATED) as fp:
            for a_file in file_list:
                fp.write(a_file, os.path.relpath(a_file, data_dir))

    logging.info('Finished saving data to {}.'.format(data_dir))

def make_spectral_network_plot(spectral_network_data, master=None,
                               show_plot=True, plot_range=None, **kwargs):
    sw_data = spectral_network_data.sw_data
    spectral_networks = spectral_network_data.spectral_networks
    spectral_network_plot_title = 'Spectral Network'

    if matplotlib.rcParams['backend'] == 'TkAgg':
        spectral_network_plot = NetworkPlotTk(
            master=master,
            title=spectral_network_plot_title,
            plot_range=plot_range,
        )
    else:
        spectral_network_plot = NetworkPlot(
            title=spectral_network_plot_title,
            plot_range=plot_range,
        )

    for spectral_network in spectral_networks:
        logging.info('Generating the plot of a spectral network '
                     '@ theta = {}...'.format(spectral_network.phase))
        spectral_network_plot.draw(
            spectral_network, 
            sw_data.branch_points,
            punctures=sw_data.punctures,
            irregular_singularities=sw_data.irregular_singularities,
            g_data=sw_data.g_data,
            **kwargs
        )

    if show_plot is True:
        spectral_network_plot.show()

    return spectral_network_plot
