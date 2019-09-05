import os
import subprocess
import time
import flask
import logging
import uuid
import zipfile
import glob
import shutil
# import pdb

from io import BytesIO
from cmath import pi

from numpy import inf

from api import (
    SpectralNetworkData,
    get_loom_dir,
#    get_logging_handler,
    set_logging,
)
from web_api import(
    WEB_APP_NAME,
    STAT_LOGGER_NAME,
    DEFAULT_NUM_PROCESSES,
    LoomDB,
    get_logging_file_path,
    get_loom_config,
    get_logger_name,
    get_full_data_dir,
    get_data_file_path_list,

)
from config import LoomConfig
from bokeh_plot import get_spectral_network_bokeh_plot
#from plotting import get_legend
from plot_api import get_sw_data_legend, SolitonTreePlot

# Flask configuration
DEBUG = True
SECRET_KEY = 'web_loom_key'

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
    ).decode('utf8').strip().split("\n")
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

    return flask.render_template(
        'load.html',
        data_names=data_names,
        n_processes=n_processes_val,
    )


def progress():
    loom_db = flask.current_app.loom_db

    loom_config = None
    prev_process_uuid = None
    new_process_uuid = None
    spectral_network_data = None
    full_data_dir = None
    event_source_url = None
    text_area_content = ''

    kwargs = {
        'n_processes': None,
        'search_radius': None,
        'saved_data': None,
        'additional_n_steps': 0,
        'new_mass_limit': None,
        'additional_iterations': 0,
        'additional_phases': None
    }

    kwargs_string_valued = {
        'task': None,
        'process_uuid': None,
        'data_name': None,
    }

    set_kwargs_from_request(kwargs, kwargs_string_valued, flask.request.form)

    # Save the previous process_uuid to load the data.
    prev_process_uuid = kwargs['process_uuid']
    data_name = kwargs['data_name']
    saved_data = kwargs['saved_data']
    task = kwargs['task']

    # Create a new process UUID for the given task.
    new_process_uuid = str(uuid.uuid4())

    if task == 'generate':
        # Generate a new spectral network.
        logger_name = get_logger_name(new_process_uuid)
        loom_config = get_loom_config(flask.request.form, logger_name)
        kwargs['saved_data'] = False
    else:
        if prev_process_uuid is not None:
            result_queue = loom_db.get_result_queue(
                prev_process_uuid, create=False,
            )
            if result_queue is not None:
                if task == 'rotate_back' or task == 'plot_two_way_streets':
                    # No need to load the data from files.
                    # Go directly to the plot page.
                    return flask.redirect(
                        flask.url_for(
                            'plot',
                            **kwargs
                        )
                    )
                elif task == 'extend':
                    spectral_network_data = result_queue.get()
        elif data_name is None:
            raise RuntimeError(
                'No data of spectral networks to load.'
            )

        if spectral_network_data is None:
            full_data_dir = get_full_data_dir(
                process_uuid=prev_process_uuid,
                data_name=data_name,
                saved_data=saved_data,
            )
        # Now either spectral_network_data is not None,
        # which is retrieved from the previous result queue,
        # or full_data_dir is not None,
        # from which spectral_network_data will be loaded.

        if task == 'extend':
            # Extend a spectral network.
            if (
                kwargs['additional_n_steps'] == 0 and
                kwargs['additional_iterations'] == 0 and
                kwargs['new_mass_limit'] is None and
                kwargs['additional_phases'] is None
            ):
                raise RuntimeError(
                    'No additional parameter for '
                    'the extension of spectral networks.'
                )
            # An extended spectral network is a new data.
            kwargs['saved_data'] = False
        elif (
            task == 'load'
            or task == 'rotate_back'
            or task == 'plot_two_way_streets'
        ):
            pass
        else:
            raise RuntimeError('Unknown task for loom: {}'.format(task))

    kwargs['process_uuid'] = new_process_uuid

    loom_db.start_loom_process(
        loom_config=loom_config,
        spectral_network_data=spectral_network_data,
        full_data_dir=full_data_dir,
        **kwargs
    )

    event_source_url = flask.url_for(
        'logging_stream', process_uuid=kwargs['process_uuid'],
    )
    text_area_content = (
        "Start loom, uuid = {}".format(kwargs['process_uuid'])
    )

    return flask.render_template(
        'progress.html',
        event_source_url=event_source_url,
        text_area_content=text_area_content,
        **kwargs
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
    task = None
    legacy = False
    spectral_network_data = None
    kwargs = {
        'saved_data': None,
        'n_processes': None,
        'search_radius': None
    }
    kwargs_string_valued = {
        'process_uuid': None,
        'data_name': None,
        'progress_log': None
    }
    if flask.request.method == 'POST':
        set_kwargs_from_request(
            kwargs, kwargs_string_valued, flask.request.form,
        )
        process_uuid = kwargs['process_uuid']
        task = flask.request.form['task']
        spectral_network_data = loom_db.get_result(process_uuid)

    elif flask.request.method == 'GET':
        set_kwargs_from_request(
            kwargs, kwargs_string_valued, flask.request.args,
        )

        # Legacy entry point.
        try:
            kwargs['data_name'] = flask.request.args['data']
            kwargs['saved_data'] = True
            legacy = True
        except KeyError:
            pass

        try:
            task = flask.request.args['task']
        except KeyError:
            pass

        try:
            process_uuid = kwargs['process_uuid']
        except KeyError:
            pass

        if legacy is True:
            kwargs['process_uuid'] = str(uuid.uuid4())
            process_uuid = kwargs['process_uuid']
            result_queue = loom_db.get_result_queue(process_uuid, create=True)
            full_data_dir = get_full_data_dir(
                data_name=kwargs['data_name'],
                saved_data=kwargs['saved_data'],
            )
            spectral_network_data = SpectralNetworkData(
                data_dir=full_data_dir,
            )

        else:
            if process_uuid is None:
                raise RuntimeError(
                    'Need process_uuid to get the data to plot.'
                )
            result_queue = loom_db.get_result_queue(
                process_uuid, create=False,
            )
            if result_queue is None:
                raise RuntimeError(
                    'There is no result queue for process {}.'
                    .format(process_uuid)
                )
            spectral_network_data = result_queue.get()

    if task == 'rotate_back':
        spectral_network_data.rotate_back()
    else:
        spectral_network_data.reset_z_rotation()

    if task == 'plot_two_way_streets':
        plot_two_way_streets = True
    else:
        plot_two_way_streets = False

    # Put back the data in the queue for a future use.
    result_queue = loom_db.get_result_queue(process_uuid, create=False)
    if result_queue is None:
        raise RuntimeError(
            'There is no result queue for process {}.'
            .format(process_uuid)
        )
    result_queue.put(spectral_network_data)

    return render_plot_template(
        spectral_network_data,
        plot_two_way_streets=plot_two_way_streets,
        **kwargs
    )


def save_data_to_server():
    if flask.request.method == 'POST':
        process_uuid = flask.request.form['process_uuid']
        data_name = flask.request.form['data_name']
        saved_data = eval(flask.request.form['saved_data'])
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
        full_data_dir = get_full_data_dir(
            process_uuid=process_uuid,
            saved_data=saved_data,
        )
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
        data_name = flask.request.form['data_name']
        saved_data = eval(flask.request.form['saved_data'])
    else:
        raise RuntimeError

    full_data_dir = get_full_data_dir(
        process_uuid=process_uuid,
        data_name=data_name,
        saved_data=saved_data,
    )
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


def download_plot(two_way_streets=False):
    loom_db = flask.current_app.loom_db

    if flask.request.method == 'POST':
        kwargs = {
            'saved_data': None,
            'n_processes': None,
            'search_radius': None,
            'plot_two_way_streets': None
        }
        kwargs_string_valued = {
            'process_uuid': None,
            'data_name': None,
            'progress_log': None
        }
        set_kwargs_from_request(
            kwargs, kwargs_string_valued, flask.request.form,
        )
    else:
        raise RuntimeError

    process_uuid = kwargs['process_uuid']
    data_name = kwargs['data_name']
    saved_data = kwargs['saved_data']
    result_queue = loom_db.get_result_queue(
        process_uuid, create=False,
    )
    if result_queue is not None:
        spectral_network_data = result_queue.get()
    else:
        full_data_dir = get_full_data_dir(
            process_uuid=process_uuid,
            data_name=data_name,
            saved_data=saved_data,
        )
        spectral_network_data = SpectralNetworkData(data_dir=full_data_dir)
    spectral_network_data.reset_z_rotation()

    data = {}
    if two_way_streets is False:
        plot_file_name = 'loom_plot_{}.html'.format(process_uuid)
        zip_file_prefix = plot_file_name
        data[plot_file_name] = render_plot_template(
            spectral_network_data,
            download=True,
            **kwargs
        )

    else:
        plot_range = eval(flask.request.form['plot_range'])
        zip_file_prefix = 'loom_streets_{}'.format(process_uuid)
        soliton_tree_data = spectral_network_data.find_two_way_streets()
        #fp = BytesIO()
        for i, trees in enumerate(soliton_tree_data):
            for j, tree in enumerate(trees):
                #fp.seek(0)
                soliton_tree_plot = SolitonTreePlot(
                    plot_range=plot_range,
                )
                # Make a plot title.
                Z = tree.Z
                title = (
                    'SN #{}, tree #{}, '.format(i, j)
                    #+ r'$Z = $'
                    + 'Z = '
                    + '({:.6}) + ({:.6})'.format(Z.real, Z.imag)
                    #+ r'$i$'
                    + 'i'
                )
                soliton_tree_plot.draw(
                    title=title,
                    sw_data=spectral_network_data.sw_data,
                    soliton_tree=soliton_tree_data[i][j],
                )
                fp = BytesIO()
                soliton_tree_plot.figure.savefig(fp, format='pdf')
                fp.seek(0)
                file_name = '{}_{}.pdf'.format(i, j)
                # XXX
                #soliton_tree_plot.figure.savefig(
                #    'trees/' + file_name,
                #    format='pdf',
                #)
                data[file_name] = fp.read()

    zip_fp = BytesIO()
    with zipfile.ZipFile(zip_fp, 'w') as zfp:
        for file_name, data_str in data.items():
            zip_info = zipfile.ZipInfo(file_name)
            zip_info.date_time = time.localtime(time.time())[:6]
            zip_info.compress_type = zipfile.ZIP_DEFLATED
            zip_info.external_attr = 0o040664 << 16
            zfp.writestr(zip_info, data_str)
    zip_fp.seek(0)

    if result_queue is not None:
        # Put back the data in the queue for a future use.
        result_queue.put(spectral_network_data)

    return flask.send_file(
        zip_fp,
        attachment_filename='{}.zip'.format(zip_file_prefix),
        as_attachment=True,
    )

def download_two_way_streets():
    return download_plot(two_way_streets=True)


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
        '/progress', 'progress', progress, methods=['POST']
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
        '/download_two_way_streets', 'download_two_way_streets',
        download_two_way_streets,
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


def render_plot_template(
    spectral_network_data,
    process_uuid=None,
    data_name=None,
    saved_data=False,
    n_processes=None,
    download=False,
    progress_log=None,
    plot_two_way_streets=False,
    search_radius=None,
):
    loom_config = spectral_network_data.config
    sw_data = spectral_network_data.sw_data
    spectral_networks = spectral_network_data.spectral_networks
    soliton_tree_data = None

    if plot_two_way_streets is True:
        soliton_tree_data = spectral_network_data.find_two_way_streets(
            search_radius=search_radius,
        )

    # Make a Bokeh plot
    plot_range=loom_config['plot_range']
    if plot_range is None:
        x_min = inf
        x_max = -inf
        y_min = inf
        y_max = -inf
        for sn in spectral_networks:
            for s_wall in sn.s_walls:
                x = s_wall.z.real
                y = s_wall.z.imag
                new_x_min = x.min()
                new_x_max = x.max()
                new_y_min = y.min()
                new_y_max = y.max()
                if new_x_min < x_min:
                    x_min = new_x_min
                if new_x_max > x_max:
                    x_max = new_x_max
                if new_y_min < y_min:
                    y_min = new_y_min
                if new_y_max > y_max:
                    y_max = new_y_max
        # Need to maintain the aspect ratio.
        range_min = min(x_min, y_min)
        range_max = max(x_max, y_max)
        plot_x_range = plot_y_range = [range_min, range_max]
        plot_range = [plot_x_range, plot_y_range]
    
    bokeh_plot_script, div = get_spectral_network_bokeh_plot(
        spectral_network_data,
        plot_range=plot_range,
        plot_two_way_streets=plot_two_way_streets,
        soliton_tree_data=soliton_tree_data,
        logger_name=get_logger_name(),
        download=download,
    )

# XXX: Uncomment the following to remove plots with no street.
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

    legend = get_sw_data_legend(sw_data)

    with open('static/bokeh_callbacks.js', 'r') as fp:
        bokeh_custom_script = fp.read()

    if len(spectral_network_data.spectral_networks) > 1:
        show_sn_slider = True
    else:
        show_sn_slider = False

    return flask.render_template(
        'plot.html',
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
        process_uuid=process_uuid,
        data_name=data_name,
        saved_data=saved_data,
        default_search_radius=loom_config['size_of_bp_neighborhood'],
        plot_two_way_streets=str(plot_two_way_streets),
        search_radius=search_radius,
        show_sn_slider=str(show_sn_slider),
        plot_range=plot_range,
    )


def set_kwargs_from_request(kwargs, kwargs_string_valued, request_dict):
    for key in kwargs.keys():
        try:
            kwargs[key] = eval(request_dict[key])
        except (KeyError, SyntaxError):
            pass

    for key in kwargs_string_valued.keys():
        try:
            value = request_dict[key]
            if value != '':
                kwargs_string_valued[key] = value
        except KeyError:
            pass

    kwargs.update(kwargs_string_valued)
