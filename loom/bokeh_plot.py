import logging
import numpy
import bokeh

from cmath import phase, pi
from copy import deepcopy
from sympy import oo
from bokeh.io import vform
from bokeh.models import CustomJS, ColumnDataSource, Slider
from bokeh.models import HoverTool
from bokeh.models.widgets import Button
from bokeh.plotting import figure


def get_spectral_network_bokeh_plot(
    spectral_network_data, plot_range=None,
    plot_joints=False, plot_data_points=False, plot_on_cylinder=False,
    plot_two_way_streets=False,
    soliton_tree_data=None,
    plot_width=800, plot_height=800,
    notebook=False,
    slide=False,
    logger_name=None,
    marked_points=[],
    without_errors=False,
    download=False,
):
    # Determine if the data set corresponds to a multi-parameter
    # configuration.
    if type(spectral_network_data.sw_data) is list:
        multi_parameter = True
    else:
        multi_parameter = False

    if without_errors is True:
        spectral_networks = [
            sn for sn in spectral_network_data.spectral_networks
            if len(sn.errors) == 0
        ]
    else:
        spectral_networks = spectral_network_data.spectral_networks

    if len(spectral_networks) == 0:
        raise RuntimeError(
            'get_spectral_network_bokeh_plot(): '
            'No spectral network to plot.'
        )

    sw_data = spectral_network_data.sw_data

    plot_x_range, plot_y_range = plot_range
    y_min, y_max = plot_y_range

    # Setup tools.
    hover = HoverTool(
        tooltips=[
            ('name', '@label'),
            ('root', '@root'),
        ]
    )

    # Prepare a bokeh Figure.
    bokeh_figure = figure(
        tools='reset,box_zoom,pan,wheel_zoom,save,tap',
        plot_width=plot_width,
        plot_height=plot_height,
        title=None,
        x_range=plot_x_range,
        y_range=plot_y_range,
    )
    bokeh_figure.add_tools(hover)
    bokeh_figure.grid.grid_line_color = None

    # Data source for marked points, which are drawn for an illustration.
    mpds = ColumnDataSource(
        {'x': [], 'y': [], 'color': [], 'label': [], 'root': []}
    )
    for mp in marked_points:
        z, color = mp
        mpds.data['x'].append(z.real)
        mpds.data['y'].append(z.imag)
        mpds.data['color'].append(color)
        mpds.data['label'].append('')
        mpds.data['root'].append('')
    bokeh_figure.circle(
        x='x', y='y', size=5, color='color', source=mpds,
    )

    # Data source for punctures.
    if multi_parameter is False:
        puncts = sw_data.regular_punctures + sw_data.irregular_punctures
    else:
        puncts = sw_data[0].regular_punctures + sw_data[0].irregular_punctures
    ppds = ColumnDataSource({'x': [], 'y': [], 'label': [], 'root': []})
    for pp in puncts:
        if pp.z == oo:
            continue
        ppds.data['x'].append(pp.z.real)
        ppds.data['y'].append(pp.z.imag)
        ppds.data['label'].append(str(pp.label))
        ppds.data['root'].append('')
    bokeh_figure.circle(
        'x', 'y', size=10, color="#e6550D", fill_color=None,
        line_width=3, source=ppds,
    )

    # Data source for branch points & cuts.
    if multi_parameter is False:
        bpds = ColumnDataSource({'x': [], 'y': [], 'label': [], 'root': []})
        for bp in sw_data.branch_points:
            if bp.z == oo:
                continue
            bpds.data['x'].append(bp.z.real)
            bpds.data['y'].append(bp.z.imag)
            bpds.data['label'].append(str(bp.label))
            positive_roots = bp.positive_roots
            if len(positive_roots) > 0:
                root_label = ''
                for root in positive_roots:
                    root_label += str(root.tolist()) + ', '
                bpds.data['root'].append(root_label[:-2])
            else:
                bpds.data['root'].append('')
        bokeh_figure.x(
            'x', 'y', size=10, color="#e6550D", line_width=3, source=bpds,
        )

        bcds = ColumnDataSource({'xs': [], 'ys': []})
        try:
            branch_cut_rotation = sw_data.branch_cut_rotation
        except AttributeError:
            branch_cut_rotation = None
        if branch_cut_rotation is not None:
            for bl in sw_data.branch_points + sw_data.irregular_singularities:
                y_r = (2j * y_max) * complex(sw_data.branch_cut_rotation)
                bcds.data['xs'].append([bl.z.real, bl.z.real + y_r.real])
                bcds.data['ys'].append([bl.z.imag, bl.z.imag + y_r.imag])

            bokeh_figure.multi_line(
                xs='xs', ys='ys', line_width=2, color='gray',
                line_dash='dashed', source=bcds,
            )

    # XXX: Need to clean up copy-and-pasted codes.
    else:
        bpds = []
        bcds = []
        for swd in sw_data:
            bpds_i = ColumnDataSource(
                {'x': [], 'y': [], 'label': [], 'root': []}
            )
            for bp in swd.branch_points:
                if bp.z == oo:
                    continue
                bpds_i.data['x'].append(bp.z.real)
                bpds_i.data['y'].append(bp.z.imag)
                bpds_i.data['label'].append(str(bp.label))
                root_label = ''
                for root in bp.positive_roots:
                    root_label += str(root.tolist()) + ', '
                bpds_i.data['root'].append(root_label[:-2])
            bpds.append(bpds_i)

            bcds_i = ColumnDataSource({'xs': [], 'ys': []})
            for bl in swd.branch_points + swd.irregular_singularities:
                y_r = (2j * y_max) * complex(swd.branch_cut_rotation)
                bcds_i.data['xs'].append([bl.z.real, bl.z.real + y_r.real])
                bcds_i.data['ys'].append([bl.z.imag, bl.z.imag + y_r.imag])
            bcds.append(bcds_i)

        # In this case the branch points and cuts will be
        # drawn differently for each spectral network.
        # Each call of the slider will deal with them.

    # Data source for the current plot
    cds = ColumnDataSource({
        'xs': [],
        'ys': [],
        'ranges': [],
        'color': [],
        'alpha': [],
        'arrow_x': [],
        'arrow_y': [],
        'arrow_angle': [],
        'label': [],
        'root': [],
    })

    # Data source for plotting data points
    dpds = ColumnDataSource({
        'x': [],
        'y': [],
    })

    # Data source for phases
    pds = ColumnDataSource({
        'phase': [],
    })
    for sn in spectral_networks:
        sn_phase = '{:.3f}'.format(sn.phase / pi)
        pds.data['phase'].append(sn_phase)

    # Data source containing all the spectral networks
    snds = ColumnDataSource({
        'spectral_networks': [],
    })

    if plot_two_way_streets is True:
        # snds['spectral_networks'] is a 2-dim array,
        # where the first index chooses a spectral network
        # and the second index chooses a soliton tree
        # of the two-way streets of the spectral network.
        for i, soliton_trees in enumerate(soliton_tree_data):
            data_entry = []
            if len(soliton_trees) == 0:
                # Fill with empty data.
                empty_data = get_s_wall_plot_data(
                    [], sw_data, logger_name,
                    spectral_networks[i].phase,
                )
                data_entry.append(empty_data)
            else:
                for tree in soliton_trees:
                    tree_data = get_s_wall_plot_data(
                        tree.streets, sw_data, logger_name,
                        spectral_networks[i].phase,
                    )
                    # The first data contains all the soliton trees
                    # of the two-way streets in a spectral network.
                    if len(data_entry) == 0:
                        data_entry.append(deepcopy(tree_data))
                    else:
                        for key in tree_data.keys():
                            data_entry[0][key] += tree_data[key]
                    data_entry.append(tree_data)

            snds.data['spectral_networks'].append(data_entry)

        init_data = snds.data['spectral_networks'][0][0]
    else:
        # snds['spectral_networks'] is a 1-dim array,
        # of spectral network data.
        for sn in spectral_networks:
            skip_plotting = False
            for error in sn.errors:
                error_type, error_msg = error
                if error_type == 'Unknown':
                    skip_plotting = True
            if skip_plotting is True:
                continue

            sn_data = get_s_wall_plot_data(
                sn.s_walls, sw_data, logger_name, sn.phase,
            )
            snds.data['spectral_networks'].append(sn_data)

        init_data = snds.data['spectral_networks'][0]

    # Initialization of the current plot data source.
    for key in cds.data.keys():
        cds.data[key] = init_data[key]

    bokeh_figure.scatter(x='x', y='y', alpha=0.5, source=dpds,)

    bokeh_figure.multi_line(
        xs='xs', ys='ys',
        color='color', alpha='alpha', line_width=1.5,
        source=cds,
    )

    bokeh_figure.triangle(
        x='arrow_x', y='arrow_y', angle='arrow_angle',
        color='color', alpha='alpha', size=8,
        source=cds,
    )

    bokeh_obj = {}
    notebook_vform_elements = []

    # XXX: Where is a good place to put the following?
    custom_js_code = ''
    if notebook is True or slide is True:
        with open('static/bokeh_callbacks.js', 'r') as fp:
            custom_js_code += fp.read()
            custom_js_code += '\n'

    # Data source for plot ranges
    if download is False and notebook is False and slide is False:
        range_callback = CustomJS(
            args={
                'x_range': bokeh_figure.x_range,
                'y_range': bokeh_figure.y_range
            },
            code=(custom_js_code + 'update_plot_range(x_range, y_range);'),
        )
        bokeh_figure.x_range.callback = range_callback
        bokeh_figure.y_range.callback = range_callback

    # 'Redraw arrows' button.
    redraw_arrows_button = Button(
        label='Redraw arrows',
        callback=CustomJS(
            args={
                'cds': cds,
                'x_range': bokeh_figure.x_range,
                'y_range': bokeh_figure.y_range
            },
            code=(custom_js_code + 'redraw_arrows(cds, x_range, y_range);'),
        ),
    )
    bokeh_obj['redraw_arrows_button'] = redraw_arrows_button
    notebook_vform_elements.append(redraw_arrows_button)

    # 'Show data points' button
    show_data_points_button = Button(
        label='Show data points',
    )
    show_data_points_button.callback = CustomJS(
        args={'cds': cds, 'dpds': dpds, 'hover': hover},
        code=(custom_js_code + 'show_data_points(cds, dpds, hover);'),
    )
    bokeh_obj['show_data_points_button'] = show_data_points_button
    notebook_vform_elements.append(show_data_points_button)

    # 'Hide data points' button
    hide_data_points_button = Button(
        label='Hide data points',
    )
    hide_data_points_button.callback = CustomJS(
        args={'cds': cds, 'dpds': dpds, 'hover': hover},
        code=(custom_js_code + 'hide_data_points(cds, dpds, hover);'),
    )
    bokeh_obj['hide_data_points_button'] = hide_data_points_button
    notebook_vform_elements.append(hide_data_points_button)

    # Prev/Next soliton tree button
    tree_idx_ds = ColumnDataSource({'j': ['0']})
    sn_idx_ds = ColumnDataSource({'i': ['0']})
    plot_options_ds = ColumnDataSource(
        {'notebook': [notebook], 'show_trees': [plot_two_way_streets]}
    )

    if plot_two_way_streets is True:
        prev_soliton_tree_button = Button(
            label='<',
        )
        prev_soliton_tree_button.callback = CustomJS(
            args={
                'cds': cds, 'snds': snds, 'sn_idx_ds': sn_idx_ds,
                'tree_idx_ds': tree_idx_ds,
                'plot_options_ds': plot_options_ds,
            },
            code=(
                custom_js_code +
                'show_prev_soliton_tree(cds, snds, sn_idx_ds, tree_idx_ds, '
                'plot_options_ds);'
            ),
        )
        bokeh_obj['prev_soliton_tree_button'] = prev_soliton_tree_button
        notebook_vform_elements.append(prev_soliton_tree_button)

        next_soliton_tree_button = Button(
            label='>',
        )
        next_soliton_tree_button.callback = CustomJS(
            args={
                'cds': cds, 'snds': snds, 'sn_idx_ds': sn_idx_ds,
                'tree_idx_ds': tree_idx_ds,
                'plot_options_ds': plot_options_ds,
            },
            code=(
                custom_js_code +
                'show_next_soliton_tree(cds, snds, sn_idx_ds, tree_idx_ds, '
                'plot_options_ds);'
            ),
        )
        bokeh_obj['next_soliton_tree_button'] = next_soliton_tree_button
        notebook_vform_elements.append(next_soliton_tree_button)

    # Slider
    num_of_plots = len(snds.data['spectral_networks'])
    if num_of_plots > 1:
        if multi_parameter is False:
            sn_slider = Slider(
                start=0, end=num_of_plots - 1,
                value=0, step=1, title="spectral network #"
            )

            sn_slider.callback = CustomJS(
                args={
                    'cds': cds, 'snds': snds, 'sn_idx_ds': sn_idx_ds,
                    'dpds': dpds, 'pds': pds, 'hover': hover,
                    'plot_options': plot_options_ds,
                    'tree_idx_ds': tree_idx_ds
                },
                code=(
                    custom_js_code +
                    'sn_slider(cb_obj, cds, snds, sn_idx_ds, dpds, pds, '
                    'hover, plot_options, tree_idx_ds);'
                ),
            )
            plot = vform(bokeh_figure, sn_slider, width=plot_width,)
            notebook_vform_elements = (
                [bokeh_figure, sn_slider] + notebook_vform_elements
            )

        else:
            # TODO: implement new js routine for sn_slider when
            # there are multiple parameters.
            # Need to draw branch points and cuts for each value of the
            # parameters.
            sn_slider = Slider(
                start=0, end=num_of_plots - 1,
                value=0, step=1, title="spectral network #"
            )

            sn_slider.callback = CustomJS(
                args={
                    'cds': cds, 'snds': snds, 'sn_idx_ds': sn_idx_ds,
                    'dpds': dpds, 'pds': pds, 'hover': hover,
                    'plot_options': plot_options_ds,
                    'tree_idx_ds': tree_idx_ds
                },
                code=(
                    custom_js_code +
                    'sn_slider(cb_obj, cds, snds, sn_idx_ds, dpds, pds, '
                    'hover, plot_options, tree_idx_ds);'
                ),
            )
            plot = vform(bokeh_figure, sn_slider, width=plot_width,)
            notebook_vform_elements = (
                [bokeh_figure, sn_slider] + notebook_vform_elements
            )

    else:
        plot = bokeh_figure
        notebook_vform_elements = (
            [bokeh_figure] + notebook_vform_elements
        )

    bokeh_obj['plot'] = plot

    if notebook is True:
        # TODO: Include phase text input
        return vform(*notebook_vform_elements, width=plot_width)
    elif slide is True:
        return plot
    else:
        return bokeh.embed.components(bokeh_obj)


def get_s_wall_plot_data(s_walls, sw_data, logger_name, sn_phase):
    logger = logging.getLogger(logger_name)

    data_dict = {}
    data_dict['xs'] = []
    data_dict['ys'] = []
    data_dict['ranges'] = []
    data_dict['color'] = []
    data_dict['alpha'] = []
    data_dict['arrow_x'] = []
    data_dict['arrow_y'] = []
    data_dict['arrow_angle'] = []
    data_dict['label'] = []
    data_dict['root'] = []

    if type(sw_data) is list:
        g_data = sw_data[0].g_data
    else:
        g_data = sw_data.g_data

    for s_wall in s_walls:
        alpha = 1.0 / s_wall.get_generation()
        z_segs = s_wall.get_segments()
        for start, stop in z_segs:
            z_r = s_wall.z[start:stop].real
            z_i = s_wall.z[start:stop].imag
            a_i = int(numpy.floor(len(z_r) / 2.0))
            # TODO: Check if the arrow is within the plot range.
            a_angle = pi
            a_angle = (
                phase(
                    (z_r[a_i] - z_r[a_i - 1]) + 1j * (z_i[a_i] - z_i[a_i - 1])
                ) - (pi / 2.0)
            )
            data_dict['xs'].append(z_r)
            data_dict['ys'].append(z_i)
            data_dict['ranges'].append(
                [z_r.min(), z_r.max(), z_i.min(), z_i.max()]
            )
            data_dict['arrow_x'].append(z_r[a_i])
            data_dict['arrow_y'].append(z_i[a_i])
            data_dict['arrow_angle'].append(a_angle)
        # XXX: temporary routine to label multiple roots.
        if s_wall.multiple_local_roots is not None:
            for roots in s_wall.multiple_local_roots:
                data_dict['label'].append(str(s_wall.label))
                root_label = ''
                for root in roots:
                    root_label += str(root.tolist()) + ', '
                data_dict['root'].append(root_label[:-2])
                color = g_data.get_root_color(roots[0])
                if color is None:
                    color = '#000000'
                    logger.warning(
                        'Unknown root color for {} '
                        'of a spectral network with phase = {}.'
                        .format(s_wall.label, sn_phase)
                    )
                data_dict['color'].append(color)
                data_dict['alpha'].append(alpha)
        elif s_wall.is_trivialized():
            for root in s_wall.local_roots:
                data_dict['label'].append(str(s_wall.label))
                data_dict['root'].append(str(root.tolist()))
                color = g_data.get_root_color(root)
                if color is None:
                    color = '#000000'
                    logger.warning(
                        'Unknown root color for {} '
                        'of a spectral network with phase = {}.'
                        .format(s_wall.label, sn_phase)
                    )
                data_dict['color'].append(color)
                data_dict['alpha'].append(alpha)
        else:
            data_dict['label'].append(str(s_wall.label))
            data_dict['root'].append([])
            # (R, G, B, A)
            # data_dict['color'].append((0, 0, 255, 1,))
            data_dict['color'].append('#0000FF')
            data_dict['alpha'].append(alpha)

    return data_dict
