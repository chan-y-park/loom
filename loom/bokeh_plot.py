import logging
import numpy
import bokeh
import pdb

from cmath import phase, pi
from sympy import oo
from bokeh.io import vform
from bokeh.models import CustomJS, ColumnDataSource, Slider
from bokeh.models import (HoverTool, BoxZoomTool, PanTool, WheelZoomTool,
                          ResetTool, PreviewSaveTool, TapTool,)
from bokeh.models.widgets import Button
# from bokeh.models.widgets import Toggle
from bokeh.plotting import figure

from misc import get_splits_with_overlap


def get_spectral_network_bokeh_plot(
    spectral_network_data, plot_range=None,
    plot_joints=False, plot_data_points=False, plot_on_cylinder=False,
    plot_two_way_streets=False,
    soliton_tree_data=None,
    plot_width=800, plot_height=800,
    notebook=False, logger_name=None,
    marked_points=[],
    without_errors=False,
):
    #logger = logging.getLogger(logger_name)

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

    plot_width = plot_width
    plot_height = plot_height

    x_min = min([min([min([z.real for z in s_wall.z])
                      for s_wall in sn.s_walls])
                 for sn in spectral_networks])
    x_max = max([max([max([z.real for z in s_wall.z])
                      for s_wall in sn.s_walls])
                 for sn in spectral_networks])
    y_min = min([min([min([z.imag for z in s_wall.z])
                      for s_wall in sn.s_walls])
                 for sn in spectral_networks])
    y_max = max([max([max([z.imag for z in s_wall.z])
                      for s_wall in sn.s_walls])
                 for sn in spectral_networks])
    if plot_range is None:
        # Need to maintain the aspect ratio.
        range_min = min(x_min, y_min)
        range_max = max(x_max, y_max)
        plot_x_range = plot_y_range = (range_min, range_max)
    else:
        plot_x_range, plot_y_range = plot_range

    # Setup tools.
    hover = HoverTool(
        tooltips=[
            ('name', '@label'),
            ('root', '@root'),
        ]
    )

    # Prepare a bokeh Figure.
    bokeh_figure = figure(
        tools=[ResetTool(), BoxZoomTool(), PanTool(), WheelZoomTool(),
               PreviewSaveTool(), TapTool(), hover],
        plot_width=plot_width,
        plot_height=plot_height,
        title=None,
        x_range=plot_x_range,
        y_range=plot_y_range,
    )
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
    ppds = ColumnDataSource({'x': [], 'y': [], 'label': [], 'root': []})
    for pp in (sw_data.regular_punctures + sw_data.irregular_punctures):
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
    bpds = ColumnDataSource({'x': [], 'y': [], 'label': [], 'root': []})
    for bp in sw_data.branch_points:
        if bp.z == oo:
            continue
        bpds.data['x'].append(bp.z.real)
        bpds.data['y'].append(bp.z.imag)
        bpds.data['label'].append(str(bp.label))
        root_label = ''
        for root in bp.positive_roots:
            root_label += str(root.tolist()) + ', '
        bpds.data['root'].append(root_label[:-2])

    bcds = ColumnDataSource({'xs': [], 'ys': []})
    for bl in sw_data.branch_points + sw_data.irregular_singularities:
        y_r = (2j * y_max) * complex(sw_data.branch_cut_rotation)
        bcds.data['xs'].append([bl.z.real, bl.z.real + y_r.real])
        bcds.data['ys'].append([bl.z.imag, bl.z.imag + y_r.imag])

    bokeh_figure.x(
        'x', 'y', size=10, color="#e6550D", line_width=3, source=bpds,
    )
    bokeh_figure.multi_line(
        xs='xs', ys='ys', line_width=2, color='gray', line_dash='dashed',
        source=bcds,
    )

    # Data source for the current plot
    cds = ColumnDataSource({
        'xs': [],
        'ys': [],
        'ranges': [],
        'color': [],
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
                # The first data contains all the soliton trees
                # of the two-way streets in a spectral network.
                for tree in soliton_trees:
                    tree_data = get_s_wall_plot_data(
                        tree.streets, sw_data, logger_name,
                        spectral_networks[i].phase,
                    )
                    if len(data_entry) > 0:
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
        xs='xs', ys='ys', color='color', line_width=1.5, source=cds,
    )

    bokeh_figure.triangle(
        x='arrow_x', y='arrow_y', color='color', angle='arrow_angle',
        size=8, source=cds,
    )

    # 'Redraw arrows' button.
    redraw_arrows_button = Button(
        label='Redraw arrows',
        callback=CustomJS(
            args={'cds': cds, 'x_range': bokeh_figure.x_range,
                  'y_range': bokeh_figure.y_range},
            code='redraw_arrows(cds, x_range, y_range);',
        ),
    )

    # 'Show data points' button
    show_data_points_button = Button(
        label='Show data points',
    )
    show_data_points_button.callback = CustomJS(
        args={'cds': cds, 'dpds': dpds, 'hover': hover},
        code="show_data_points(cds, dpds, hover);",
    )

    # 'Hide data points' button
    hide_data_points_button = Button(
        label='Hide data points',
    )
    hide_data_points_button.callback = CustomJS(
        args={'cds': cds, 'dpds': dpds, 'hover': hover},
        code="hide_data_points(cds, dpds, hover);",
    )
    bokeh_obj = {
        'redraw_arrows_button': redraw_arrows_button,
        'show_data_points_button': show_data_points_button,
        'hide_data_points_button': hide_data_points_button,
    }

    tree_idx_ds = ColumnDataSource({'j': ['0']})
    sn_idx_ds = ColumnDataSource({'i': ['0']})
    plot_options_ds = ColumnDataSource(
        {'notebook': [notebook], 'show_trees': [plot_two_way_streets]}
    )

    if plot_two_way_streets is True:
        next_soliton_tree_button = Button(
            label='>',
        )
        next_soliton_tree_button.callback = CustomJS(
            args={
                'cds': cds, 'snds': snds, 'sn_idx_ds': sn_idx_ds,
                'tree_idx_ds': tree_idx_ds,
            },
            code=(
                'show_next_soliton_tree(cds, snds, sn_idx_ds, tree_idx_ds);'
            ),
        )
        bokeh_obj['next_soliton_tree_button'] = next_soliton_tree_button

    num_of_plots = len(snds.data['spectral_networks'])
    if num_of_plots > 1:
        # Add a slider.
        sn_slider = Slider(
            start=0, end=num_of_plots - 1,
            value=0, step=1, title="spectral network #"
        )

        custom_js_code = ''
        with open('static/bokeh_callbacks.js', 'r') as fp:
            custom_js_code += fp.read()
            custom_js_code += '\n'
        custom_js_code += (
            'sn_slider(cb_obj, cds, snds, sn_idx_ds, dpds, pds, hover, '
            'plot_options, tree_idx_ds);'
        )

        custom_js_args = {
            'cds': cds, 'snds': snds, 'sn_idx_ds': sn_idx_ds,
            'dpds': dpds, 'pds': pds, 'hover': hover,
            'plot_options': plot_options_ds, 'tree_idx_ds': tree_idx_ds
        }

        sn_slider.callback = CustomJS(args=custom_js_args, code=custom_js_code)
        plot = vform(bokeh_figure, sn_slider, width=plot_width,)
    else:
        plot = bokeh_figure

    bokeh_obj['plot'] = plot

    if notebook is True:
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
    data_dict['arrow_x'] = []
    data_dict['arrow_y'] = []
    data_dict['arrow_angle'] = []
    data_dict['label'] = []
    data_dict['root'] = []

    for s_wall in s_walls:
        z_segs = get_splits_with_overlap(s_wall.get_splits())
        for start, stop in z_segs:
            z_r = s_wall.z[start:stop].real
            z_i = s_wall.z[start:stop].imag
            a_i = int(numpy.floor(len(z_r) / 2.0))
            # TODO: Check if the arrow is within the plot range.
            a_angle = pi
            a_angle = (
                phase((z_r[a_i] - z_r[a_i - 1]) +
                      1j * (z_i[a_i] - z_i[a_i - 1]))
                - (pi / 2.0)
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
                color = sw_data.g_data.get_root_color(roots[0])
                if color is None:
                    color = '#000000'
                    logger.warning(
                        'Unknown root color for {} '
                        'of a spectral network with phase = {}.'
                        .format(s_wall.label, sn_phase)
                    )
                data_dict['color'].append(color)
        else:
            for root in s_wall.local_roots:
                data_dict['label'].append(str(s_wall.label))
                data_dict['root'].append(str(root.tolist()))
                color = sw_data.g_data.get_root_color(root)
                if color is None:
                    color = '#000000'
                    logger.warning(
                        'Unknown root color for {} '
                        'of a spectral network with phase = {}.'
                        .format(s_wall.label, sn_phase)
                    )
                data_dict['color'].append(color)

    return data_dict
