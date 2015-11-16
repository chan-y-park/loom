import numpy
import bokeh

from cmath import phase, pi
from sympy import oo
from bokeh.io import vform
from bokeh.models import CustomJS, ColumnDataSource, Slider
from bokeh.models import (HoverTool, BoxZoomTool, PanTool, WheelZoomTool,
                          ResetTool, PreviewSaveTool)
from bokeh.models.widgets import Button, Toggle
from bokeh.plotting import figure

from misc import get_splits_with_overlap

def get_spectral_network_bokeh_plot(
    spectral_network_data, plot_range=None,
    plot_joints=False, plot_data_points=False, plot_on_cylinder=False,
):
    spectral_networks = spectral_network_data.spectral_networks
    sw_data = spectral_network_data.sw_data
    plot_width = 800
    plot_height = plot_width

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
    plot_idx_ds = ColumnDataSource({'i': ['0']})
    bokeh_figure = figure(
        tools=[ResetTool(), BoxZoomTool(), PanTool(), WheelZoomTool(), 
               PreviewSaveTool(), hover],
        plot_width=plot_width,
        plot_height=plot_height,
        title=None, 
        x_range=plot_x_range,
        y_range=plot_y_range,
    )
    bokeh_figure.grid.grid_line_color = None

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
    bcds = ColumnDataSource({'xs': [], 'ys': []})
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
        bcds.data['xs'].append([bp.z.real, bp.z.real])
        bcds.data['ys'].append([bp.z.imag, y_max])

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

    # Data source containing all the spectral networks
    snds = ColumnDataSource({
        'spectral_networks': [],
    })

    for sn in spectral_networks:
        sn_data = {}
        sn_data['xs'] = []
        sn_data['ys'] = []
        sn_data['ranges'] = []
        sn_data['color'] = []
        sn_data['arrow_x'] = []
        sn_data['arrow_y'] = []
        sn_data['arrow_angle'] = []
        sn_data['label'] = []
        sn_data['root'] = []

        for s_wall in sn.s_walls:
            z_segs = get_splits_with_overlap(s_wall.get_splits())
            for start, stop in z_segs:
                z_r = s_wall.z[start:stop].real
                z_i = s_wall.z[start:stop].imag
                a_i = int(numpy.floor(len(z_r) / 2.0))
                a_angle = pi
                a_angle = (
                    phase((z_r[a_i] - z_r[a_i - 1]) + 
                          1j*(z_i[a_i] - z_i[a_i - 1]))
                    - (pi / 2.0)
                )
                sn_data['xs'].append(z_r)
                sn_data['ys'].append(z_i)
                sn_data['ranges'].append(
                    [z_r.min(), z_r.max(), z_i.min(), z_i.max()]
                )
                sn_data['arrow_x'].append(z_r[a_i])
                sn_data['arrow_y'].append(z_i[a_i])
                sn_data['arrow_angle'].append(a_angle)
            for root in s_wall.local_roots:
                sn_data['label'].append(str(s_wall.label))
                sn_data['root'].append(str(root.tolist()))
                sn_data['color'].append(sw_data.g_data.get_root_color(root))

        snds.data['spectral_networks'].append(sn_data)

    # Initialization of the current plot data source.
    for key in cds.data.keys():
        cds.data[key] = snds.data['spectral_networks'][0][key]

    #dpds.data['x'] = numpy.concatenate([x for x in cds.data['xs']]) 
    #dpds.data['y'] = numpy.concatenate([y for y in cds.data['ys']]) 

    bokeh_figure.scatter(x='x', y='y', alpha=0.5, source=dpds,)

    bokeh_figure.multi_line(
        xs='xs', ys='ys', color='color', line_width=1.5, source=cds,
    )

    bokeh_figure.triangle(
        x='arrow_x', y='arrow_y', color='color', angle='arrow_angle',
        size=8, source=cds,
    )

    # Callbacks after resizing the plot.
    #bokeh_figure.x_range.callback = bokeh_figure.y_range.callback

    # 'Redraw arrows' button.
    callback_code = open('loom/javascripts/bokeh_redraw_arrows_callback.js',
                         'r').read()
    redraw_arrows_callback = CustomJS(
        args={'cds': cds, 'x_range': bokeh_figure.x_range,
              'y_range': bokeh_figure.y_range,},
        code=callback_code,
    )
    redraw_arrows_button = Button(
        label='Redraw arrows',
        callback=CustomJS(
            args={'cds': cds, 'x_range': bokeh_figure.x_range,
                  'y_range': bokeh_figure.y_range,},
            code='redraw_arrows(cds, x_range, y_range);',
        ),
    )

    # 'Show data points' button
    show_data_points_button = Toggle(
        label='Show data points',
    )
    show_data_points_button.callback = CustomJS(
        args={'cds': cds, 'dpds': dpds, 'toggle': show_data_points_button,},
        code="show_data_points(cds, dpds, toggle);",
    )

    bokeh_obj = {
        'redraw_arrows_button': redraw_arrows_button,
        'show_data_points_button': show_data_points_button,
    } 

    if len(spectral_networks) > 1:
        # Add a slider.
        slider = Slider(start=0, end=len(spectral_networks)-1, 
                        value=0, step=1, title="plot index",)
        slider.callback = CustomJS(
            args={'cds': cds, 'snds': snds, 'plot_idx_ds': plot_idx_ds,
                  'dpds': dpds, 'toggle': show_data_points_button,}, 
            code="slider(cb_obj, cds, snds, plot_idx_ds, dpds, toggle);",
        )
        plot = vform(bokeh_figure, slider, width=plot_width,)
    else:
        plot = bokeh_figure

    bokeh_obj['plot'] = plot

    return bokeh.embed.components(bokeh_obj)
