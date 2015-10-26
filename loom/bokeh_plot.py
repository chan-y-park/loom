import numpy

from cmath import phase, pi, exp
from bokeh.io import vform
from bokeh.models import CustomJS, ColumnDataSource, Slider
from bokeh.models import (HoverTool, BoxZoomTool, PanTool, WheelZoomTool,
                          ResetTool, PreviewSaveTool)
from bokeh.plotting import figure, output_file, show

from misc import split_with_overlap

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
        plot_x_range = (x_min, x_max)
        plot_y_range = (y_min, y_max)
    else:
        plot_x_range, plot_y_range = plot_range

    # Setup tools.
    hover = HoverTool(
        #tooltips=[('label', '@label'),]
        tooltips=[
            ('name', '@label'),
            ('root', '@root'),
        ]
        #tooltips='@label'
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

    # Data source for branch points & cuts.
    bpds = ColumnDataSource({'x': [], 'y': [], 'label': [], 'root': []})
    bcds = ColumnDataSource({'xs': [], 'ys': []})
    for bp in sw_data.branch_points:
        bpds.data['x'].append(bp.z.real)
        bpds.data['y'].append(bp.z.imag)
        bpds.data['label'].append(bp.label)
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
        'color': [],
        'arrow_x': [],
        'arrow_y': [],
        'arrow_angle': [],
        'label': [],
        'root': [],
    })

    # Data source containing all the spectral networks
    snds = ColumnDataSource({
        'spectral_networks': [],
    })

    for sn in spectral_networks:
        sn_data = {}
        sn_data['xs'] = []
        sn_data['ys'] = []
        sn_data['color'] = []
        sn_data['arrow_x'] = []
        sn_data['arrow_y'] = []
        sn_data['arrow_angle'] = []
        sn_data['label'] = []
        sn_data['root'] = []

        for s_wall in sn.s_walls:
            z_segs = split_with_overlap(s_wall.z, s_wall.get_splittings())
            for z_seg in z_segs:
                z_r = z_seg.real
                z_i = z_seg.imag
                a_i = int(numpy.floor(len(z_r) / 2.0))
                a_angle = (
                    phase((z_r[a_i] - z_r[a_i - 1]) + 
                          1j*(z_i[a_i] - z_i[a_i - 1]))
                    - (pi / 2.0)
                )
                sn_data['xs'].append(z_r)
                sn_data['ys'].append(z_i)
                sn_data['arrow_x'].append(z_r[a_i])
                sn_data['arrow_y'].append(z_i[a_i])
                sn_data['arrow_angle'].append(a_angle)
            for root in s_wall.local_roots:
                #sn_data['label'].append(s_wall.label + '\n' + str(root))
                sn_data['label'].append(s_wall.label)
                sn_data['root'].append(str(root.tolist()))
                sn_data['color'].append(sw_data.g_data.root_color(root))

        snds.data['spectral_networks'].append(sn_data)

    # Initialization of the current plot data source.
    for key in cds.data.keys():
        cds.data[key] = snds.data['spectral_networks'][0][key]

    bokeh_figure.multi_line(
        xs='xs', ys='ys', color='color', line_width=1.5, source=cds,
    )

    bokeh_figure.triangle(
        x='arrow_x', y='arrow_y', color='color', angle='arrow_angle',
        size=8, source=cds,
    )

    callback_code = open('loom/bokeh_slider_callback.js', 'r').read()
    callback = CustomJS(
        args={'cds': cds, 'snds': snds, 'plot_idx_ds': plot_idx_ds}, 
        code=callback_code,
    )
    slider = Slider(start=0, end=len(spectral_networks), 
                    value=0, step=1, title="plot index",
                    callback=callback)

    layout = vform(bokeh_figure, slider, width=plot_width,)

    return layout
