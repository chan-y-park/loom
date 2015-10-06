from bokeh.io import vform
from bokeh.models import CustomJS, ColumnDataSource, Slider
from bokeh.models import HoverTool, BoxZoomTool, PanTool, WheelZoomTool
from bokeh.plotting import figure, output_file, show

from misc import split_with_overlap

def get_spectral_network_bokeh_plot(
    spectral_networks, sw_data, plot_range=None,
    plot_joints=False, plot_data_points=False, plot_on_cylinder=False,
):
    # Data source for branch points.
    bpds = ColumnDataSource({'x': [], 'y': [], 'label': [],})

    for bp in sw_data.branch_points:
        bpds.data['x'].append(bp.z.real)
        bpds.data['y'].append(bp.z.imag)
        bpds.data['label'].append(bp.label)


    # Data source for the current plot
    cds = ColumnDataSource({
        'xs': [],
        'ys': [],
        'color': [],
        #'arrow_xs': [],
        #'arrow_ys': [],
        'label': [],
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
        sn_data['arrow_xs'] = []
        sn_data['arrow_ys'] = []
        sn_data['label'] = []

        for s_wall in sn.s_walls:
            z_segs = split_with_overlap(s_wall.z, s_wall.get_splittings())
            for z_seg in z_segs:
                sn_data['xs'].append(z_seg.real)
                sn_data['ys'].append(z_seg.imag)
            for root in s_wall.local_roots:
                sn_data['label'].append(str(root))
                sn_data['color'].append(sw_data.g_data.root_color(root))

        snds.data['spectral_networks'].append(sn_data)

    # Initialization of the current plot data source.
    for key in cds.data.keys():
        cds.data[key] = snds.data['spectral_networks'][0][key]

    # Setup tools.
    hover = HoverTool(
        tooltips=[('label', '@label'),]
    )

    # Start plotting data. 
    if plot_range is None:
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
        plot_x_range = (x_min, x_max)
        plot_y_range = (y_min, y_max)

    else:
        plot_x_range, plot_y_range = plot_range

    bokeh_figure = figure(
        tools=[BoxZoomTool(), PanTool(), WheelZoomTool(), hover],
        #plot_width=None,
        #plot_height=None,
        #title=None, 
        x_range=plot_x_range,
        y_range=plot_y_range,
    )

    bokeh_figure.x(
        'x', 'y', size=10, color="#e6550D", line_width=2, source=bpds
    )

    bokeh_figure.multi_line(
        xs='xs', ys='ys', color='color', source=cds
    )

    callback = CustomJS(
        args={'cds': cds, 'snds': snds}, 
        code="""
        var cd = cds.get('data');
        var snd = snds.get('data');
        var plot_idx = cb_obj.get('value');
        
        cd['xs'] = snd['spectral_networks'][plot_idx]['xs'];
        cd['ys'] = snd['spectral_networks'][plot_idx]['ys'];
        cd['label'] = snd['spectral_networks'][plot_idx]['label'];
        cd['color'] = snd['spectral_networks'][plot_idx]['color'];
        
        cds.trigger('change');
    """)

    slider = Slider(start=0, end=len(spectral_networks), 
                    value=0, step=1, title="plot index",
                    callback=callback)
    layout = vform(slider, bokeh_figure)

    return layout
