import numpy

from cmath import phase, pi, exp
from math import cos, sin
from bokeh.plotting import (figure, output_notebook, ColumnDataSource,
                            output_file)
from bokeh.palettes import Spectral10
from bokeh.models import (HoverTool, BoxZoomTool, PanTool, WheelZoomTool, 
                          Callback, ColumnDataSource, DataTable, TableColumn,)
from bokeh.io import show, vform
from IPython.html.widgets import IntSlider, interactive
from IPython.display import display
from random import randint

from misc import delete_duplicates

DELTA = 0.2


def dynamic_offset(angle):
    if cos(angle) >= 0.707 and 0.0 < sin(angle) <= 0.707:
        return DELTA * exp(1j * (angle + pi / 2))
    elif cos(angle) >= 0.707 and -0.707 < sin(angle) <= 0.0:
        return DELTA * exp(1j * (angle + pi / 2))
    elif 0.707 > cos(angle) >= 0.0 and 0.0 <= sin(angle):
        return DELTA * exp(1j * (angle + pi / 2))
    elif 0.0 > cos(angle) >= -0.707 and 0.0 <= sin(angle):
        return DELTA * exp(1j * (angle - pi / 2))
    elif 0.707 > cos(angle) >= 0.0 and 0.0 > sin(angle):
        return DELTA * exp(1j * (angle - pi / 2))
    elif 0.0 > cos(angle) >= -0.707 and 0.0 > sin(angle):
        return DELTA * exp(1j * (angle + pi / 2))
    elif cos(angle) < -0.707 and 0.0 < sin(angle) <= 0.707:
        return DELTA * exp(1j * (angle - pi / 2))
    elif cos(angle) < -0.707 and -0.707 < sin(angle) <= 0.0:
        return DELTA * exp(1j * (angle - pi / 2))
    

def dynamic_angle(angle):
    if cos(angle) >= 0 and sin(angle) > 0:
        return angle
    elif cos(angle) >= 0 and sin(angle) < 0:
        return angle
    elif cos(angle) < 0 and sin(angle) >= 0:
        return angle + pi
    elif cos(angle) < 0 and sin(angle) < 0:
        return angle + pi


def dynamic_alignment(angle):
    if cos(angle) >= 0.707:
        return 'left'
    elif 0.707 > cos(angle) >= 0.0:
        return 'right'
    elif 0.0 > cos(angle) >= -0.707:
        return 'left'
    elif cos(angle) < -0.707:
        return 'right'


def plot_s_wall(s, figure, label, g_data):
    splittings = s.get_splittings()
    n_sec = len(splittings) + 1

    if n_sec == 1:
        sections = [s.z]
    else:
        sections = [s.z[0:(splittings[0]+1)]] + \
                    [s.z[splittings[i]:(splittings[i+1]+1)] \
                        for i in range(len(splittings)-2)] + \
                    [s.z[splittings[-1]:-1]]

    roots = s.local_roots

    for i, sec in enumerate(sections):
        z_r = [z.real for z in sec]
        z_i = [z.imag for z in sec]

        ### Plot the lines.
        sec_root = roots[i]
        color = g_data.root_color(sec_root)
        figure.line(z_r, z_i, line_width=2, line_color=color)

    ### Plot arrows.
    z_r = [z.real for z in s.z]
    z_i = [z.imag for z in s.z]
    hl = int(numpy.floor(len(z_r) / 2.0))
    hl_angle = phase((z_r[hl+1] - z_r[hl]) + 1j*(z_i[hl+1] - z_i[hl]))
    # fl = len(z_r) - 1
    # fl_angle = phase((z_r[fl] - z_r[fl-1]) + 1j*(z_i[fl] - z_i[fl-1]))
    # text_z = s.z[fl] + dynamic_displacement(fl_angle)
    figure.triangle(
        x=z_r[hl], y=z_i[hl], size=7, line_width=2, 
        angle=(hl_angle - pi/2), line_color='black', fill_color='black',
    )

        
    ### Now plot the labels
    hl_angle = phase((z_r[hl+1] - z_r[hl]) + 1j*(z_i[hl+1] - z_i[hl]))
    text_offset = dynamic_offset(hl_angle)
    text_z = s.z[hl] + text_offset
    figure.text(x=[text_z.real], y=[text_z.imag], 
                # text=[s.label],
                text = [label],
                # angle=dynamic_angle(hl_angle),
                angle=0.0,
                text_align=dynamic_alignment(hl_angle),
                text_font='times', text_font_style='italic',
                text_font_size='12pt',)

    #pass


def plot_branch_points(bps, figure, y_max):
    loci = [bp.z for bp in bps]
    for i, z in enumerate(loci):
        figure.cross(x=z.real, y=z.imag, size=20, color="#E6550D",
                     line_width=2, angle=0.5)
        figure.line(x=[z.real,z.real], y=[z.imag, y_max], line_width=2,
                    color="#E6550D", line_dash='dashed')
        figure.text(x=[z.real], y=[z.imag], 
                    text = ['['+str(i)+']'],
                    text_align='left',
                    text_baseline='top',
                    text_font='times',
                    text_font_style='italic', 
                    text_font_size='12pt',)
    return None


#def plot_swn(swn, sw_data, figure, y_max, root_color_map):
def plot_swn(swn, sw_data, figure, y_max):
    g_data = sw_data.g_data
    ### Add S-walls
    for i, s in enumerate(swn.s_walls):
        plot_s_wall(s, figure, str(i), g_data)

    ### Add branch points
    plot_branch_points(sw_data.branch_points, figure, y_max)


def plot_multi_swn(spectral_network_data):
    swn_list = spectral_network_data.spectral_networks
    sw_data = spectral_network_data.sw_data
    #g_data = sw_data.g_data
    #root_color_map = create_root_color_map(g_data)
    plots = []
    swn_tables = []
    ### TO BE GIVEN BY HAND, CLIPPING BOUNDARY
    y_max = max([max([max([z.imag for z in s.z]) 
                for s in swn.s_walls]) 
                for swn in swn_list]) 
    for swn in swn_list:
        ### Create a figure object
        p = figure(width=600, height=600)
        #plot_swn(swn, sw_data, p, y_max, root_color_map)
        plot_swn(swn, sw_data, p, y_max)
        plot_tables = swn_data_tables(swn, sw_data)
        plots.append(vform(p, plot_tables))
    return plots


def plots_with_slider(data):
    plots = plot_multi_swn(data)
    
    def show_plot(plot_index):
        show(plots[plot_index])

    plot_slider = interactive(
        show_plot,
        plot_index=IntSlider(min=0, max=len(plots)-1, step=1, value=0)
    )
    return display(plot_slider)



def dummy_weight():
    return [randint(0, 1) for j in range(5)]


def dummy_weight_label():
    return 'mu_'+str(randint(0,5))


def s_wall_data_table(swn, sw_data):
    g_data = sw_data.g_data
    g_roots = list(g_data.roots)

    root_dictionary = {'alpha_' + str(i) : rt for i, rt in enumerate(g_roots)}

    data = dict(
        s_wall_label=[s.label for s in swn.s_walls],
        root_label=[
            [
                [k for k, v in root_dictionary.iteritems() 
                 if numpy.array_equal(v, rt)][0] for rt in s.local_roots
            ] for s in swn.s_walls
        ],
        root=[map(str, s.local_roots) for s in swn.s_walls]
    )
    source = ColumnDataSource(data)

    columns = [
        TableColumn(field="s_wall_label", title="S walls"),
        TableColumn(field="root_label", title="Root name"),
        TableColumn(field="root", title="Root type"),
    ]
    data_table = DataTable(source=source, columns=columns, width=800)

    return data_table

def branch_point_data_table(sw_data):
    ### ASSUMING square-root type branch points for now
    ### TODO : Generalize
    bpts = sw_data.branch_points
    g_data = sw_data.g_data
    g_roots = list(g_data.roots)

    root_dictionary = {'alpha_' + str(i) : rt for i, rt in enumerate(g_roots)}
    data = dict(
            bp_label=['['+str(i)+'], '+bp.label for i,bp in enumerate(bpts)],
            bp_position=[str(bp.z) for bp in bpts],
            bp_root_label=[
                [k for k, v in root_dictionary.iteritems()
                 if numpy.array_equal(v, bp.positive_roots[0])][0] 
                 for bp in bpts
            ],
            bp_root=[str(bp.positive_roots[0]) for bp in bpts]
        )
    source = ColumnDataSource(data)

    columns = [
        TableColumn(field="bp_label", title="branch points"),
        TableColumn(field="bp_position", title="position"),
        TableColumn(field="bp_root_label", title="Root label"),
        TableColumn(field="bp_root", title="Root type"),
    ]
    data_table = DataTable(source=source, columns=columns, width=800)

    return data_table


def weight_data_table(weight_dictionary):
    weight_labels = weight_dictionary.keys()
    weights = weight_dictionary.values()
    
    data = dict(weight_label=[wl for wl in weight_labels],
                weight=[str(w) for w in weights],)
    source = ColumnDataSource(data)

    columns = [
        TableColumn(field="weight_label", title="Weights"),
        TableColumn(field="weight", title=""),
    ]
    data_table = DataTable(source=source, columns=columns,
                           width=800, height=300)

    return data_table

def root_data_table(root_dictionary, sw_data):
    root_labels = root_dictionary.keys()
    roots = root_dictionary.values()
    g_data = sw_data.g_data

    data = dict(
        root_label=[rl for rl in root_labels],
        root=[str(r) for r in roots],
        weight_pairs=[
            [('(mu_'+str(p[0])+',mu_'+str(p[1])+')') 
             for p in g_data.ordered_weight_pairs(rt)] for rt in roots
        ]
    )
    source = ColumnDataSource(data)

    columns = [
        TableColumn(field="root_label", title="Roots"),
        TableColumn(field="root", title=""),
        TableColumn(field="weight_pairs", title="Pairs of weights"),
    ]
    data_table = DataTable(source=source, columns=columns,
                           width=800, height=300)

    return data_table

def swn_data_tables(swn, sw_data):
    return s_wall_data_table(swn, sw_data)

def g_data_tables(spectral_network_data):
    swn_list = spectral_network_data.spectral_networks
    sw_data = spectral_network_data.sw_data
    g_data = sw_data.g_data
    g_roots = list(g_data.roots)
    g_weights = list(g_data.weights)
    weight_dictionary = {'mu_' + str(i) : w for i, w in enumerate(g_weights)}
    root_dictionary = {'alpha_' + str(i) : rt for i, rt in enumerate(g_roots)}
    show(vform(weight_data_table(weight_dictionary))) 
    show(vform(root_data_table(root_dictionary, sw_data)))
    show(vform(branch_point_data_table(sw_data)))