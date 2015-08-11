import numpy
from cmath import phase, pi, exp
from math import cos, sin

from bokeh.plotting import (figure, output_notebook, ColumnDataSource,
                            output_file)
from bokeh.palettes import Spectral10
from bokeh.models import HoverTool, BoxZoomTool, PanTool, WheelZoomTool, Callback, ColumnDataSource, DataTable, TableColumn
from bokeh.io import show, vform

from IPython.html.widgets import IntSlider, interactive
from IPython.display import display

from random import randint


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

def plot_s_wall(s, figure, label, root_color_map):
    n_sec = len(s.splittings) + 1

    if n_sec == 1:
        sections = [s.z]
    else:
        sections = [s.z[0:(s.splittings[0]+1)]] + \
                    [s.z[s.splittings[i]:(s.splittings[i+1]+1)] \
                        for i in range(len(s.splittings)-2)] + \
                    [s.z[s.splittings[-1]:-1]]

    roots = s.local_roots

    for i, sec in enumerate(sections):
        z_r = [z.real for z in sec]
        z_i = [z.imag for z in sec]

        # plot the lines
        sec_root =roots[i]
        color = [k for k, v in root_color_map.iteritems() \
                                if numpy.array_equal(v, sec_root)][0]
        figure.line(z_r, z_i, line_width=2, line_color=color)

    ### Plot arrows
    z_r = [z.real for z in s.z]
    z_i = [z.imag for z in s.z]
    hl = int(numpy.floor(len(z_r) / 2.0))
    hl_angle = phase((z_r[hl+1] - z_r[hl]) + 1j*(z_i[hl+1] - z_i[hl]))
    # fl = len(z_r) - 1
    # fl_angle = phase((z_r[fl] - z_r[fl-1]) + 1j*(z_i[fl] - z_i[fl-1]))
    # text_z = s.z[fl] + dynamic_displacement(fl_angle)
    figure.triangle(x=z_r[hl], y=z_i[hl], size=7, line_width=2, 
        angle=(hl_angle - pi/2), line_color='black', fill_color='black')

        
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
                text_font='times', text_font_style='italic', text_font_size='12pt'
               )

    pass


def delete_duplicates(l):
    seen = set()
    uniq = []
    for x in l:
        if x not in seen:
            uniq.append(x)
            seen.add(x)
    return uniq


def plot_branch_points(ramification_points, figure, y_max):
    loci = delete_duplicates([r.z for r in ramification_points])
    for z in loci:
        figure.cross(x=z.real, y=z.imag, size=20, color="#E6550D", line_width=2, angle=0.5)
        figure.line(x=[z.real,z.real], y=[z.imag, y_max], line_width=2, color="#E6550D", line_dash='dashed')
    return None


def plot_swn(swn, figure, y_max, root_color_map):
    # add the S-walls
    for i, s in enumerate(swn.s_walls):
        plot_s_wall(s, figure, str(i), root_color_map)

    # add the branch points
    plot_branch_points(swn.sw_data.branch_points, figure, y_max)


def plot_multi_swn(data):
    swn_list = data.spectral_networks
    g_data = swn_list[0].sw_data.g_data
    root_color_map = create_root_color_map(g_data)
    plots = []
    swn_tables = []
    #TO BE GIVEN BY HAND, CLIPPING BOUNDARY
    y_max = max([max([max([z.imag for z in s.z]) 
                for s in swn.s_walls]) 
                for swn in swn_list]) 
    for swn in swn_list:
        # create a figure object
        p = figure(width=600, height=600)
        plot_swn(swn, p, y_max, root_color_map)
        plot_tables = swn_data_tables(swn)
        plots.append(vform(p, plot_tables))
    return plots


def plots_with_slider(data):
    plots = plot_multi_swn(data)
    
    def show_plot(plot_index):
        show(plots[plot_index])

    plot_slider = interactive(show_plot, plot_index=IntSlider(min=0, max=len(plots)-1, step=1, value=0))
    return display(plot_slider)



def dummy_weight():
    return [randint(0, 1) for j in range(5)]

def dummy_weight_label():
    return 'mu_'+str(randint(0,5))


def s_wall_data_table(swn):
    g_data = swn.sw_data.g_data
    g_roots = list(g_data.roots)

    root_dictionary = {'alpha_' + str(i) : rt for i, rt in enumerate(g_roots)}

    data = dict(
            s_wall_label=[s.label for s in swn.s_walls],
            root_label=[ \
                        [\
                            [k for k, v in root_dictionary.iteritems() \
                                if numpy.array_equal(v, rt)][0] \
                            for rt in s.local_roots\
                        ] \
                        for s in swn.s_walls\
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


def weight_data_table(weight_dictionary):
    weight_labels = weight_dictionary.keys()
    weights = weight_dictionary.values()
    
    data = dict(
            weight_label=[wl for wl in weight_labels],
            weight=[str(w) for w in weights],
        )
    source = ColumnDataSource(data)

    columns = [
            TableColumn(field="weight_label", title="Weights"),
            TableColumn(field="weight", title=""),
        ]
    data_table = DataTable(source=source, columns=columns, width=800, height=300)

    return data_table

def root_data_table(root_dictionary, data):
    root_labels = root_dictionary.keys()
    roots = root_dictionary.values()
    g_data = data.spectral_networks[0].sw_data.g_data

    data = dict(
            root_label=[rl for rl in root_labels],
            root=[str(r) for r in roots],
            weight_pairs=[
                            [
                                ('(mu_'+str(p[0])+',mu_'+str(p[1])+')') \
                                for p in g_data.ordered_weight_pairs(rt)] \
                        for rt in roots]
        )
    source = ColumnDataSource(data)

    columns = [
            TableColumn(field="root_label", title="Roots"),
            TableColumn(field="root", title=""),
            TableColumn(field="weight_pairs", title="Pairs of weights"),
        ]
    data_table = DataTable(source=source, columns=columns, width=800, height=300)

    return data_table

def swn_data_tables(swn):
    return s_wall_data_table(swn)

def g_data_tables(data):
    swn_list = data.spectral_networks
    g_data = swn_list[0].sw_data.g_data
    g_roots = list(g_data.roots)
    g_weights = list(g_data.weights)
    weight_dictionary = {'mu_' + str(i) : w for i, w in enumerate(g_weights)}
    root_dictionary = {'alpha_' + str(i) : rt for i, rt in enumerate(g_roots)}
    show(vform(weight_data_table(weight_dictionary))) 
    show(vform(root_data_table(root_dictionary, data)))



def create_root_color_map(g_data):
    g_roots = list(g_data.roots)
    n_rts = len(g_roots)
    x = numpy.random.random(size=n_rts) * 100
    y = numpy.random.random(size=n_rts) * 100
    z = numpy.random.random(size=n_rts) * 100
    colors = [
    "#%02x%02x%02x" % (r, g, b) for r, g, b in \
            zip(numpy.floor(50+2*x), numpy.floor(50+2*y), numpy.floor(50+2*z))
    ]
    return {colors[i] : rt for i, rt in enumerate(g_roots)}

