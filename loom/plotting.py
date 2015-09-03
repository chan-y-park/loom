import os
import numpy
import pdb
import logging
import Tkinter as tk
#import mpldatacursor

import matplotlib
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg as FigureCanvas,
    NavigationToolbar2TkAgg as NavigationToolbar,
)
from matplotlib.widgets import Button
from matplotlib import pyplot
from math import pi

from network_plot import NetworkPlotBase
#from misc import PSL2C, put_on_cylinder, split_with_overlap
from misc import put_on_cylinder, split_with_overlap
#from bokeh_plot import create_root_color_map, root_color

class SpectralNetworkPlotBase(NetworkPlotBase):
    def draw(
        self,
        spectral_network,
        branch_points,
        punctures=None, 
        plot_range=None, 
        plot_joints=False,
        plot_data_points=False,
        plot_on_cylinder=False,
        C=[[1, 0], [0, 1]],
        g_data=None
    ):
        
        labels = {'branch_points': [], 'joints': [], 'walls': []}
        if plot_range is None:
            if plot_on_cylinder is True:
                plot_range = [[-pi, pi], [-5, 5]]

        branch_points_z = []
        for i, bp in enumerate(branch_points):
            if plot_on_cylinder is True:
                bp_z = put_on_cylinder(bp.z, C)
            else:
                bp_z = bp.z
            branch_points_z.append([bp_z.real, bp_z.imag])
            labels['branch_points'].append(bp.label)
   
        joints = []
        for i, jp in enumerate(spectral_network.joints):
            if plot_on_cylinder is True:
                jp_z = put_on_cylinder(jp.z, C)
            else:
                jp_z = jp.z
            joints.append([jp_z.real, jp_z.imag])
            labels['joints'].append(jp.label)

        walls = []
        walls_roots = []
        for i, s_wall in enumerate(spectral_network.s_walls):
            segments = []
            seg_labels = []
            split_at = []

            if plot_on_cylinder is True:
                zs_on_cylinder = numpy.fromfunction(
                    lambda i: put_on_cylinder(s_wall.z[i], C),
                    (len(s_wall.z)),
                )
                for j, delta_z in enumerate(numpy.diff(zs_on_cylinder)):
                    if abs(delta_z) > pi:
                        split_at.append(j)
                z_segs = numpy.split(zs_on_cylinder, split_at)
            else:
                z_segs = split_with_overlap(s_wall.z, s_wall.get_splittings())
                
            seg_labels = [s_wall.label + '\n' + lab
                          for lab in map(str, s_wall.local_roots)]

            for z_seg in z_segs:
                segments.append([z_seg.real, z_seg.imag])
                # seg_labels.append(s_wall.label)

            walls.append(segments)
            walls_roots.append(s_wall.local_roots)
            #root_color_map = create_root_color_map(g_data)
            walls_colors = [
                [g_data.root_color(root) for root in w_roots]
                for w_roots in walls_roots
            ]
            labels['walls'].append(seg_labels)


        logging.info(
                '\n------------------------'+
                '\nphase : {}'.format(spectral_network.phase) +
                '\n------------------------\n'
            )
        print_legend(g_data)

        print_spectral_network_data(
                spectral_network.s_walls, 
                branch_points,
                g_data
            )

        super(SpectralNetworkPlotBase, self).draw(
            phase=spectral_network.phase,
            branch_points=branch_points_z,
            joints=joints,
            walls=walls,
            walls_colors=walls_colors,
            labels=labels,
            plot_range=plot_range,
            plot_joints=plot_joints,
            plot_data_points=plot_data_points,
        )


class NetworkPlot(SpectralNetworkPlotBase):
    """
    This class implements UIs using matplotlib widgets
    so that it can be backend-independent.

    The content of this class is independent of the parent class.
    It only depends on the grandparent class, which should be
    'NetworkPlotBase'. Therefore this class can inherit any class
    whose parent is 'NetworkPlotBase'; just change the name of the
    parent in the definition of this class.
    """
    def __init__(
        self,
        title=None,
    ):
        super(NetworkPlot, self).__init__(
            matplotlib_figure=pyplot.figure(title),
        )

        self.axes_button_prev = None
        self.axes_button_next = None
        self.index_text = None

    def save(self, plot_dir, file_prefix=''):
        # TODO: change the current figure to plot_id.
        digits = len(str(len(self.plots)-1))
        for i, axes in enumerate(self.plots):
            self.change_current_plot(i)
            plot_file_path = os.path.join(
                plot_dir, file_prefix + str(i).zfill(digits) + '.png'
            )
            self.figure.savefig(plot_file_path)


    def change_current_plot(self, new_plot_idx):
        super(NetworkPlot, self).change_current_plot(new_plot_idx)
        if self.index_text is not None:
            self.index_text.set_text(
                "{}/{}".format(self.current_plot_idx, len(self.plots)-1)
            )



    def show_prev_plot(self, event):
        super(NetworkPlot, self).show_prev_plot(event)


    def show_next_plot(self, event):
        super(NetworkPlot, self).show_next_plot(event)


    def show(self):
        plot_idx = 0
        self.current_plot_idx = plot_idx
        self.plots[plot_idx].set_visible(True)
        self.set_data_cursor()

        if(len(self.plots) > 1):
            button_width = .05
            index_width = .03*len(str(len(self.plots)-1))
            button_height = .05
            button_bottom = .025
            center = .5
            margin = .005

            axes_prev_rect = [center - index_width/2 - margin - button_width,
                              button_bottom, button_width, button_height]
            axes_prev = self.figure.add_axes(axes_prev_rect)
            self.button_prev = Button(axes_prev, '<')
            self.button_prev.on_clicked(self.show_prev_plot)

            self.index_text = self.figure.text(
                center - index_width/2, (button_bottom+button_height)/2,
                "{}/{}".format(self.current_plot_idx, len(self.plots)-1)
            )

            axes_next_rect = [center + index_width/2 + margin,
                              button_bottom, button_width, button_height]
            axes_next = self.figure.add_axes(axes_next_rect)
            self.button_next = Button(axes_next, '>')
            self.button_next.on_clicked(self.show_next_plot)

        self.figure.show()


class NetworkPlotTk(SpectralNetworkPlotBase):
    """
    This class implements UIs using Tkinter.

    The content of this class is independent of the parent class.
    It only depends on the grandparent class, which should be
    'NetworkPlotBase'. Therefore this class can inherit any class
    whose parent is 'NetworkPlotBase'; just change the name of the
    parent in the definition of this class.
    """
    def __init__(self,
        master=None,
        title=None,
    ):
        super(NetworkPlotTk, self).__init__(
            #matplotlib_figure=pyplot.figure(title),
            matplotlib_figure=matplotlib.figure.Figure(),
        )

        #self.root = None
        if master is None:
            root = tk.Tk()
            root.withdraw()
            self.root = root
        else:
            self.root = master
        self.master = master

        # Create a Toplevel widget, which is a child of GUILoom
        # and contains plots,
        self.toplevel = tk.Toplevel(self.root)
        self.toplevel.wm_title(title)
        self.toplevel.protocol("WM_DELETE_WINDOW", self.on_closing)

        self.plot_idx_scale = None

        self.plot_idx_entry = None
        self.plot_idx_entry_var = tk.StringVar()
        self.plot_idx_entry_var.trace('w', self.plot_idx_entry_change)

        self.canvas = FigureCanvas(
            self.figure,
            master=self.toplevel,
            resize_callback=self.canvas_resize_callback
        )
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar(self.canvas, self.toplevel)
        toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
   


    def on_closing(self):
        self.toplevel.destroy()
        if self.master is None:
            self.root.destroy()


    def scale_action(self, scale_value):
        new_plot_idx = int(scale_value)
        self.update_current_plot(new_plot_idx)
        self.plot_idx_entry_var.set(new_plot_idx)


    def plot_idx_entry_change(self, *args):
        try:
            new_plot_idx = int(self.plot_idx_entry_var.get())

            if new_plot_idx == self.current_plot_idx:
                return None
            elif new_plot_idx < 0:
                new_plot_idx = 0
            elif new_plot_idx > len(self.plots) - 1:
                new_plot_idx = len(self.plots) - 1

            self.plot_idx_scale.set(new_plot_idx)
            self.update_current_plot(new_plot_idx)

        except ValueError:
            pass

        return None

    def update_current_plot(self, new_plot_idx):
        if self.data_cursor is not None:
            self.data_cursor.hide().disable()

        self.plots[self.current_plot_idx].set_visible(False)
        self.plots[new_plot_idx].set_visible(True)
        # Update the index variable for the currently displayed plot.
        self.current_plot_idx = new_plot_idx
        self.set_data_cursor()
        self.canvas.draw_idle()
        self.canvas.get_tk_widget().focus_set()

        return None


    def canvas_resize_callback(self, event):
        self.set_data_cursor()


    def show(self):
        plot_idx = 0
        self.current_plot_idx = plot_idx
        self.plots[plot_idx].set_visible(True)
        self.set_data_cursor()

        if(len(self.plots) > 1):
            tk.Label(
                self.toplevel,
                text='Plot #',
            ).pack(side=tk.LEFT)

            self.plot_idx_entry_var.set(plot_idx)
            self.plot_idx_entry = tk.Entry(
                self.toplevel,
                textvariable=self.plot_idx_entry_var,
                width=len(str(len(self.plots)-1)),
            )
            self.plot_idx_entry.pack(side=tk.LEFT)

            tk.Label(
                self.toplevel,
                text='/{}'.format(len(self.plots)-1),
            ).pack(side=tk.LEFT)

            self.plot_idx_scale = tk.Scale(
                self.toplevel,
                command=self.scale_action,
                #length=100*len(self.plots),
                orient=tk.HORIZONTAL,
                showvalue=0,
                to=len(self.plots)-1,
                variable=self.current_plot_idx,
            )
            self.plot_idx_scale.pack(
                expand=True,
                fill=tk.X,
                side=tk.LEFT,
            )


def make_root_dictionary(g_data):
    g_roots = list(g_data.roots)
    root_dictionary = (
                {'alpha_' + str(i) : rt for i, rt in enumerate(g_roots)}
            )
    return root_dictionary

def make_weight_dictionary(g_data):
    g_weights = list(g_data.weights)
    weight_dictionary = (
                {'mu_' + str(i) : w for i, w in enumerate(g_weights)}
            )
    return weight_dictionary


def print_spectral_network_data(s_walls, branch_points, g_data):
    root_dictionary = make_root_dictionary(g_data)
    #weight_dictionary = make_weight_dictionary(g_data)

    logging.info('\n\t--- The S-Wall Data ---\n')
    for s in s_walls:
        rt_labels = [get_label(rt, root_dictionary) for rt in s.local_roots]
        wt_labels = [
                [
                    [
                    'mu_' + str(pair[0]), 
                    'mu_' + str(pair[1])
                    ]
                    for pair in loc_wts
                ]
                for loc_wts in s.local_weight_pairs
            ]
        logging.info(
                '\n' + s.label + 
                '\troot types : {}\n'.format(rt_labels) +
                '\t\tsheet pairs : {}\n'.format(wt_labels)
            )

    logging.info('\n\t--- The Branch Points ---\n')
    for bp in branch_points:
        rt_labels = [get_label(rt, root_dictionary) for rt in bp.positive_roots]
        logging.info(
                '\n' + bp.label + 
                '\tposition : {}\n'.format(bp.z) +
                '\t\troot type : {}\n'.format(rt_labels)
            )


def get_label(value, dictionary):
    return [k for k, v in dictionary.iteritems() 
            if numpy.array_equal(v, value)][0]


def print_legend(g_data):
    root_dictionary = make_root_dictionary(g_data)
    weight_dictionary = make_weight_dictionary(g_data)
    root_labels = root_dictionary.keys()
    roots = root_dictionary.values()
    weight_pairs=[
                [
                    str('(mu_'+str(p[0])+', mu_'+str(p[1])+')') 
                    for p in g_data.ordered_weight_pairs(rt)
                ] for rt in roots
            ]
    weight_labels = weight_dictionary.keys()
    weights = weight_dictionary.values()

    logging.info('\n\t--- The Root System ---\n')
    for i in range(len(roots)):
        logging.info(
                root_labels[i] + 
                '\t\t{}\n'.format(list(roots[i])) +
                'ordered weight pairs : \t{}\n'.format(weight_pairs[i])
            )

    logging.info('\n\t--- The Weight System ---\n')
    for i in range(len(weights)):
        logging.info(
                weight_labels[i] + 
                '\t\t{}\n'.format(list(weights[i])) 
            )
    
# def plot_s_walls(
#     s_walls,
#     ramification_points=[],
#     joints=[],
#     plot_range=[[-5, 5], [-5, 5]],
#     plot_data_points=False,
#     marked_points=[],
#     colors=['b', 'g', 'r', 'c', 'm', 'y'], 
# ):
#     [x_min, x_max], [y_min, y_max] = plot_range

#     pyplot.figure(1)
#     pyplot.title('S-walls')

#     # z-plane
#     zax = pyplot.subplot(
#         121,
#         label='z-plane',
#         xlim=(x_min, x_max),
#         ylim=(y_min, y_max),
#         aspect='equal',
#     )
#     # Plot branch points
#     for rp in ramification_points:
#         zax.plot(rp.z.real, rp.z.imag, 'x', markeredgewidth=2, markersize=8, 
#                  color='k', label=rp.label,)
#     for jp in joints:
#         zax.plot(jp.z.real, jp.z.imag, '+', markeredgewidth=2, markersize=8, 
#                  color='k', label=jp.label,)
#     for p in marked_points:
#         zax.plot(p[0].real, p[0].imag, 'o', markeredgewidth=2,
#                  markersize=4, color='k')

#     # x-plane
#     xax = pyplot.subplot(
#         122,
#         label='x-plane',
#         xlim=(x_min, x_max),
#         ylim=(y_min, y_max),
#         aspect='equal',
#     )
#     for rp in ramification_points:
#         xax.plot(rp.x.real, rp.x.imag, 'x', markeredgewidth=2, markersize=8, 
#                  color='k', label=rp.label,)
#     for jp in joints:
#         for j, x_j in enumerate(jp.x):
#             xax.plot(
#                 x_j.real, x_j.imag,
#                 '+', markeredgewidth=2, markersize=8, color='k', 
#                 label=(jp.label + ', {}'.format(j)),
#             )
#     for p in marked_points:
#         xax.plot(p[1].real, p[1].imag, 'o', markeredgewidth=2,
#                  markersize=4, color='k')

#     for i, s_wall in enumerate(s_walls):
#         s_wall_color = colors[i % len(colors)]
#         # z-plane
#         zrs = s_wall.z.real
#         zis = s_wall.z.imag
#         zax.plot(zrs, zis, '-', color=s_wall_color, 
#                  label=s_wall.label)
#         if(plot_data_points == True):
#             zax.plot(zrs, zis, 'o', color=s_wall_color, label=s_wall.label)

#         # x-plane
#         xs = s_wall.x.T
#         for j, x_j in enumerate(xs):
#             xrs = x_j.real
#             xis = x_j.imag
#             xax.plot(
#                 xrs, xis, '-', color=s_wall_color, 
#                 label=(s_wall.label + ',{}'.format(j))
#             )
#             if(plot_data_points == True):
#                 xax.plot(
#                     xrs, xis, 'o', color=s_wall_color, 
#                     label=(s_wall.label + ',{}'.format(j))
#                 )
            
#     mpldatacursor.datacursor(
#         formatter='{label}'.format,
#         #tolerance=2,
#         #hover=True,
#         display='multiple',
#     )

#     pyplot.show()

# def plot_segments(
#     segments, 
#     marked_points=[],
#     plot_range=[-5, 5, -5, 5],
#     plot_data_points=False
# ):
#     """
#     Plot the given segments for debugging.
#     """
#     [x_min, x_max], [y_min, y_max] = plot_range

#     # Plot setting.
#     pyplot.xlim(x_min, x_max)
#     pyplot.ylim(y_min, y_max)
#     pyplot.axes().set_aspect('equal')

#     for segment in segments:
#         xs, ys = segment
#         pyplot.plot(xs, ys, '-')
#         if(plot_data_points == True):
#             pyplot.plot(xs, ys, 'o', color='b')

#     for p in marked_points:
#         pyplot.plot(p[0], p[1], 'x', markeredgewidth=2, markersize=8,
#                     color='k')

#     pyplot.show()


# def plot_s_walls(s_walls, plot_range=None, plot_data_points=False,):
#     """
#     Plot S-walls for debugging purpose.
#     """
#     pyplot.figure()
#     pyplot.axes().set_aspect('equal')

#     for s_wall in s_walls:
#         xs = s_wall.z.real 
#         ys = s_wall.z.imag 
#         pyplot.plot(xs, ys, '-', label=s_wall.label)

#         if(plot_data_points == True):
#             pyplot.plot(xs, ys, 'o', color='k', markersize=4)

#     if plot_range is None:
#         pyplot.autoscale(enable=True, axis='both', tight=None)
#     else:
#         [[x_min, x_max], [y_min, y_max]] = plot_range
#         pyplot.xlim(x_min, x_max)
#         pyplot.ylim(y_min, y_max)

#     mpldatacursor.datacursor(
#         formatter='{label}'.format,
#         hover=True,
#     )

#     pyplot.show()



