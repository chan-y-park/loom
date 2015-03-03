import numpy
import pdb
import logging
import Tkinter as tk
import mpldatacursor

import matplotlib
# use() directive must be called before importing matplotlib.pyplot
matplotlib.use('TkAgg')     

from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg as FigureCanvas,
    NavigationToolbar2TkAgg as NavigationToolbar,
)
from matplotlib import pyplot
from math import pi
from misc import PSL2C, put_on_cylinder

class SpectralNetworkPlot:
    def __init__(self, 
        master=None,
        config=None,
        plot_on_cylinder=False,
        plot_bins=False, 
        plot_joints=False,
        plot_data_points=False,
        plot_segments=False,
    ):
        self.config = config
        self.plot_on_cylinder = plot_on_cylinder
        self.plot_bins = plot_bins
        self.plot_joints = plot_joints
        self.plot_data_points = plot_data_points
        self.plot_segments = plot_segments

        # Create a Toplevel widget, which is a child of GUILoom 
        # and contains plots,
        self.toplevel = tk.Toplevel(master)
        self.toplevel.wm_title('Spectral Network Plot')

        self.plots = []
        self.data_cursor = None
        self.current_plot_idx = None 

        self.plot_idx_scale = None

        self.plot_idx_entry = None
        self.plot_idx_entry_var = tk.StringVar() 
        self.plot_idx_entry_var.trace('w', self.plot_idx_entry_change)

        self.figure = matplotlib.figure.Figure()
        self.canvas = FigureCanvas(
            self.figure,
            master=self.toplevel,
            resize_callback=(
                lambda event: self.set_data_cursor()
            )
        )
        self.canvas.show()
        self.canvas.get_tk_widget().pack()

        toolbar = NavigationToolbar(self.canvas, self.toplevel)
        toolbar.update()
        self.canvas.get_tk_widget().pack()

    
    def draw(self, spectral_network,):
        theta = spectral_network.phase
        ramification_points = spectral_network.ramification_points
        s_walls = spectral_network.s_walls

        C = self.config['mt_params']

        z_range_limits = self.config['z_range_limits']
        if z_range_limits is None:
            if self.plot_on_cylinder is True:
                z_range_limits = [-pi, pi, -5, 5]
            else:
                z_range_limits = [-5, 5, -5, 5] 
        x_min, x_max, y_min, y_max = z_range_limits
        rect = [0.125, 0.15, .8, 0.75]

        axes = self.figure.add_axes(
            rect,
            label=theta,
            xlim=(x_min, x_max),
            ylim=(y_min, y_max),
            aspect='equal',
        )

        # Draw a lattice of bins for visualization.
        if(self.plot_bins is True):
            bin_size = config['size_of_bin']
            xv = x_min
            while xv < x_max:
                xv += bin_size
                axes.axvline(x = xv, linewidth=0.5, color='0.75')

            yh = y_min  
            while yh < y_max:
                yh += bin_size
                axes.axhline(y = yh, linewidth=0.5, color='0.75')
        # End of drawing the bin lattice.

        # Plot branch points
        for rp in ramification_points:
            if (self.plot_on_cylinder is True):
                rp_z = put_on_cylinder(rp.z, C)
            else:
                rp_z = rp.z
            bpx = rp_z.real 
            bpy = rp_z.imag 
            axes.plot(bpx, bpy, 'x', markeredgewidth=2, markersize=8, 
                      color='k', label=rp.label,)
        # End of plotting branch points
   
        # Plot joints
        if(self.plot_joints is True):
            for jp in spectral_network.joints:
                if (self.plot_on_cylinder is True):
                    jp_z = put_on_cylinder(rp.z, C)
                else:
                    jp_z = jp.z
                jpx = jp_z.real 
                jpy = jp_z.imag 
                axes.plot(jpx, jpy, '+', markeredgewidth=2, markersize=8, 
                          color='k', label=jp.label,)
        # End of plotting joints

        # If we have segments of curves, draw them in different colors.
        if(self.plot_segments is True):
            hit_table = spectral_network.hit_table
            for bin_key in hit_table:
                for curve_index in hit_table[bin_key]:
                    for t_i, t_f in hit_table[bin_key][curve_index]:
                        z_seg = s_walls[curve_index].z[t_i:t_f]
                        axes.plot(z_seg.real, z_seg.imag, '-')
                        if(plot_data_points == True):
                            axes.plot(z_seg.real, z_seg.imag,
                                        'o', color='b')
        else:
            for s_wall in s_walls:
                if(self.plot_data_points is True):
                    axes.plot(s_wall.z.real, s_wall.z.imag, 'o', color='k')
                if (self.plot_on_cylinder is True):
                    result = numpy.empty(len(s_wall.z), complex)
                    split_at = []
                    result[0] = put_on_cylinder(s_wall.z[0], C)
                    for i, z in enumerate(s_wall.z[1:], start=1):
                        result[i] = put_on_cylinder(z, C)
                        if abs(result[i-1] - result[i]) > pi:
                            split_at.append(i)
                    z_segs = numpy.split(result, split_at)
                    for z_seg in z_segs:
                        axes.plot(z_seg.real, z_seg.imag, '-', color='b',
                                  label=s_wall.label,)
                else:
                    axes.plot(s_wall.z.real, s_wall.z.imag, '-', color='b',
                              label=s_wall.label,)
        axes.set_visible(False)
        self.plots.append(axes)

        return None


    def set_data_cursor(self):
        if self.current_plot_idx is None:
            return None

        # Use a DataCursor to interactively display the label
        # for artists of the current axes.
        self.data_cursor = mpldatacursor.datacursor(
            axes=self.plots[self.current_plot_idx],
            formatter='{label}'.format,
            tolerance=2,
            #hover=True,
            display='single',
        )
    
        return None 


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
            self.data_cursor.hide()

        self.plots[self.current_plot_idx].set_visible(False)
        self.plots[new_plot_idx].set_visible(True)
        # Update the index variable for the currently displayed plot.
        self.current_plot_idx = new_plot_idx
        self.set_data_cursor()
        self.canvas.draw_idle()

        return None


    def show(self):
        plot_idx = 0
        self.current_plot_idx = plot_idx
        self.plots[plot_idx].set_visible(True)
        self.set_data_cursor()

        if(len(self.plots) > 1):
            self.plot_idx_scale = tk.Scale(
                self.toplevel,
                orient=tk.HORIZONTAL,
                to=len(self.plots)-1,
                label='Plot #',
                variable=self.current_plot_idx,
                command=self.scale_action,
            ) 
            self.plot_idx_scale.pack(fill=tk.X)

            self.plot_idx_entry_var.set(plot_idx)
            self.plot_idx_entry = tk.Entry(
                master=self.toplevel,
                textvariable=self.plot_idx_entry_var,
            )
            self.plot_idx_entry.pack()


def plot_segments(
    segments, 
    marked_points=[],
    plot_range=[-5, 5, -5, 5],
    plot_data_points=False
):
    """
    Plot the given segments for debugging.
    """
    x_min, x_max, y_min, y_max = plot_range

    # Plot setting.
    pyplot.xlim(x_min, x_max)
    pyplot.ylim(y_min, y_max)
    pyplot.axes().set_aspect('equal')

    for segment in segments:
        xs, ys = segment
        pyplot.plot(xs, ys, '-')
        if(plot_data_points == True):
            pyplot.plot(xs, ys, 'o', color='b')

    for p in marked_points:
        pyplot.plot(p[0], p[1], 'x', markeredgewidth=2, markersize=8,
                    color='k')

    pyplot.show()
    
