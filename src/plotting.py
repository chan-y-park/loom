import numpy
import pdb
import logging
import Tkinter as tk

import matplotlib
# use() directive must be called before importing matplotlib.pyplot
matplotlib.use('TkAgg')     

import mpldatacursor

from matplotlib import pyplot
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg as FigureCanvas,
    NavigationToolbar2TkAgg as NavigationToolbar,
)
from matplotlib.widgets import Slider
from math import pi

from misc import PSL2C, put_on_cylinder

class SpectralNetworkPlot:
    def __init__(self, 
        config,
        plot_on_cylinder=False,
        plot_bins=False, 
        plot_joints=False,
        plot_data_points=False,
        plot_segments=False,
    ):
        # Give an identifier to the figure we are goint to produce
        self.config = config
        self.plot_on_cylinder = plot_on_cylinder
        self.plot_bins = plot_bins
        self.plot_joints = plot_joints
        self.plot_data_points = plot_data_points
        self.plot_segments = plot_segments
        self.scroll_bar = None
        self.slider_width = None

        self.root = tk.Tk()
        self.root.wm_title('Spectral Network Plot')

        self.figure = matplotlib.figure.Figure()
        self.plots = []
        self.current_plot = 0

        self.canvas = FigureCanvas(self.figure, master=self.root)
        #self.canvas.show()
        self.canvas.get_tk_widget().pack()

        toolbar = NavigationToolbar(self.canvas, self.root)
        toolbar.update()
        self.canvas._tkcanvas.pack()

    def draw(
        self,
        spectral_network,
    ):
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

        axes = self.figure.add_axes(rect,
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

    def scroll(self, *args):
        num_plots = len(self.plots)
        if args[0] == tk.MOVETO:
            f = args[1]
            new_plot = int(f * len(self.plots))
            self.scroll_bar.set(f, f + self.slider_width)
        elif args[0] == tk.SCROLL:
            step = eval(args[1])
            new_plot = self.current_plot + step
            slider_left, slider_right = self.scroll_bar.get()
            self.scroll_bar.set(
                slider_left + step*(self.slider_width),
                slider_right + step*(self.slider_width),
            )
        if(self.current_plot != new_plot):
            self.plots[self.current_plot].set_visible(False)
            self.plots[new_plot].set_visible(True)
            self.current_plot = new_plot

            current_axes = self.plots[self.current_plot]
            # Use a DataCursor to interactively display the label
            # for artists of the current axes.
            mpldatacursor.datacursor(axes=current_axes,
                                     formatter='{label}'.format)
            #self.figure.canvas.draw_idle()
            self.canvas.draw_idle()

        return None

    def show(self):
        current_axes = self.plots[self.current_plot]
        current_axes.set_visible(True)

        # Use a DataCursor to interactively display the label 
        # for artists of the current axes.
        mpldatacursor.datacursor(axes=current_axes, formatter='{label}'.format)

        if(len(self.plots) > 1):
            self.slider_width = 1.0/len(self.plots)
            self.scroll_bar = tk.Scrollbar(
                self.root, 
                command=self.scroll,
                orient=tk.HORIZONTAL,
            )
            self.scroll_bar.pack(fill=tk.X)
            self.scroll_bar.set(0.0, self.slider_width)
        self.canvas.show()
        return None

def plot_segments(segments, 
                  marked_points=[],
                  plot_range=[-5, 5, -5, 5],
                  plot_data_points=False):
    x_min, x_max, y_min, y_max = plot_range

    # Plot setting.
    self.figure.plot.xlim(x_min, x_max)
    self.figure.plot.ylim(y_min, y_max)
    self.figure.plot.axes().set_aspect('equal')

    for segment in segments:
        xs, ys = segment
        self.figure.plot.plot(xs, ys, '-')
        if(plot_data_points == True):
            self.figure.plot.plot(xs, ys, 'o', color='b')

    for p in marked_points:
        self.figure.plot.plot(p[0], p[1], 'x', markeredgewidth=2, markersize=8,
                    color='k')

    self.figure.plot.show()
    
