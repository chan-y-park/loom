import numpy
import pdb
import logging
#
import time
#

from matplotlib import pyplot
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

        label = 'spectral_network'
        self.figure = pyplot.figure(label)
        self.plots = []
        self.current_plot = 0

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
                      color='k')
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
                          color='k')
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
                        axes.plot(z_seg.real, z_seg.imag, '-', color='b')
                else:
                    axes.plot(s_wall.z.real, s_wall.z.imag, '-', color='b')
        axes.set_visible(False)
        self.plots.append(axes)

    def update(self, slider_val):
        new_plot = int(slider_val)
        if(self.current_plot != new_plot):
            self.plots[self.current_plot].set_visible(False)
            self.plots[new_plot].set_visible(True)
            self.current_plot = new_plot
            self.figure.canvas.draw_idle()

    def show(self):
        self.plots[self.current_plot].set_visible(True)
        if(len(self.plots) > 1):
            theta_axes = self.figure.add_axes([.125, .05, .8, .05])
            theta_slider = Slider(theta_axes, 'theta', 0, len(self.plots),
                                  valinit=self.current_plot, valfmt='%d',
                                  closedmax=False)
            theta_slider.on_changed(self.update)
        pyplot.show()

def plot_segments(segments, 
                  marked_points=[],
                  plot_range=[-5, 5, -5, 5],
                  plot_data_points=False):
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
    
