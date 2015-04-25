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
        self.master = master
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
            resize_callback=self.canvas_resize_callback
        )
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        toolbar = NavigationToolbar(self.canvas, self.toplevel)
        toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    

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

        axes.set_title('phase = ({:.4f})pi'.format(theta/pi))

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
                        axes.plot(z_seg.real, z_seg.imag, '-', 
                                  #color='b',
                                  label=s_wall.label,)
                else:
                    axes.plot(s_wall.z.real, s_wall.z.imag, '-',
                              #color='b',
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


def plot_s_walls(
    s_walls,
    ramification_points=[],
    joints=[],
    plot_range=[-5, 5, -5, 5],
    plot_data_points=False,
    marked_points=[],
    colors=['b', 'g', 'r', 'c', 'm', 'y'], 
):
    x_min, x_max, y_min, y_max = plot_range

    pyplot.figure(1)
    pyplot.title('S-walls')

    # z-plane
    zax = pyplot.subplot(
        121,
        label='z-plane',
        xlim=(x_min, x_max),
        ylim=(y_min, y_max),
        aspect='equal',
    )
    # Plot branch points
    for rp in ramification_points:
        zax.plot(rp.z.real, rp.z.imag, 'x', markeredgewidth=2, markersize=8, 
                 color='k', label=rp.label,)
    for jp in joints:
        zax.plot(jp.z.real, jp.z.imag, '+', markeredgewidth=2, markersize=8, 
                 color='k', label=jp.label,)
    for p in marked_points:
        zax.plot(p[0].real, p[0].imag, 'o', markeredgewidth=2,
                 markersize=4, color='k')

    # x-plane
    xax = pyplot.subplot(
        122,
        label='x-plane',
        xlim=(x_min, x_max),
        ylim=(y_min, y_max),
        aspect='equal',
    )
    for rp in ramification_points:
        xax.plot(rp.x.real, rp.x.imag, 'x', markeredgewidth=2, markersize=8, 
                 color='k', label=rp.label,)
    for jp in joints:
        for j, x_j in enumerate(jp.x):
            xax.plot(
                x_j.real, x_j.imag,
                '+', markeredgewidth=2, markersize=8, color='k', 
                label=(jp.label + ', {}'.format(j)),
            )
    for p in marked_points:
        xax.plot(p[1].real, p[1].imag, 'o', markeredgewidth=2,
                 markersize=4, color='k')

    for i, s_wall in enumerate(s_walls):
        s_wall_color = colors[i % len(colors)]
        # z-plane
        zrs = s_wall.z.real
        zis = s_wall.z.imag
        zax.plot(zrs, zis, '-', color=s_wall_color, 
                 label=s_wall.label)
        if(plot_data_points == True):
            zax.plot(zrs, zis, 'o', color=s_wall_color, label=s_wall.label)

        # x-plane
        xs = s_wall.x.T
        for j, x_j in enumerate(xs):
            xrs = x_j.real
            xis = x_j.imag
            xax.plot(
                xrs, xis, '-', color=s_wall_color, 
                label=(s_wall.label + ',{}'.format(j))
            )
            if(plot_data_points == True):
                xax.plot(
                    xrs, xis, 'o', color=s_wall_color, 
                    label=(s_wall.label + ',{}'.format(j))
                )
            
    mpldatacursor.datacursor(
        formatter='{label}'.format,
        #tolerance=2,
        #hover=True,
        display='multiple',
    )

    pyplot.show()

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
    
