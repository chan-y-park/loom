import pdb
from matplotlib import pyplot
from matplotlib.widgets import Slider


class SpectralNetworkPlot:
    def __init__(self, 
        config,
        plot_bins=False, 
        plot_joints=False,
        plot_data_points=False,
        plot_segments=False
    ):
        # Give an identifier to the figure we are goint to produce
        self.config = config
        self.plot_bins = plot_bins
        self.plot_joints = plot_joints
        self.plot_data_points = plot_data_points
        self.plot_segments = plot_segments

        label = 'spectral_network'
        self.figure = pyplot.figure(label)
        self.plots = []
        self.current_plot = 0

    def set_data(
        self,
        spectral_network_data,
    ):
        theta = spectral_network_data['phase']
        ramification_points = spectral_network_data['ramification_points']
        s_walls = spectral_network_data['s_walls']

        z_range_limits = self.config['z_range_limits']
        if z_range_limits is None:
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
            bpx = rp.z.real 
            bpy = rp.z.imag 
            axes.plot(bpx, bpy, 'x', markeredgewidth=2, markersize=8, 
                      color='k')
        # End of plotting branch points

   
        # Plot joints
        if(self.plot_joints is True):
            for jp in spectral_network_data['joints']:
                jpx = jp.z.real 
                jpy = jp.z.imag 
                axes.plot(jpx, jpy, '+', markeredgewidth=2, markersize=8, 
                            color='k')
        # End of plotting joints

        # If we have segments of curves, draw them in different colors.
        if(self.plot_segments is True):
            hit_table = spectral_network_data['hit_table']
            for bin_key in hit_table:
                for curve_index in hit_table[bin_key]:
                    for t_i, t_f in hit_table[bin_key][curve_index]:
                        seg_xcoords = [z.real for z in 
                                       s_walls[curve_index].get_zs(t_i, t_f)] 
                        seg_ycoords = [z.imag for z in 
                                       s_walls[curve_index].get_zs(t_i, t_f)] 
                        axes.plot(seg_xcoords, seg_ycoords, '-')
                        if(plot_data_points == True):
                            axes.plot(seg_xcoords, seg_ycoords,
                                        'o', color='b')
        else:
            for s_wall in s_walls:
                xcoords = [z.real for z in s_wall.get_zs()] 
                ycoords = [z.imag for z in s_wall.get_zs()] 
                axes.plot(xcoords, ycoords, '-', color='b')
                if(self.plot_data_points is True):
                    axes.plot(xcoords, ycoords, 'o', color='b')
       
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
        xs = []
        ys = []
        for p in segment:
            xs.append(p[0])
            ys.append(p[1])
        pyplot.plot(xs, ys, '-')
        if(plot_data_points == True):
            pyplot.plot(xs, ys, 'o', color='b')

    for p in marked_points:
        pyplot.plot(p[0], p[1], 'x', markeredgewidth=2, markersize=8,
                    color='k')

    pyplot.show()
    
