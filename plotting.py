import pdb
from matplotlib import pyplot

class SpectralNetworkPlot:
    def __init__(self, 
        spectral_network_data,
        #plot_range=[[-5, 5], [-5, 5]], 
        plot_bins=False, 
        plot_joints=False, 
        plot_data_points=False,
        plot_segments=False,
    ):
        self.plot = pyplot
        self.axis_limits = spectral_network_data['config_data'].z_range_limits

        s_walls = spectral_network_data['s_walls']

        # Give an identifier to the figure we are goint to produce
        label = 'spectral_network_theta_{}'.format(
            spectral_network_data['phase']
        )
        pyplot.figure(label)

        x_min, x_max, y_min, y_max = self.axis_limits

        # Plot setting.
        pyplot.xlim(x_min, x_max)
        pyplot.ylim(y_min, y_max)
        pyplot.axes().set_aspect('equal')

        # Draw a lattice of bins for visualization.
        if(plot_bins is True):
            bin_size = spectral_network_data['hit_table'].get_bin_size()
            xv = x_min
            while xv < x_max:
                xv += bin_size
                pyplot.axvline(x = xv, linewidth=0.5, color='0.75')

            yh = y_min  
            while yh < y_max:
                yh += bin_size
                pyplot.axhline(y = yh, linewidth=0.5, color='0.75')
        # End of drawing the bin lattice.

        # Plot branch points
        for rp in spectral_network_data['sw_curve'].ramification_points:
            bpx = rp.z.real 
            bpy = rp.z.imag 
            pyplot.plot(bpx, bpy, 'x', markeredgewidth=2, markersize=8, 
                        color='k')
        # End of plotting branch points

        # Plot joints
        if(plot_joints == True):
            for jp in spectral_network_data['joints']:
                jpx = jp.z.real 
                jpy = jp.z.imag 
                pyplot.plot(jpx, jpy, '+', markeredgewidth=2, markersize=8, 
                            color='k')
        # End of plotting joints

        # If we have segments of curves, draw them in different colors.
        if(plot_segments == True):
            hit_table = spectral_network_data['hit_table']
            for bin_key in hit_table:
                for curve_index in hit_table[bin_key]:
                    for t_i, t_f in hit_table[bin_key][curve_index]:
                        seg_xcoords = [z.real for z in 
                                       s_walls[curve_index].get_zs(t_i, t_f)] 
                        seg_ycoords = [z.imag for z in 
                                       s_walls[curve_index].get_zs(t_i, t_f)] 
                        pyplot.plot(seg_xcoords, seg_ycoords, '-')
                        if(plot_data_points == True):
                            pyplot.plot(seg_xcoords, seg_ycoords,
                                        'o', color='b')
        else:
            for s_wall in s_walls:
                xcoords = [z.real for z in s_wall.get_zs()] 
                ycoords = [z.imag for z in s_wall.get_zs()] 
                pyplot.plot(xcoords, ycoords, '-', color='b')
                if(plot_data_points == True):
                    pyplot.plot(xcoords, ycoords, 'o', color='b')
        

    def show(self):
        self.plot.show()

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
    
