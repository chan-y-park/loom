from matplotlib import pyplot

prepare_spectral_network_plot(
    spectral_network,
    plot_range=[[-5, 5], [-5, 5]], 
    plot_bins=False, 
    plot_joints=False, 
    plot_data_points=False,
    plot_segments=False,): 

    hit_table = spectral_network.hit_table
    s_walls = spectral_network.s_walls

    # Give an identifier to the figure we are goint to produce
    label = 'spectral_network_theta_{}'.format(spectral_network.theta)
    pyplot.figure(label)
    pyplot.figure(label).clear()

    [[x_min, x_max], [y_min, y_max]] = plot_range

    # Plot setting.
    pyplot.xlim(x_min, x_max)
    pyplot.ylim(y_min, y_max)
    pyplot.axes().set_aspect('equal')

    # Draw a lattice of bins for visualization.
    if(plot_bins is True):
        bin_size = hit_table.get_bin_size()
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
    for rp in spectral_network.sw_curve.ramification_points:
        bpx = rp.z.real 
        bpy = rp.z.imag 
        pyplot.plot(bpx, bpy, 'x', markeredgewidth=2, markersize=8, 
                    color='k')
    # End of plotting branch points
# TODO:
    # Plot intersection points
    if(plot_joints == True):
        for jp in spectral_network.joints:
            ipx = jp.z.real 
            ipy = jp.z.imag 
            pyplot.plot(ipx, ipy, '+', markeredgewidth=2, markersize=8, 
                        color='k')
    # End of plotting intersection points

    # If we have segments of curves, draw them in different colors.
    if(plot_segments == True):
        for bin_key in hit_table:
            for curve_index in hit_table[bin_key]:
                for t_i, t_f in hit_table[bin_key][curve_index]:
                    seg_xcoords, seg_ycoords = \
                        [list(c) for c in \
                            zip(*k_walls[curve_index].coordinates[t_i: t_f+1])
                        ] 
                    pyplot.plot(seg_xcoords, seg_ycoords, '-')
                    if(plot_data_points == True):
                        pyplot.plot(seg_xcoords, seg_ycoords, 'o', color='b')
    else:
        for k_wall in k_walls:
            xcoords, ycoords = [list(c) for c in zip(*k_wall.coordinates)] 
            pyplot.plot(xcoords, ycoords, '-', color='b')
        if(plot_data_points == True):
            pyplot.plot(xcoords, ycoords, 'o', color='b')
    
    return pyplot.figure(label)



