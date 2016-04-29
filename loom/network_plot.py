# import pdb
import mpldatacursor


class NetworkPlotBase(object):
    def __init__(self, matplotlib_figure=None, plot_range=None,):
        self.plots = []
        self.data_cursor = None
        self.current_plot_idx = None
        self.plot_range = plot_range

        self.figure = matplotlib_figure
        if self.figure is not None:
            self.figure.clf()

    def draw(
        self, phase=None, branch_points=None, joints=None, punctures=None,
        irregular_singularities=None, walls=None, wall_segments=None,
        wall_colors=None, labels=None, plot_joints=False,
        plot_data_points=False, branch_cut_rotation=1,
    ):
        """
        branch_points = [[bpx, bpy], ...]
        joints = [[jpx, jpy], ...]
        punctures = [[px, py], ...]
        walls = [s_wall[0].z, s_wall[1].z, ...]
        wall_segments = [[[wall_i_seg_j_start, wall_i_seg_j_stop], ...], ...]
        labels = {'branch_points': [bp1_label, ...],
                  'joints': [jp1_label, ...],
                  'walls': [[wall_i_seg_j_label], ...]}
        """
        rect = [.1, 0.15, .8, .8]

        axes = self.figure.add_axes(
            rect,
            label="Network #{}".format(len(self.plots)),
            aspect='equal',
        )

        if self.plot_range is not None:
            [[x_min, x_max], [y_min, y_max]] = self.plot_range
            axes.set_xlim(x_min, x_max)
            axes.set_ylim(y_min, y_max)
        else:
            axes.autoscale(enable=True, axis='both', tight=None)
            x_min, x_max = axes.get_xlim()
            y_min, y_max = axes.get_ylim()
            self.plot_range = [[x_min, x_max], [y_min, y_max]]

        # Plot wall segments.
        for i, wall in enumerate(wall_segments):
            for j, segment in enumerate(wall):
                start, stop = segment
                zs = walls[i][start:stop]

                if plot_data_points is True:
                    axes.plot(zs.real, zs.imag, 'o', color='k')

                seg_color = wall_colors[i][j]
                axes.plot(zs.real, zs.imag, '-',
                          color=seg_color,
                          label=labels['walls'][i][j],)

        # Plot branch points.
        for i, bp in enumerate(branch_points):
            bpx, bpy = bp
            axes.plot(bpx, bpy, 'x', markeredgewidth=2, markersize=8,
                      color='k',
                      label=labels['branch_points'][i],)

        # Plot irregular singularities
        for i, irr_sing in enumerate(irregular_singularities):
            isx, isy = irr_sing
            axes.plot(isx, isy, 'o', markeredgewidth=2, markersize=8,
                      color='k', markerfacecolor='none',
                      label=labels['irregular_singularities'][i],)

        # Plot branch cuts according to the z-plane rotation.
        y_r = complex((2j * axes.get_ylim()[1]) * branch_cut_rotation)

        for i, bp in enumerate(branch_points):
            bpx, bpy = bp
            axes.plot([bpx, bpx + y_r.real], [bpy, bpy + y_r.imag],
                      ':', color='k',
                      label='Cut of ' + labels['branch_points'][i],)

        for i, irr_sing in enumerate(irregular_singularities):
            isx, isy = irr_sing
            axes.plot([isx, isx + y_r.real], [isy, isy + y_r.imag],
                      ':', color='k',
                      label='Cut of ' + labels['irregular_singularities'][i],)

        # Plot joints.
        if plot_joints is True:
            for i, jp in enumerate(joints):
                jpx, jpy = jp
                axes.plot(jpx, jpy, '+', markeredgewidth=2,
                          markersize=8, color='k', label=labels['joints'][i],)

        # Plot puncturess.
        for i, p in enumerate(punctures):
            px, py = p
            axes.plot(px, py, 'o', markeredgewidth=2, markersize=8,
                      color='k', markerfacecolor='none',
                      label=labels['punctures'][i],)

        axes.set_visible(False)
        self.plots.append(axes)

    def autoscale(self):
        min_x_min = None
        max_x_max = None
        min_y_min = None
        max_y_max = None

        for axes in self.plots:
            x_min, x_max = axes.get_xlim()
            if min_x_min is None or min_x_min > x_min:
                min_x_min = x_min
            if max_x_max is None or max_x_max < x_max:
                max_x_max = x_max

            y_min, y_max = axes.get_ylim()
            if min_y_min is None or min_y_min > y_min:
                min_y_min = y_min
            if max_y_max is None or max_y_max < y_max:
                max_y_max = y_max

        for axes in self.plots:
            axes.set_xlim(min_x_min, max_x_max)
            axes.set_ylim(min_y_min, max_y_max)

    def set_data_cursor(self):
        if self.current_plot_idx is None:
            return None

        # Use a DataCursor to interactively display the label
        # for artists of the current axes.
        self.data_cursor = mpldatacursor.datacursor(
            axes=self.plots[self.current_plot_idx],
            formatter='{label}'.format,
            tolerance=4,
            # hover=True,
            # display='single',
            display='multiple',
            draggable=True,
        )

    def change_current_plot(self, new_plot_idx):
        if new_plot_idx == self.current_plot_idx:
            return None
        elif new_plot_idx < 0:
            new_plot_idx = 0
        elif new_plot_idx > len(self.plots) - 1:
            new_plot_idx = len(self.plots) - 1

        if self.current_plot_idx is not None:
            self.plots[self.current_plot_idx].set_visible(False)
        self.plots[new_plot_idx].set_visible(True)
        # Update the index variable for the currently displayed plot.
        self.current_plot_idx = new_plot_idx

    def show_prev_plot(self, event):
        if self.data_cursor is not None:
            self.data_cursor.hide().disable()
        self.change_current_plot(self.current_plot_idx - 1)
        self.set_data_cursor()
        self.figure.show()

    def show_next_plot(self, event):
        if self.data_cursor is not None:
            self.data_cursor.hide().disable()
        self.change_current_plot(self.current_plot_idx + 1)
        self.set_data_cursor()
        self.figure.show()
