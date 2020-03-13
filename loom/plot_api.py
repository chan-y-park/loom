import numpy
import logging
# import pdb

import matplotlib

from math import pi

from sympy import oo
from matplotlib.backends.backend_agg import FigureCanvas

from misc import put_on_cylinder, get_splits_with_overlap

class NetworkPlot(object):
    def __init__(self, matplotlib_figure=None, plot_range=None,):
        self.plot_range = plot_range
        self.figure = matplotlib_figure
        if self.figure is not None:
            self.figure.clf()

    def draw(
        self,
        axes=None,
        #phase=None,
        branch_points=None,
        joints=None,
        punctures=None,
        irregular_singularities=None,
        walls=None,
        wall_segments=None,
        wall_colors=None,
        labels=None,
        plot_joints=False,
        plot_data_points=False,
        branch_cut_rotation=1,
        markersize=2,
        markeredgewidth=1,
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
            axes.plot(
                bpx, bpy, 'x',
                color='k',
                #markeredgewidth=2, markersize=8,
                label=labels['branch_points'][i],
            )

        # Plot irregular singularities
        for i, irr_sing in enumerate(irregular_singularities):
            isx, isy = irr_sing
            axes.plot(
                isx, isy, 'o',
                #markeredgewidth=2, markersize=8,
                color='k', markerfacecolor='none',
                label=labels['irregular_singularities'][i],
            )

        # Plot branch cuts according to the z-plane rotation.
        y_r = complex((2j * axes.get_ylim()[1]) * branch_cut_rotation)

        for i, bp in enumerate(branch_points):
            bpx, bpy = bp
            axes.plot(
                [bpx, bpx + y_r.real], [bpy, bpy + y_r.imag],
                ':', color='k',
                label='Cut of ' + labels['branch_points'][i],
            )

        for i, irr_sing in enumerate(irregular_singularities):
            isx, isy = irr_sing
            axes.plot([
                isx, isx + y_r.real], [isy, isy + y_r.imag],
                ':', color='k',
                label='Cut of ' + labels['irregular_singularities'][i],
            )

        # Plot joints.
        if plot_joints is True:
            for i, jp in enumerate(joints):
                jpx, jpy = jp
                axes.plot(
                    jpx, jpy, '+',
                    #markeredgewidth=2, markersize=8,
                    color='k', label=labels['joints'][i],
                )

        # Plot puncturess.
        for i, p in enumerate(punctures):
            px, py = p
            axes.plot(
                px, py, 'o',
                #markeredgewidth=2, markersize=8,
                color='k', markerfacecolor='none',
                label=labels['punctures'][i],
            )

        self.figure.tight_layout()


class SpectralNetworkPlot(NetworkPlot):
    def draw(
        self,
        axes=None,
        sw_data=None,
        spectral_network=None,
        soliton_tree=None,
        logger_name='loom',
        plot_joints=False,
        plot_data_points=False,
        plot_on_cylinder=False,
        C=[[1, 0], [0, 1]],
    ):
        logger = logging.getLogger(logger_name)

        labels = {'branch_points': [], 'joints': [], 'punctures': [],
                  'walls': [], 'irregular_singularities': []}

        branch_points = sw_data.branch_points
        punctures = sw_data.regular_punctures + sw_data.irregular_punctures
        irregular_singularities = sw_data.irregular_singularities
        g_data = sw_data.g_data
        branch_cut_rotation = sw_data.branch_cut_rotation

        if self.plot_range is None:
            if plot_on_cylinder is True:
                self.plot_range = [[-pi, pi], [-5, 5]]

        branch_points_z = []
        for i, bp in enumerate(branch_points):
            if plot_on_cylinder is True:
                bp_z = put_on_cylinder(bp.z, C)
            else:
                bp_z = bp.z
            branch_points_z.append([bp_z.real, bp_z.imag])
            labels['branch_points'].append(bp.label)

        joints_z = []
        if plot_joints is True:
            for i, jp in enumerate(spectral_network.joints):
                if plot_on_cylinder is True:
                    jp_z = put_on_cylinder(jp.z, C)
                else:
                    jp_z = jp.z
                joints_z.append([jp_z.real, jp_z.imag])
                labels['joints'].append(jp.label)

        punctures_z = []
        for i, p in enumerate(punctures):
            if p.z == oo:
                continue
            if plot_on_cylinder is True:
                p_z = put_on_cylinder(p.z, C)
            else:
                p_z = p.z
            punctures_z.append([p_z.real, p_z.imag])
            labels['punctures'].append(p.label)

        irregular_singularities_z = []
        for i, irs in enumerate(irregular_singularities):
            if irs.z == oo:
                continue
            if plot_on_cylinder is True:
                irs_z = put_on_cylinder(irs.z, C)
            else:
                irs_z = irs.z
            irregular_singularities_z.append([irs_z.real, irs_z.imag])
            labels['irregular_singularities'].append(irs.label)

        walls = []
        wall_segments = []
        wall_roots = []

        if soliton_tree is not None:
            s_walls = soliton_tree.streets
        else:
            spectral_network.s_walls

        for i, s_wall in enumerate(s_walls):
            seg_labels = []

            if plot_on_cylinder is True:
                # XXX: Incomplete
                # Need to include local root data.
                raise NotImplementedError

                zs_on_cylinder = numpy.fromfunction(
                    lambda i: put_on_cylinder(s_wall.z[i], C),
                    (len(s_wall.z)),
                )
                walls.append(zs_on_cylinder)

                splits = []
                for j, delta_z in enumerate(numpy.diff(zs_on_cylinder)):
                    if abs(delta_z) > pi:
                        splits.append(j)

                wall_segments.append(get_splits_with_overlap(splits))

            else:
                walls.append(s_wall.z)
                wall_segments.append(
                    get_splits_with_overlap(s_wall.get_splits())
                )
                wall_roots.append(s_wall.local_roots)

                seg_labels = [s_wall.label + '\n' + root_str
                              for root_str in map(str, s_wall.local_roots)]
                labels['walls'].append(seg_labels)

        wall_colors = []
        for seg_roots in wall_roots:
            colors = []
            for root in seg_roots:
                color = g_data.get_root_color(root)
                if color is None:
                    color = '#000000'
                    logger.warning(
                        'Unknown root color for {} '
                        'of a spectral network with phase = {}.'
                        .format(s_wall.label, spectral_network.phase)
                    )
                colors.append(color)
            wall_colors.append(colors)

        super(SpectralNetworkPlot, self).draw(
            axes=axes,
            #phase=spectral_network.phase,
            branch_points=branch_points_z,
            joints=joints_z,
            punctures=punctures_z,
            irregular_singularities=irregular_singularities_z,
            walls=walls,
            wall_segments=wall_segments,
            wall_colors=wall_colors,
            labels=labels,
            plot_joints=plot_joints,
            plot_data_points=plot_data_points,
            branch_cut_rotation=branch_cut_rotation,
        )

    def get_legend(
        self,
        sw_data=None,
        spectral_network=None,
    ):
        plot_legend = (
            '\n'
            '------------------------\n'
            'phase : {}\n'.format(spectral_network.phase) +
            '------------------------\n'\
        )
        plot_legend += get_sw_data_legend(sw_data)
        return plot_legend


class SolitonTreePlot(SpectralNetworkPlot):
    def __init__(self, plot_range=None):
        super(SolitonTreePlot, self).__init__(
            matplotlib_figure=matplotlib.figure.Figure(),
            plot_range=plot_range,
        )
        self.canvas = FigureCanvas(
            self.figure
        )

    def show(self):
        self.figure.show()

    def draw(
        self,
        title=None,
        **kwargs
    ):
        #rect = [.1, 0.15, .8, .8]
        #rect = [0, 0, 1, 1]

        #axes = self.figure.add_axes(
        #    rect,
        #    aspect='equal',
        #)
        axes = self.figure.add_subplot(1, 1, 1, aspect='equal')
        axes.set_title(title)

        super(SolitonTreePlot, self).draw(
            axes=axes,
            **kwargs
        )


def get_label(value, dictionary):
    """
    Returns the key of a dictionary entry.
    Values of the dictionary must be numpy arrays.
    """
    is_in_dictionary = False
    for v in dictionary.values():
        if numpy.array_equal(v, value):
            is_in_dictionary = True

    if is_in_dictionary:
        return [
            k for k, v in dictionary.items() if numpy.array_equal(v, value)
        ][0]
    else:
        return 'Not a root: {}'.format(value)


def make_root_dictionary(g_data):
    g_roots = list(g_data.roots)
    root_dictionary = {'alpha_' + str(i): rt for i, rt in enumerate(g_roots)}
    return root_dictionary


def make_weight_dictionary(g_data):
    g_weights = list(g_data.weights)
    weight_dictionary = {'mu_' + str(i): w for i, w in enumerate(g_weights)}
    return weight_dictionary


def get_sw_data_legend(sw_data):
    branch_points = sw_data.branch_points
    regular_punctures = sw_data.regular_punctures
    irregular_singularities = sw_data.irregular_singularities
    g_data = sw_data.g_data

    root_dictionary = make_root_dictionary(g_data)
    weight_dictionary = make_weight_dictionary(g_data)
    root_labels = list(root_dictionary.keys())
    roots = list(root_dictionary.values())
    weight_pairs = [
        [str('(mu_' + str(p[0]) + ', mu_' + str(p[1]) + ')')
         for p in g_data.ordered_weight_pairs(rt)]
        for rt in roots
    ]
    weight_labels = list(weight_dictionary.keys())
    weights = list(weight_dictionary.values())

    legend = ('\t--- Root System ---\n')
    for i in range(len(roots)):
        legend += (
            root_labels[i] + ' : {}\n'.format(list(roots[i])) +
            'ordered weight pairs : {}\n'.format(weight_pairs[i])
        )

    legend += ('\t--- Weight System ---\n')
    for i in range(len(weights)):
        legend += (weight_labels[i] + ' : {}\n'.format(list(weights[i])))

    if regular_punctures is not None:
        legend += ('\t--- Regular punctures ---\n')
        for p in regular_punctures:
            legend += (
                p.label +
                '\tposition : {}\n'.format(p.z)
            )

    if branch_points is not None:
        legend += ('\t--- Branch Points ---\n')
        for bp in branch_points:
            rt_labels = [get_label(rt, root_dictionary)
                         for rt in bp.positive_roots]
            legend += (
                bp.label +
                '\tposition : {}\n'.format(bp.z) +
                '\t\troot type : {}\n'.format(rt_labels) +
                '\t\tmonodromy matrix : \n{}\n'.format(bp.monodromy)
            )

    if irregular_singularities is not None:
        legend += ('\t--- Irregular Singularities ---\n')
        for irs in irregular_singularities:
            legend += (
                irs.label +
                '\tposition : {}\n'.format(irs.z) +
                '\tmonodomry matrix : \n{}\n'.format(irs.monodromy)
            )

    return legend
