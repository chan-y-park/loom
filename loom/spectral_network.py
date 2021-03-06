import numpy
import ctypes
import logging
import itertools
import json
import cmath

import constants

from cmath import exp
from sympy import oo

from s_wall import (
    SWall, Joint, get_s_wall_seeds, MIN_NUM_OF_DATA_PTS,
)
from s_wall import GrowLibs
from misc import ctor2, r2toc
from misc import nearest_index
from misc import (
    n_nearest_indices, get_turning_points, get_splits_with_overlap,
    get_descendant_roots, sort_roots, get_delta,
)
from intersection import (
    NoIntersection, find_intersection_of_segments,
)
from geometry import BranchPoint
from geometry import get_xs_along_zs
from trivialization import SWDataWithTrivialization
from ctypes_api import libcgal_get_intersections


class Street(SWall):
    def __init__(
        self,
        s_wall=None,
        end_z=None,
        end_t=None,
        tree_label=None,
        label=None,
        logger_name='loom',
    ):
        logger = logging.getLogger(logger_name)

        super(Street, self).__init__(logger_name=logger_name)
        self.label = label

        if s_wall is not None:
            if end_t is None:
                end_t = numpy.argmin(abs(s_wall.z - end_z))

            if end_t == 0:
                logger.warning(
                    '{}, {}: __init__(): end_t == 0.'
                    .format(tree_label, self.label)
                )
                # XXX: To draw this street, include at least two points.
                end_t = 1

            # XXX: Because a joint is not back-inserted into S-walls,
            # neither [:end_t] nor [:end_t+1] is always correct.
            # However, inserting a joint is an expensive operation.
            self.z = s_wall.z[:end_t + 1]
            self.x = s_wall.x[:end_t + 1]
            self.M = s_wall.M[:end_t + 1]
            self.parents = s_wall.parents

            if s_wall.is_trivialized():
                self.trivialize(s_wall)

    def trivialize(self, s_wall):
        if not s_wall.is_trivialized():
            return None

        end_t = len(self.z)
        self.parent_roots = s_wall.parent_roots
        self.cuts_intersections = [
            [br_loc, t, d] for br_loc, t, d in s_wall.cuts_intersections
            if t < end_t
        ]
        n_segs = len(self.cuts_intersections) + 1
        self.local_roots = s_wall.local_roots[:n_segs]
        if s_wall.multiple_local_roots is not None:
            self.multiple_local_roots = (
                s_wall.multiple_local_roots[:n_segs]
            )
        self.local_weight_pairs = s_wall.local_weight_pairs[:n_segs]


class SolitonTree:
    def __init__(
        self,
        root_s_wall=None,
        root_s_wall_end_t=None,
        root_branch_point=None,
        phase=None,
        label=None,
        logger_name='loom',
    ):
        self.phase = phase
        self.logger_name = logger_name
        self.root_branch_point = root_branch_point
        self.label = label

        if root_s_wall is None:
            self.streets = []
        else:
            root_street = Street(
                s_wall=root_s_wall,
                end_t=root_s_wall_end_t,
                tree_label='{} @ {}'.format(self.label, self.phase),
                label=root_s_wall.label,
            )
            self.streets = [root_street]
            self.grow(root_street)

            street_dict = {}
            for street in self.streets:
                street_dict[street.label] = street
            for street in self.streets:
                if len(street.parents) == 1:
                    parent = street.parents[0]
                    if isinstance(parent, BranchPoint):
                        continue
                else:
                    street.parents = [
                        street_dict[parent.label]
                        for parent in street.parents
                    ]
            for i, street in enumerate(self.streets):
                street.label = 'Street #{}'.format(i)
        # stability is the number of streets after an improvement.
        # If it is greater than one, it means an improvement will result
        # in multiple streets; if it is zero, it means the street
        # will disappear after an improvement
        self.stability = None

    def Z(self):
        # Get the value of Z from those of its streets.
        Z = 0
        for street in self.streets:
            Z += street.Z()
        return Z

    def get_json_data(self):
        json_data = {
            'phase': self.phase,
            'root_branch_point': self.root_branch_point.label,
            'streets': [street.get_json_data() for street in self.streets],
            'stability': self.stability
        }

        return json_data

    def set_from_json_data(self, json_data=None, sw_data=None,):
        self.phase = json_data['phase']
        self.root_branch_point = sw_data[json_data['root_branch_point']]

        obj_dict = sw_data.obj_dict.copy()

        self.streets = []
        for street_data in json_data['streets']:
            a_street = Street()
            a_street.set_from_json_data(street_data)
            self.streets.append(a_street)
            obj_dict[a_street.label] = a_street

        for street in self.streets: 
            street.set_refs(obj_dict)

        self.stability = json_data['stability']

    def save(self, file_path):
        with open(file_path, 'wb') as fp:
            json_data = self.get_json_data()
            json.dump(json_data, fp,)

    def load(self, file_path, sw_data):
        with open(file_path, 'r') as fp:
            json_data = json.load(fp)
            self.set_from_json_data(json_data, sw_data)

    def grow(self, street=None):
        # TODO: use misc.build_family_tree()
        for parent in street.parents:
            if isinstance(parent, BranchPoint):
                # Branch points are leaves of this tree.
                continue
            parent_street = Street(
                s_wall=parent,
                end_z=street.z[0],
                tree_label='{} @ {}'.format(self.label, self.phase),
                label=parent.label,
            )
            self.streets.append(parent_street)
            self.grow(parent_street)
                
        return None

    def get_improved_tree(
        self,
        config=None,
        sw_data=None,
        search_radius=None,
        max_n_iters=constants.SOLITON_TREE_MAX_N_ITERS,
        cache_file_path=None,
        logger_name='loom',
        label=None,
        # XXX
        debug=False,
    ):
        logger = logging.getLogger(self.logger_name)
        step_size = config['size_of_small_step']
        n_steps = config['num_of_steps']
        accuracy = config['accuracy']
        if search_radius is None:
            search_radius = config['size_of_bp_neighborhood']

        N, phi_k_n_czes, phi_k_d_czes = sw_data.ffr_curve.get_phi_k_czes()

        tree = self
        root_street = tree.streets[0]
        z_start = root_street.z[-1]
        z_end = self.root_branch_point.z
        min_D_z = abs(z_end - z_start)

        max_gen = root_street.get_generation()
        theta_0 = self.phase
        seed_data = []
        for street in self.streets:
            rp_i = None
            for parent in street.parents:
                if isinstance(parent, BranchPoint):
                    bp = parent
                    rp_i = bp.ffr_ramification_points[0].i
                    for rp in bp.ffr_ramification_points[1:]:
                        if rp_i != rp.i:
                            raise NotImplementedError

            if rp_i is None:
                continue

            z_0 = street.z[0]
            xs_0 = street.x[0]
            seed_data.append([z_0, xs_0, rp_i, bp])

        for nth in range(max_n_iters):
            prev_min_D_z = min_D_z

            n_pts = int(min_D_z / step_size)
            if n_pts < constants.SOLITON_TREE_MIN_N_DZ_STEPS:
                n_pts = constants.SOLITON_TREE_MIN_N_DZ_STEPS

            dz = (z_end - z_start) / (n_pts - 1) 
            zs = numpy.array(
                [z_start + dz * i for i in range(n_pts)],
                dtype=numpy.complex128,
            )
            xs = numpy.empty((n_pts, 2), dtype=numpy.complex128)
            xs[0] = root_street.x[-1]
            get_xs_along_zs(
                N, phi_k_n_czes, phi_k_d_czes, zs, xs, accuracy
            )
#            Zs = numpy.empty(n_pts, dtype=numpy.complex128)
#            Zs[0] = root_street.M[-1] * cmath.exp(tree.phase * 1j)
#            for i in range(len(zs[1:])):
#                Zs[i + 1] = Zs[i] + (xs[i][0] - xs[i][1]) * dz
#
#            root_street.z = numpy.concatenate((root_street.z, zs[1:]))
#            root_street.x = numpy.concatenate((root_street.x, xs[1:]))
#            root_street.M = numpy.concatenate((root_street.M, abs(Zs[1:])))

            leg = SWall()
            leg.z = zs
            leg.x = xs

            Z = tree.Z() + leg.Z()
            theta_n = cmath.phase(Z)
            Delta_theta = theta_n - tree.phase
            logger.debug('Delta_theta = {}'.format(Delta_theta))

            seed_s_walls = []
            for z_0, xs_0, rp_i, bp in seed_data:
                z_n = (
                    exp(rp_i * 1.0j * (theta_n - theta_0) / (rp_i + 1)) * 
                    (z_0 - bp.z)
                ) + bp.z
                zs = numpy.array([z_0, z_n], dtype=numpy.complex128)
                xs_n = numpy.empty((2, 2), dtype=numpy.complex128)
                xs_n[0] = xs_0
                get_xs_along_zs(
                    N, phi_k_n_czes, phi_k_d_czes, zs, xs_n, accuracy
                )
                seed_s_walls.append(
                    SWall(
                        z_0=z_n,
                        x_0=xs_n[1],
                        M_0=0,
                        parents=[bp],
                        parent_roots=[],
                        label='Street #{}'.format(len(seed_s_walls)),
                        n_steps=n_steps,
                        logger_name=logger_name,
                    )
                )

            # Grow a mini spectral network.
            mini_sn = SpectralNetwork(
                phase=theta_n,
                logger_name=self.logger_name,
            )
            mini_sn.grow(
                config=config,
                sw_data=sw_data,
                seed_s_walls=seed_s_walls,
                num_of_iterations=max_gen,
                s_wall_label_prefix='Street',
            )
            trees = mini_sn.find_two_way_streets(
                config=config,
                sw_data=sw_data,
                search_radius=search_radius,
            )
            stability = len(trees)
            if stability == 1:
                new_tree = trees[0]
            elif stability > 1:
                logger.warning(
                    'More than one soliton tree are obtained '
                    'while improving a soliton tree '
                    'with root street = {}, root branch point = {} '
                    'from a spectral network @ phase = {}.'
                    .format(
                        self.streets[0].label,
                        self.root_branch_point.label,
                        self.phase,
                    )
                )
                if debug:
                    return (mini_sn, trees)

                min_Delta_Z = constants.P_INF
                for i, a_tree in enumerate(trees):
                    Delta_Z = abs(a_tree.Z() - Z)
                    if Delta_Z < min_Delta_Z:
                        new_tree = a_tree
                        min_Delta_Z = Delta_Z
                    
            elif stability == 0:
                logger.warning(
                    'Failed at improving a soliton tree '
                    'with root street = {}, root branch point = {} '
                    'from a spectral network @ phase = {}.'
                    .format(
                        self.streets[0].label,
                        self.root_branch_point.label,
                        self.phase,
                    )
                )
                if debug:
                    return (mini_sn, None)

                break

            z_start = new_tree.streets[0].z[-1]
            min_D_z = abs(z_end - z_start)
            logger.debug('prev_min_D_z = {}, min_D_z = {}.'
                         .format(prev_min_D_z, min_D_z))
            #if min_D_z > prev_min_D_z:
            if min_D_z > prev_min_D_z and min_D_z < search_radius:
                break
            else:
                tree = new_tree

            root_street = tree.streets[0]

            if (abs(Delta_theta) < accuracy) and min_D_z < step_size:
                break
#            if min_D_z < step_size:
#                break

        tree.stability = stability
        tree.label = label

        if cache_file_path is not None:
            logger.info('Saving cache data to {}.'.format(cache_file_path))
            tree.save(cache_file_path)

        return tree

    def set_z_rotation(self, z_rotation):
        for street in self.streets:
            street.set_z_rotation(z_rotation)

    def trivialize(self):
        for street in self.streets:
            street.trivialize()

    def is_trivialized(self):
        for street in self.streets:
            if not street.is_trivialized():
                return False
        return True

    def draw_graph(self, file_name='soliton_tree.pdf'):
        # XXX: Move the following import to the beginning of this module
        # if draw() becomes a permanent feature.
        import pygraphviz as pgv
        soliton_tree_graph = pgv.AGraph(directed=True)
        soliton_tree_graph.add_edge(
            self.streets[0].label,
            'root: ' + self.root_branch_point.label,
        )
        for street in self.streets:
            for parent in street.parents:
                soliton_tree_graph.add_edge(
                    street.label,
                    parent.label,
                )
        soliton_tree_graph.layout(args='-Goverlap=false')
        soliton_tree_graph.draw(file_name)

    def get_bps(self):
        my_bps = set([self.root_branch_point])
        for street in self.streets:
            if street.is_primary():
                my_bps.add(street.parents[0])
        return my_bps

    def simeq(self, other, Z_r):
        my_bps = self.get_bps()
        your_bps = other.get_bps()
        
        if my_bps != your_bps:
            return False
        
        #d = self.diff(other)
        delta_Z = abs(self.Z() - other.Z())
        if delta_Z > Z_r:
            return False
        else:
            return True

    def diff(self, other):
        delta_theta = abs((self.phase - other.phase) / cmath.pi )
        delta_Z = abs((self.Z() - other.Z()) / self.Z())
        return (delta_theta, delta_Z)

    def d_l1(self, other):
        d_l1 = 0
        for x in self.diff(other):
            d_l1 += x
        return d_l1
            

class SpectralNetwork:
    def __init__(
        self,
        phase=None,
        logger_name='loom',
    ):
        self.phase = phase
        self.s_walls = []
        # TODO: Currently SpectralNetwork.joints are used
        # only when growing a spectral network.
        # Decide whether to save them as data
        # or discard them after generating a spectral network.
        self.joints = []
        # errors is a list of (error type string, error value tuples).
        self.errors = []
        self.n_finished_s_walls = None
        self.soliton_trees = None
        self.data_attributes = [
            'phase', 's_walls', 'joints', 'errors',
            'n_finished_s_walls',
            'soliton_trees',
        ]

        self.logger_name = logger_name

    def set_z_rotation(self, z_rotation):
        for s_wall in self.s_walls:
            s_wall.set_z_rotation(z_rotation)
        for joint in self.joints:
            joint.set_z_rotation(z_rotation)
        if self.soliton_trees is not None:
            for tree in self.soliton_trees:
                tree.set_z_rotation(z_rotation)

    def save(self, file_path):
        with open(file_path, 'wb') as fp:
            json_data = self.get_json_data()
            json.dump(json_data, fp,)

    def load(self, file_path, sw_data):
        with open(file_path, 'r') as fp:
            json_data = json.load(fp)
            self.set_from_json_data(json_data, sw_data)

    def get_json_data(self):
        """
        Prepare the spectral network data in a JSON-compatible file.
        """
        json_data = {}
        json_data['phase'] = self.phase
        json_data['n_finished_s_walls'] = self.n_finished_s_walls
        json_data['s_walls'] = [s_wall.get_json_data()
                                for s_wall in self.s_walls]
        json_data['joints'] = [joint.get_json_data()
                               for joint in self.joints]
        json_data['errors'] = self.errors
        if self.soliton_trees is not None:
            json_data['soliton_trees'] = [
                tree.get_json_data() for tree in self.soliton_trees
            ]
        return json_data

    def set_from_json_data(self, json_data, sw_data):
        """
        Load the spectral network data from a JSON-compatible file.
        """
        branch_loci = sw_data.branch_points + sw_data.irregular_singularities
        obj_dict = {}
        for p in branch_loci:
            obj_dict[p.label] = p

        self.phase = json_data['phase']
        try:
            self.n_finished_s_walls = json_data['n_finished_s_walls']
        except KeyError:
            pass
        try:
            self.errors = json_data['errors']
        except KeyError:
            pass

        for s_wall_data in json_data['s_walls']:
            an_s_wall = SWall(logger_name=self.logger_name,)
            an_s_wall.set_from_json_data(s_wall_data)
            self.s_walls.append(an_s_wall)
            obj_dict[an_s_wall.label] = an_s_wall

        # Substitute labels with objects
        for s_wall in self.s_walls:
            s_wall.set_refs(obj_dict)

        for joint_data in json_data['joints']:
            a_joint = Joint()
            a_joint.set_from_json_data(joint_data)
            a_joint.parents = [
                obj_dict[parent_label]
                for parent_label in a_joint.parents
            ]
            self.joints.append(a_joint)

        try:
            self.soliton_trees = []
            for tree_data in json_data['soliton_trees']:
                a_tree = SolitonTree()
                a_tree.set_from_json_data(
                    json_data=tree_data, sw_data=sw_data,
                )
                self.soliton_trees.append(a_tree)
        except KeyError:
            self.soliton_trees = None

    def downsample(self, ratio=None):
        for s_wall in self.s_walls:
            s_wall.downsample(ratio=ratio)

    def grow(
        self, config=None, sw_data=None,
        additional_iterations=0,
        additional_n_steps=0,
        new_mass_limit=None,
        cache_file_path=None,
        method=None,
        downsample=False,
        downsample_ratio=None,
        seed_s_walls=None,
        num_of_iterations=None,
        s_wall_label_prefix='S-wall'
    ):
        """
        Grow the spectral network by seeding SWall's
        and then calling SWall.grow() for each S-wall.

        Finding joints and growing S-walls from the joints
        is done for a given iteration, therefore it is possible
        that there is a joint from which an S-wall is not grown
        if the depth of the joint is too deep.
        """
        logger = logging.getLogger(self.logger_name)

        # Determine the intersection-finding algorithm
        # according to the availability of CGAL.
        try:
            get_intersections = libcgal_get_intersections()
            logger.info('Use CGAL to find intersections.')
            use_cgal = True
        except OSError:
            logger.warning('CGAL not available; use interpolation '
                           'to find intersections.')
            get_intersections = find_intersections_of_curves
            use_cgal = False

        s_wall_grow_libs = GrowLibs(
            config=config,
            sw_data=sw_data,
            phase=self.phase,
            logger_name=self.logger_name,
        )

        use_default_msg = None
        if (
            method == constants.LIB_C and
            not s_wall_grow_libs.ctypes_s_wall.is_available()
        ):
            use_default_msg = 'C libraries'
        elif(
            method == constants.LIB_NUMBA and
            s_wall_grow_libs.numba_grow is None
        ):
            use_default_msg = 'Numba'
        if use_default_msg is not None:
            logger.warning(
                '{} not available, use a default library instead.'
                .format(use_default_msg)
            )
            method = s_wall_grow_libs.default_lib

        # Gather z-coordinates of punctures and branch points.
        ppzs = [
            p.z for p in
            sw_data.regular_punctures + sw_data.irregular_punctures
            if p.z != oo
        ]
        bpzs = [
            p.z for p in
            #sw_data.ffr_ramification_points + sw_data.irregular_punctures
            sw_data.ffr_ramification_points
            if p.z != oo
        ]

        if num_of_iterations is None:
            num_of_iterations = config['num_of_iterations']
        n_steps = config['num_of_steps']
        if(
            additional_n_steps == 0 and
            additional_iterations == 0 and
            new_mass_limit is None
        ):
            self.n_finished_s_walls = 0
        else:
            # self.n_finished_s_walls is loaded from data
            # and is not None.
            num_of_iterations = 1 + additional_iterations

        iteration = 1
        # Start iterations.
        while(iteration <= num_of_iterations):
            logger.info('Start iteration #{}...'.format(iteration))
            """
            Iterate until there is no new joint
            or for a specified number of iterations.

            Each iteration starts with seeding S-walls by either
            a) seeding them around each branch point,
               if this spectral network is a brand new one, or
            b) using the endpoints of S-walls as seeds
               if this is extending an existing S-walls.
            """

            # Seed S-walls.
            if (seed_s_walls is not None and iteration == 1):
                # Use given seed S-walls.
                new_s_walls = seed_s_walls
            elif (
                additional_n_steps == 0 and additional_iterations == 0 and
                new_mass_limit is None and iteration == 1
            ):
                # S-walls that will be grown at this iteration.
                new_s_walls = []
                # Seed S-walls around each branch point.
                logger.info('Seed S-walls at branch points...')
                for bp in sw_data.branch_points:
                    try:
                        s_wall_seeds = get_s_wall_seeds(
                            sw_data, self.phase, bp, config, self.logger_name,
                        )
                        for z_0, x_0, M_0 in s_wall_seeds:
                            label = '{} #{}'.format(
                                s_wall_label_prefix, len(new_s_walls)
                            )
                            new_s_walls.append(
                                SWall(
                                    z_0=z_0,
                                    x_0=x_0,
                                    M_0=M_0,
                                    parents=[bp],
                                    parent_roots=[root for root
                                                  in bp.positive_roots],
                                    label=label,
                                    n_steps=n_steps,
                                    logger_name=self.logger_name,
                                )
                            )
                    except RuntimeError as e:
                        error_msg = (
                            'Error while seeding S-walls at {}: {}\n'
                            'Skip the seeding.'
                            .format(bp.label, e)
                        )
                        logger.error(error_msg)
                        self.errors.append(('RuntimeError', error_msg))
                        continue

            elif (
                (additional_n_steps > 0 or new_mass_limit is not None) and
                iteration == 1
            ):
                # Use the endpoint of each S-wall as a seed.
                msg = 'Extending S-walls'
                if additional_n_steps > 0:
                    msg += (
                        ' by {} steps'
                        .format(additional_n_steps)
                    )
                if new_mass_limit is not None:
                    msg += (
                        ' with new mass limit {}'
                        .format(new_mass_limit)
                    )
                logger.info(msg + '.')

                new_s_walls = []
                for s_wall in self.s_walls:
                    prev_array_size = len(s_wall.z)
                    if prev_array_size == 0:
                        logger.warning(
                            '{} has zero data.'.format(s_wall)
                        )
                        continue
                    new_n_steps = (
                        n_steps - prev_array_size + additional_n_steps
                    )
                    if new_n_steps > MIN_NUM_OF_DATA_PTS:
                        new_s_walls.append(
                            SWall(
                                z_0=s_wall.z[-1],
                                x_0=s_wall.x[-1],
                                M_0=s_wall.M[-1],
                                parents=s_wall.parents,
                                parent_roots=s_wall.parent_roots,
                                label=s_wall.label,
                                n_steps=new_n_steps,
                                logger_name=self.logger_name,
                            )
                        )
            elif additional_iterations > 0 and iteration == 1:
                logger.info(
                    'Do {} additional iteration(s) to find joints '
                    'and grow S-walls from them.'
                    .format(additional_iterations)
                )
                new_s_walls = []

            # Grow each newly-seeded S-wall.
            i = 0
            while (i < len(new_s_walls)):
                s_i = new_s_walls[i]
                try:
                    s_i.grow(
                        branch_point_zs=bpzs,
                        puncture_point_zs=ppzs,
                        config=config,
                        libs=s_wall_grow_libs,
                        use_scipy_ode=config['use_scipy_ode'],
                        twist_lines=sw_data.twist_lines,
                        method=method,
                    )

                    if len(s_i.z) < MIN_NUM_OF_DATA_PTS:
                        logger.warning(
                            '{} has only {} data point(s); remove this S-wall.'
                            .format(s_i, len(s_i).z)
                        )
                        new_s_walls.pop(i)
                        continue

                except RuntimeError as e:
                    error_msg = (
                        'Error while growing {}: {}\n'
                        'Stop growing this S-wall.'
                        .format(s_i.label, e)
                    )
                    logger.error(error_msg)
                    self.errors.append(('RuntimeError', error_msg))

                if sw_data.is_trivialized():
                    # Cut the grown S-walls
                    # at the intersetions with branch cuts
                    # and decorate each segment with its root data.
                    try:
                        root_types = s_i.determine_root_types(
                            sw_data,
                            cutoff_radius=config['size_of_small_step'],
                        )
                        if (
                            root_types == 'Rebuild S-wall' and
                            config['use_scipy_ode'] is True
                        ):
                            logger.info(
                                'Grow this S-wall again, '
                                'using manual integration.'
                            )
                            s_i.grow(
                                branch_point_zs=bpzs,
                                puncture_point_zs=ppzs,
                                config=config,
                                libs=s_wall_grow_libs,
                                use_scipy_ode=False,
                                method=method,
                            )
                            root_types = s_i.determine_root_types(
                                sw_data,
                                cutoff_radius=config['size_of_small_step'],
                            )
                            if root_types == 'Rebuild S-wall':
                                logger.warning(
                                    'Could not determine the root '
                                    'types of this wall even manually. Likely '
                                    'numerical failure.'
                                )
                            raise RuntimeError
                    except RuntimeError as e:
                        error_msg = (
                            'Error while determining root types of {}: {}\n'
                            'Remove this S-wall.'
                            .format(s_i.label, e)
                        )
                        logger.error(error_msg)
                        self.errors.append(
                            ('RuntimeError', error_msg)
                        )
                        # Remove the S-wall.
                        new_s_walls.pop(i)
                        continue

                # End of growing the i-th new S-wall.
                i += 1

            logger.info(
                'Growing S-walls in iteration #{} finished.'
                .format(iteration)
            )

            # Find joints between the new S-wall and the previous S-walls,
            # and among the new S-walls.
            new_joints = []     # New joints found in each iteration.
            if iteration < num_of_iterations:
                # S-walls that are already searched for joints.
                finished_s_walls = self.s_walls[:self.n_finished_s_walls]
                # S-walls that are not searched for joints.
                unfinished_s_walls = (
                    self.s_walls[self.n_finished_s_walls:] + new_s_walls
                )

                all_s_walls = unfinished_s_walls + finished_s_walls
                # Now for each of the new (unfinished) walls, we check its
                # joints with other unfinished S-walls that come after it,
                # as well as with all the old (finished) walls.
                # This corresponds to the list slicing
                # all_s_walls[m + 1:]
                for m, unfinished_s_wall in enumerate(unfinished_s_walls):
                    try:
                        new_joints += self.get_new_joints(
                            unfinished_s_wall, all_s_walls[m + 1:],
                            config, sw_data, get_intersections, use_cgal,
                        )
                    except RuntimeError as e:
                        error_msg = (
                            'Error while finding joints from {}: {}\n'
                            'Stop finding joints from this S-wall.'
                            .format(unfinished_s_wall, e)
                        )
                        logger.error(error_msg)
                        self.errors.append(
                            ('RuntimeError', error_msg)
                        )
                    finished_s_walls.append(unfinished_s_wall)

                if (
                    (additional_n_steps > 0 or new_mass_limit is not None) and
                    iteration == 1
                ):
                    self.n_finished_s_walls = len(self.s_walls)
                else:
                    self.n_finished_s_walls += len(unfinished_s_walls)

            # Add the new S-walls to the spectral network
            if (
                (additional_n_steps > 0 or new_mass_limit is not None) and
                iteration == 1
            ):
                prev_s_walls = {}
                for s_wall in self.s_walls:
                    prev_s_walls[s_wall.label] = s_wall

                for nsw in new_s_walls:
                    # Attach new S-walls to existing S-walls
                    try:
                        psw = prev_s_walls[nsw.label]
                    except KeyError:
                        raise RuntimeError(
                            'S-wall extension mismatch: '
                            'cannot attach a new {}.'
                            .format(nsw.label)
                        )

                    psw_n_t = len(psw.z)
                    psw.z = numpy.concatenate((psw.z, nsw.z[1:]))
                    psw.x = numpy.concatenate((psw.x, nsw.x[1:]))
                    psw.M = numpy.concatenate((psw.M, nsw.M[1:]))

                    if sw_data.is_trivialized():
                        psw.local_roots += nsw.local_roots[1:]
                        psw.multiple_local_roots += (
                            nsw.multiple_local_roots[1:]
                        )
                        psw.local_weight_pairs += nsw.local_weight_pairs[1:]

                        for ci in nsw.cuts_intersections:
                            bp, t, d = ci
                            psw.cuts_intersections.append(
                                [bp, psw_n_t + t, d]
                            )

                    logger.info(
                        'Extended {} by {} steps.'
                        .format(nsw.label, len(nsw.z))
                    )
            else:
                self.s_walls += new_s_walls

            if(len(new_joints) == 0):
                if iteration < num_of_iterations:
                    logger.info(
                        'No additional joint found: '
                        'Stop growing this spectral network '
                        'at iteration #{}.'.format(iteration)
                    )
                    break
            else:
                logger.info(
                    'Found {} new joints.'
                    .format(len(new_joints))
                )

            new_s_walls = []
            # Seed an S-wall for each new joint.
            for joint in new_joints:
                joint.label = 'joint #{}'.format(len(self.joints))
                self.joints.append(joint)
                label = '{} #{}'.format(
                    s_wall_label_prefix,
                    len(self.s_walls) + len(new_s_walls),
                )
                if (
                    config['mass_limit'] is None or
                    joint.M < config['mass_limit']
                ):
                    new_s_walls.append(
                        SWall(
                            z_0=joint.z,
                            # The numerics of S-walls involves sheets
                            # from the first fundamental cover.
                            x_0=joint.ode_xs,
                            M_0=joint.M,
                            parents=joint.parents,
                            # XXX: parent_roots != (roots of parents)
                            parent_roots=joint.roots,
                            label=label,
                            n_steps=n_steps,
                            logger_name=self.logger_name,
                        )
                    )

            logger.info('Iteration #{} finished.'.format(iteration))
            iteration += 1

        logger.info('Finished growing a spectral network at phase = {}'
                    .format(self.phase))

        if downsample:
            self.downsample(ratio=downsample_ratio)

        if cache_file_path is not None:
            logger.info('Saving cache data to {}.'.format(cache_file_path))
            self.save(cache_file_path)

    def get_new_joints(
        self, new_s_wall, prev_s_walls, config, sw_data,
        get_intersections, use_cgal,
    ):
        """
        Find intersections between S-walls using
        either CGAL 2d curve intersection or scipy interpolation
        according to the availability, then form joints
        from the intersection points.
        """

        logger = logging.getLogger(self.logger_name)
        accuracy = config['accuracy']

        if (config['root_system'] in ['A1', ]):
            logger.info(
                'There is no joint for the given root system {}.'
                .format(config['root_system'])
            )
            return []

        new_joints = []
        for prev_s_wall in prev_s_walls:
            # First check if the two S-walls are compatible
            # for forming a joint.

            # 1. Check if the new S-wall is a descendant
            # of an existing S-wall.
            if prev_s_wall in new_s_wall.parents:
                continue

            # 2. Split the two S-walls into segments
            # according to the trivialization, then
            # check the compatibility of a pair
            # of segments.
            n_z_splits = new_s_wall.get_splits(endpoints=True)
            num_n_z_segs = len(n_z_splits) - 1
            p_z_splits = prev_s_wall.get_splits(endpoints=True)
            num_p_z_segs = len(p_z_splits) - 1

            for n_z_seg_i, p_z_seg_i in itertools.product(
                range(num_n_z_segs), range(num_p_z_segs)
            ):
                n_z_i = n_z_splits[n_z_seg_i]
                n_z_f = n_z_splits[n_z_seg_i + 1]
                p_z_i = p_z_splits[p_z_seg_i]
                p_z_f = p_z_splits[p_z_seg_i + 1]

                if sw_data.is_trivialized():
                    descendant_roots = get_descendant_roots(
                        (prev_s_wall.multiple_local_roots[p_z_seg_i] +
                         new_s_wall.multiple_local_roots[n_z_seg_i]),
                        sw_data.g_data,
                    )
                    if len(descendant_roots) == 0:
                        # The two segments are not compatible for
                        # forming a joint.
                        continue

                # Find an intersection of the two segments on the z-plane.
                if use_cgal is True:
                    buffer_size = 10
                    intersection_search_finished = False
                    # logger.debug(
                    #     'Finding intersections between {} [{}:{}] '
                    #     'and {} [{}:{}].'
                    #     .format(
                    #         new_s_wall.label, n_z_i, n_z_f + 1,
                    #         prev_s_wall.label, p_z_i, p_z_f + 1,
                    #     )
                    # )

                    while not intersection_search_finished:
                        intersections = numpy.empty((buffer_size, 2),
                                                    dtype=numpy.float64)
                        num_of_intersections = get_intersections(
                            new_s_wall.z[n_z_i:n_z_f + 1],
                            ctypes.c_long(n_z_f + 1 - n_z_i),
                            prev_s_wall.z[p_z_i:p_z_f + 1],
                            ctypes.c_long(p_z_f + 1 - p_z_i),
                            intersections, buffer_size
                        )
                        if num_of_intersections == 0:
                            intersections = []
                            intersection_search_finished = True
                        elif num_of_intersections > buffer_size:
                            logger.info('Number of intersections larger than '
                                        'the buffer size; increase its size '
                                        'to {} and find intersections again.'
                                        .format(num_of_intersections))
                            buffer_size = num_of_intersections
                        else:
                            intersections.resize((num_of_intersections, 2))
                            intersection_search_finished = True

                else:
                    intersections = get_intersections(
                        new_s_wall.z[n_z_i:n_z_f + 1],
                        prev_s_wall.z[p_z_i:p_z_f + 1],
                        accuracy,
                    )

                for ip_x, ip_y in intersections:
                    ip_z = ip_x + 1j * ip_y

                    # Discard apparent intersections of sibling S-walls
                    # that emanate from the same branch point, if they occur
                    # at the beginning of the S-walls
                    # Also discard any intersectins that occur near the
                    # beginning of an S-wall
                    if (
                        # XXX: Do we need to check if the parent
                        # is a branch point?
                        prev_s_wall.parents == new_s_wall.parents and
                        abs(ip_z - prev_s_wall.z[0]) < accuracy and
                        abs(ip_z - new_s_wall.z[0]) < accuracy
                    ):
                        continue
                    elif (
                        nearest_index(prev_s_wall.z, ip_z) == 0 or
                        nearest_index(new_s_wall.z, ip_z) == 0
                    ):
                        continue
                    else:
                        # t_n: index of new_s_wall.z nearest to ip_z
                        t_n = get_nearest_point_index(
                            new_s_wall.z, ip_z, sw_data.branch_points,
                            accuracy,
                        )

                        # t_p: index of z_seg_p nearest to ip_z
                        t_p = get_nearest_point_index(
                            prev_s_wall.z, ip_z, sw_data.branch_points,
                            accuracy,
                        )

                    # TODO: need to put the joint into the parent
                    # S-walls?

                    # logger.debug('Intersection at z = {}'.format(ip_z))

                    dxs = (get_delta(new_s_wall.x, t_n).tolist() +
                           get_delta(prev_s_wall.x, t_p).tolist())
                    if sw_data.is_trivialized():
                        # TODO: check if the following descendant-roots
                        # finding is necessary, note that we calculate
                        # descendant roots above.
                        descendant_roots = get_descendant_roots(
                            (prev_s_wall.get_roots_at_t(t_p) +
                             new_s_wall.get_roots_at_t(t_n)),
                            sw_data.g_data,
                        )
                        joint_data_groups = get_joint_data_groups_by_roots(
                            descendant_roots, ip_z, sw_data, accuracy,
                        )
                    else:
                        joint_data_groups = get_joint_data_groups_from_xs(
                            new_s_wall.x[t_n],
                            prev_s_wall.x[t_p],
                            sw_data.g_data.type,
                            accuracy,
                            dxs,
                        )

                    joint_M = prev_s_wall.M[t_p] + new_s_wall.M[t_n]
                    joint_parents = [prev_s_wall, new_s_wall]

                    for roots, ode_xs in joint_data_groups:
                        if len(roots) > 0:
                            joint_roots = sort_roots(roots, sw_data.g_data)
                        else:
                            joint_roots = []

                        a_joint = Joint(
                            z=ip_z,
                            M=joint_M,
                            ode_xs=ode_xs,
                            parents=joint_parents,
                            roots=joint_roots,
                        )
                        joint_is_new = True
                        for old_joint in self.joints:
                            if a_joint.is_equal_to(
                                old_joint,
                                accuracy,
                                max([abs(dx) for dx in dxs]),
                            ):
                                joint_is_new = False
                                break
                        if joint_is_new:
                            new_joints.append(a_joint)

        return new_joints

    def trivialize(
        self, config, sw_data,
        cache_file_path=None,
        z_plane_rotation=None,
        two_way_streets_only=False,
    ):
        logger = logging.getLogger(self.logger_name)
        accuracy = config['accuracy']

        # Check the trivialized sw data.
        if isinstance(sw_data, SWDataWithTrivialization) is False:
            raise RuntimeError('Need a trivalized Seiberg-Witten data.')

        # Rotate the network properly.
        if z_plane_rotation is not None:
            self.set_z_rotation(1 / z_plane_rotation)

        ppzs = [
            p.z for p in
            sw_data.regular_punctures + sw_data.irregular_punctures
            if p.z != oo
        ]
        bpzs = [
            p.z for p in
            #sw_data.ffr_ramification_points + sw_data.irregular_punctures
            sw_data.ffr_ramification_points
            if p.z != oo
        ]

        # Set the parent roots of S-walls and
        # determine the root type of each S-wall.
        s_wall_grow_libs = GrowLibs(
            config=config,
            sw_data=sw_data,
            phase=self.phase,
            logger_name=self.logger_name,
        )

        if two_way_streets_only:
            if self.soliton_trees is None:
                raise RuntimeError

            s_wall_labels = []
            for tree in self.soliton_trees:
                for street in tree.streets:
                    if street.label not in s_wall_labels:
                        s_wall_labels.append(street.label)
            # NOTE: In principle each S-wall should be trivialized
            # after their parents are trivialized. This is implicitly done
            # in the following.
            s_walls = [
                s_wall for s_wall in self.s_walls
                if s_wall.label in s_wall_labels
            ]
        else:
            s_walls = self.s_walls

        for s_i in s_walls:
            if s_i.is_trivialized():
                continue

            if len(s_i.parents) == 1:
                parent = s_i.parents[0]
                if isinstance(parent, BranchPoint):
                    s_i.parent_roots = [root for root in parent.positive_roots]
                else:
                    raise RuntimeError(
                        '{} has an invalid parent: {}.'
                        .format(s_i.label, parent.label)
                    )
            elif len(s_i.parents) == 2:
                parent_roots = []
                ip_z = s_i.z[0]
                for parent in s_i.parents:
                    if not(isinstance(parent, SWall)):
                        raise RuntimeError(
                            '{} has an invalid parent: {}.'
                            .format(s_i.label, parent.label)
                        )
                    if not parent.is_trivialized():
                        continue
                    t = get_nearest_point_index(
                        parent.z, ip_z, sw_data.branch_points, accuracy,
                    )
                    parent_roots += parent.get_roots_at_t(t)
                if len(parent_roots) == 0:
                    continue
                descendant_roots = get_descendant_roots(
                    parent_roots, sw_data.g_data,
                )
                if len(descendant_roots) != 1:
                    raise NotImplementedError
                else:
                    # NOTE: parent_roots != (roots of parents)
                    s_i.parent_roots = descendant_roots
            else:
                raise RuntimeError(
                    '{} has invalid parents: {}.'
                    .format(
                        s_i.label,
                        [parent.label for parent in s_i.parents]
                    )
                )

            try:
                root_types = s_i.determine_root_types(
                    sw_data,
                    cutoff_radius=config['size_of_small_step'],
                )
                if (
                    root_types == 'Rebuild S-wall' and
                    config['use_scipy_ode'] is True
                ):
                    logger.info(
                        'Grow this S-wall again, '
                        'using manual integration.'
                    )
                    # TODO: Need to reset numpy arrays of s_i.
                    s_i.grow(
                        branch_point_zs=bpzs,
                        puncture_point_zs=ppzs,
                        config=config,
                        libs=s_wall_grow_libs,
                        use_scipy_ode=False,
                    )
                    root_types = s_i.determine_root_types(
                        sw_data,
                        cutoff_radius=config['size_of_small_step'],
                    )
                    if root_types == 'Rebuild S-wall':
                        logger.warning(
                            'Could not determine the root '
                            'types of this wall even manually. Likely '
                            'numerical failure.'
                        )
                    raise RuntimeError
            except RuntimeError as e:
                error_msg = (
                    'Error while determining root types of {}: {}'
                    .format(s_i.label, e)
                )
                logger.error(error_msg)
                self.errors.append(
                    ('RuntimeError', error_msg)
                )
                continue

        # XXX: is this necessary?
        # determine the root type of each joint.

        if two_way_streets_only:
            for tree in self.soliton_trees:
                tree.trivialize()

        if cache_file_path is not None:
            logger.info('Saving cache data to {}.'.format(cache_file_path))
            self.save(cache_file_path)

    def find_two_way_streets(
        self, config=None, sw_data=None,
        search_radius=None,
        cache_file_path=None,
    ):
        """
        Find soliton trees that make up two-way streets
        of the spectral network, set the streets attribute,
        and return the soliton trees.
        """
        logger = logging.getLogger(self.logger_name)

        if search_radius is None:
            search_radius = config['size_of_bp_neighborhood']

        soliton_trees = []
        # Search for the root of a soliton tree.
        for s_wall in self.s_walls:
            for bp in sw_data.branch_points:
                tree = None

                if bp not in s_wall.parents:
                    min_t = numpy.argmin(abs(s_wall.z - bp.z))
                else:
                    # Skip if the S-wall is a child of the branch point.
                    # unless it forms a loop and comes back
                    # to the branch point.
                    tps = get_turning_points(s_wall.z)
                    if len(tps) >= 3:
                        min_t = (
                            (numpy.argmin(abs(s_wall.z[tps[2]:] - bp.z)) +
                                tps[2])
                        )
                    else:
                        continue

                z_t = s_wall.z[min_t]

                if (abs(z_t - bp.z) > search_radius):
                    continue

                if s_wall.is_trivialized():
                    bp_roots = (
                        bp.positive_roots.tolist() +
                        (-bp.positive_roots).tolist()
                    )
                    s_wall_roots = s_wall.get_roots_at_t(min_t)
                    if(len(s_wall_roots) > 1):
                        raise NotImplementedError
                    s_wall_root = s_wall_roots[0].tolist()
                    # XXX: If a branch cut goes through the z[min_t]
                    # then the following criterion may not work.
                    # Consider checking ODE x's?
                    if s_wall_root not in bp_roots:
                        continue
                else:
                    found = False
                    x1_t, x2_t = s_wall.x[min_t]
                    curve_xs = sw_data.ffr_curve.get_xs(z_t)

                    sheet_1 = abs(curve_xs - x1_t).argmin()
                    sheet_2 = abs(curve_xs - x2_t).argmin()
                    if sheet_1 == sheet_2:
                        logger.warning(
                            'Failed in finding sheets of {} '
                            'at {} near {} at phase ={}:\n'
                            '\txs[t = {}] = {}\n'
                            '\tcurve_xs = {}.'
                            .format(
                                s_wall.label, z_t, bp.label, self.phase,
                                min_t, (x1_t, x2_t), curve_xs,
                            )
                        )
                    for rp in bp.ffr_ramification_points:
                        rp_type = rp.ramification_type
                        if not(
                            rp_type == 'type_I'
                            # or rp_type == 'type_II'
                            # or rp_type == 'type_III'
                            # or rp_type == 'type_AD'
                        ):
                            raise NotImplementedError
                        rp_x = rp.x
                        sheets = n_nearest_indices(
                            sw_data.ffr_curve.get_xs(z_t), rp_x, 2,
                        )
                        if (sheet_1 in sheets and sheet_2 in sheets):
                            if found:
                                raise RuntimeError(
                                    'S-wall matched with more than '
                                    'one ramification point.'
                                )
                            else:
                                found = True
                    if not found:
                        continue

                # Found a root of this S-wall tree.
                try:
                    tree = SolitonTree(
                        root_s_wall=s_wall,
                        root_s_wall_end_t=min_t,
                        root_branch_point=bp,
                        phase=self.phase,
                        label='Soliton Tree #{}'.format(len(soliton_trees)),
                    )
                except RuntimeError as e:
                    logger.warning(str(e))
                    logger.warning(
                        'Finding two-way streets '
                        'for a spectral network at phase = {} '
                        'using {} and {} failed.'
                        .format(self.phase, s_wall.label, bp.label)
                    )
                    continue

                if tree is not None:
                    soliton_trees.append(tree)

        self.soliton_trees = soliton_trees

        if cache_file_path is not None:
            logger.info('Saving cache data to {}.'.format(cache_file_path))
            self.save(cache_file_path)

        return soliton_trees


# End of class SpectralNetwork


def get_nearest_point_index(s_wall_z, p_z, branch_points, accuracy,
                            logger_name='loom',):
    """
    Get the index of the point on the S-wall that is nearest to
    the given point on the z-plane.

    When the point found is within the accuracy limit from a branch cut,
    look fot the next nearest point and return its index.
    """
    logger = logging.getLogger(logger_name)

    t_0 = nearest_index(s_wall_z, p_z)

    t = t_0

    t_max = len(s_wall_z) - 1

    t_before = t_0

    t_after = t_0

    for bp in branch_points:
        if bp.z.real - p_z.real == 0.0 and bp.z.imag <= p_z.imag:
            raise RuntimeError(
                'get_nearest_point_index(): '
                'Intersection point lies exactly above {}. '
                'Try perturbing the phase.'
                .format(bp.label)
            )

        elif (
            abs(s_wall_z[t].real - bp.z.real) < accuracy and
            s_wall_z[t].imag - bp.z.imag > 0
        ):
            logger.info(
                'The nearest point is too close to a branch cut, '
                'find the next nearest point.'
            )

            # Check the points before & after the given point
            t_m = t_p = t
            while t_m > 0 or t_p < t_max:
                t_m -= 1
                if t_m >= 0:
                    z_m = s_wall_z[t_m]
                    if (
                        (p_z.real - bp.z.real) * (z_m.real - bp.z.real) > 0 and
                        abs(s_wall_z[t_m].real - bp.z.real) > accuracy and
                        t_m <= t_before
                    ):
                        t_before = t_m

                        break

                t_p += 1
                if t_p <= t_max:
                    z_p = s_wall_z[t_p]
                    if (
                        (p_z.real - bp.z.real) * (z_p.real - bp.z.real) > 0 and
                        abs(s_wall_z[t_p].real - bp.z.real) > accuracy and
                        t_p >= t_after
                    ):
                        t_after = t_p
                        break

            if t_m == 0 and t_p == t_max:
                logger.warning(
                    'Unable to find the next nearest point '
                    'that is on the same side from the branch cut '
                    'as the reference point.'
                )
                break

    if t_before == 0 and t_after == t_max:
        raise RuntimeError(
            'get_nearest_point_index(): '
            'Could not find a suitable nearest point, '
            'due to the presence of some branch cut.'
        )

    elif t_before >= 0:
        t = t_before
    elif t_after <= t_max:
        t = t_max

    # Final check.
    for bp in branch_points:
        if abs(s_wall_z[t].real - bp.z.real) < accuracy:
            logger.warning('Unable to find the nearest point on the S-wall '
                           'that is outside the accuracy limit '
                           'from branch cuts of {}.'.format(bp.label))
    return t


def find_intersections_of_curves(a_zs, b_zs, accuracy):

    a_tps = get_turning_points(a_zs)
    a_z_segs = get_splits_with_overlap(a_tps)

    b_tps = get_turning_points(b_zs)
    b_z_segs = get_splits_with_overlap(b_tps)

    intersections = []

    for a_start, a_stop in a_z_segs:
        a_seg = a_zs[a_start:a_stop]
        for b_start, b_stop in b_z_segs:
            b_seg = a_zs[b_start:b_stop]
            # Find an intersection on the z-plane.
            try:
                ip_x, ip_y = find_intersection_of_segments(
                    (a_seg.real, a_seg.imag),
                    (b_seg.real, b_seg.imag),
                    accuracy,
                )
                intersections.append((ip_x, ip_y))

            except NoIntersection:
                pass

    return intersections


def get_joint_data_groups_from_xs(
    new_s_wall_xs, prev_s_wall_xs, g_type, x_accuracy, dxs=None,
):
    ode_xs = None

    # dxs = (delta(n_x1), delta(n_x2), delta(p_x1), delta(p_x2)),
    # which shows the change of xi's near the intersection along the s-walls.
    if dxs is None:
        dnx1, dnx2, dpx1, dpx2 = (x_accuracy,) * 4
    else:
        dnx1, dnx2, dpx1, dpx2 = map(abs, dxs)

    nx1, nx2 = new_s_wall_xs
    px1, px2 = prev_s_wall_xs

    if g_type == 'A':
        if(
            (abs(nx1 - px2) < max(dnx1, dpx2)) and
            (abs(px1 - nx2) < max(dpx1, dnx2))
        ):
            # Do not form a joint for anti-parallel S-walls.
            return []

        if (abs(nx2 - px1) < max(dnx2, dpx1)):
            ode_xs = nx1, px2
        elif (abs(px2 - nx1) < max(dpx2, dnx1)):
            ode_xs = px1, nx2

    elif g_type == 'D':
        # NOTE: For a D-type network, each S-wall is lifted
        # to a pair of trajectories on the spectral cover.
        m_nx1, m_nx2 = -nx2, -nx1

        if(
            ((abs(nx1 - px2) < max(dnx1, dpx2)) and
             (abs(px1 - nx2) < max(dpx1, dnx2))) or
            ((abs(m_nx1 - px2) < max(dnx2, dpx2)) and
             (abs(px1 - m_nx2) < max(dpx1, dnx1))) 
        ):
            # Do not form a joint for anti-parallel S-walls.
            return []
        if(
            (abs(nx2 - px1) < max(dnx2, dpx1)) and
            not ((abs(nx1 - (-px2)) < max(dnx1, dpx2)))
        ):
            ode_xs = nx1, px2
        elif(
            (abs(px2 - nx1) < max(dpx2, dnx1)) and
            not ((abs(px1 - (-nx2)) < max(dpx1, dnx2)))
        ):
            ode_xs = px1, nx2
        elif(
            (abs(m_nx2 - px1) < max(dnx2, dpx1)) and
            not ((abs(m_nx1 - (-px2)) < max(dnx2, dpx2)))
        ):
            ode_xs = m_nx1, px2
        elif(
            (abs(px2 - m_nx1) < max(dpx2, dnx1)) and
            not ((abs(px1 - (-m_nx2)) < max(dpx1, dnx1)))
        ):
            ode_xs = px1, m_nx2

    else:
        raise NotImplementedError

    if ode_xs is None:
        return []
    else:
        return [([], ode_xs)]


def get_joint_data_groups_by_roots(descendant_roots, z, sw_data, accuracy):
    # From a set of descendant roots find a set of joint data,
    #    [(root, ode_xs), ...].
    joint_data = []
    for root in descendant_roots:
        # TODO: The following, including self.ode_xs, can be removed
        # once the seeding of an S-wall is done by using a root.
        ffr_xs_at_z = sw_data.get_sheets_at_z(z, ffr=True).values()
        ffr_new_wall_weight_pairs = (
            sw_data.g_data.ordered_weight_pairs(root, ffr=True)
        )
        ffr_w_p_0 = ffr_new_wall_weight_pairs[0]
        ode_x1 = ffr_xs_at_z[ffr_w_p_0[0]]
        ode_x2 = ffr_xs_at_z[ffr_w_p_0[1]]
        ode_xs = [ode_x1, ode_x2]
        joint_data.append((root, ode_xs))

    # The following assumes that \lambda = x dz.
    # Group joint data according to x_1 - x_2,
    # which determines the trajectory of an S-wall.

    # joint_data_groups = [(roots, ode_xs), ...]
    joint_data_groups = []
    for root, ode_xs in joint_data:
        new_joint_data_group = True
        Dx = ode_xs[0] - ode_xs[1]
        for joint_data_group in joint_data_groups:
            group_roots, group_ode_xs = joint_data_group
            group_Dx = group_ode_xs[0] - group_ode_xs[1]
            if abs(Dx - group_Dx) < accuracy:
                new_joint_data_group = False
                group_roots.append(root)
                break
        if new_joint_data_group is True:
            joint_data_groups.append(([root], ode_xs))

    return joint_data_groups
