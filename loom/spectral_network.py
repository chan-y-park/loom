import platform
import os
import numpy
import ctypes
import logging
import signal
import multiprocessing
import subprocess
import json

import pdb

from cmath import exp
from itertools import combinations

from geometry import RamificationPoint, SWData
from s_wall import SWall, Joint, get_s_wall_seeds, get_joint
from misc import (n_nearest, n_nearest_indices, find_xs_at_z_0, get_ode)
from hit_table import HitTable
from intersection import (NoIntersection,
                          find_intersection_of_segments, 
                          find_curve_range_intersection)


class SpectralNetwork:
    def __init__(
        self,
        phase=None,
        sw_data=None,
    ):
        self.phase = phase
        self.s_walls = []
        self.joints = []
        self.sw_data =sw_data


    def grow(self, config):
        """
        Grow the spectral network by seeding SWall's
        and then calling SWall.grow() for each S-wall.

        Finding joints and growing S-walls from the joints
        is done for a given iteration, therefore it is possible
        that there is a joint from which an S-wall is not grown
        if the depth of the joint is too deep.
        """
        accuracy = config['accuracy']
        n_steps=config['num_of_steps']
        logging.info('Start growing a new spectral network...')

        logging.info('Seed S-walls at ramification points...')
        
        for bp in self.sw_data.branch_points:
            s_wall_seeds = get_s_wall_seeds(self.sw_data, self.phase, 
                                                        bp, config)
            for z_0, x_0 in s_wall_seeds:
                label = 'S-wall #{}'.format(len(self.s_walls))
                self.s_walls.append(
                    SWall(
                        z_0=z_0,
                        x_0=x_0,
                        parents=[bp.label],
                        label=label,
                        n_steps=n_steps,
                        spectral_network=self,
                    )
                )

        logging.info('Setup the ODE integrator...')
        ode = get_ode(self.sw_data, self.phase, accuracy)

        logging.info('Start growing a new spectral network...')
        ppzs = self.sw_data.punctures

        bpzs = [bp.z for bp in self.sw_data.branch_points]
        
        n_finished_s_walls = 0 
        iteration = 0
        while(iteration < config['num_of_iterations']):
            """
            Iterate until there is no new joint.
            Each S-wall is grown only once.
            """
            new_joints = []     # number of new joints found in each iteration
            for i in range(n_finished_s_walls, len(self.s_walls)):
                logging.info('Growing S-wall #{}...'.format(i))
                self.s_walls[i].grow(ode, bpzs, ppzs, config,)
                new_joints += self.get_new_joints(i, config)

            n_finished_s_walls = len(self.s_walls)
            if(len(new_joints) == 0):
                logging.info('No additional joint found: '
                             'Stop growing this spectral network '
                             'at iteration #{}.'.format(iteration))
                break
            else:
                logging.info('Growing S-walls in iteration #{} finished; '
                             'start finding joints of S-walls.'
                             .format(iteration))

            # Seed an S-wall for each new joint.
            for joint in new_joints:
                joint_is_new = True
                # check if this is not in the previous joint list
                for old_joint in self.joints:
                    if joint.is_equal_to(old_joint, accuracy):
                        joint_is_new = False
                        break
                if not joint_is_new:
                    continue
                joint.label = 'joint #{}'.format(len(self.joints))
                self.joints.append(joint)
                label = 'S-wall #{}'.format(len(self.s_walls))
                self.s_walls.append(
                    SWall(
                        z_0=joint.z,
                        x_0=joint.x, 
                        parents=joint.parents,
                        label=label,
                        n_steps=n_steps,
                        spectral_network=self,
                    )
                )
            iteration += 1

            logging.info('Iteration #{} finished.'.format(iteration))


    def save_json_data(self, file_object, **kwargs):
        """
        Save the spectral network data in a JSON-compatible file.
        """
        json_data = {}

        json_data['phase'] = self.phase
        json_data['s_walls'] = [s_wall.get_json_data()
                                for s_wall in self.s_walls]
        json_data['joints'] = [joint.get_json_data()
                               for joint in self.joints]
        json.dump(json_data, file_object, **kwargs)
    

    def set_from_json_data(self, file_object, **kwargs):
        """
        Load the spectral network data from a JSON-compatible file.
        """
        json_data = json.load(file_object, **kwargs)

        self.phase = json_data['phase']

        for s_wall_data in json_data['s_walls']:
            an_s_wall = SWall()
            an_s_wall.set_json_data(s_wall_data)
            self.s_walls.append(an_s_wall)

        for joint_data in json_data['joints']:
            a_joint = Joint()
            a_joint.set_json_data(joint_data)
            self.joints.append(a_joint)


    def get_new_joints(self, new_s_wall_index, config):
        try:
            linux_distribution = platform.linux_distribution()[0]
            if linux_distribution != '':
                return self.get_new_joints_using_cgal(
                    new_s_wall_index, config,
                    linux_distribution=linux_distribution,
                )
            else:
                raise OSError
        except OSError:
            logging.warning('CGAL not available; switch from '
                            'get_new_joints_using_cgal() to '
                            'get_new_joints_using_interpolation().')
            return self.get_new_joints_using_interpolation(new_s_wall_index,
                                                           config,)

    def get_new_joints_using_cgal(self, new_s_wall_index, config,
                                  linux_distribution=None):
        """
        Find new wall-wall intersections using CGAL 2d curve intersection.
        """
        new_joints = []
        if (config['root_system'] in ['A1',]):
            logging.info('There is no joint for the given root system {}.'
                         .format(config['root_system']))
            return new_joints
        lib_name = 'libcgal_intersection'
        if linux_distribution == 'Ubuntu':
            lib_name += '_ubuntu'
        elif linux_distribution == 'Scientific Linux':
            lib_name += '_het-math2'
        else:
            raise OSError

        logging.info('Using CGAL to find intersections.')

        # Load CGAL shared library.
        libcgal_intersection = numpy.ctypeslib.load_library(
            lib_name, 
            os.path.dirname(os.path.realpath(__file__)) + '/cgal_intersection/'
        )
        cgal_find_intersections_of_curves = (libcgal_intersection.
                                             find_intersections_of_curves)
        # Prepare types for CGAL library.
        array_2d_float = numpy.ctypeslib.ndpointer(
            dtype=numpy.float64,
            ndim=2,
            flags=['C_CONTIGUOUS', 'ALIGNED'],
        )
        array_2d_complex = numpy.ctypeslib.ndpointer(
            dtype=numpy.complex128,
            ndim=1,
            flags=['C_CONTIGUOUS', 'ALIGNED'],
        )

        cgal_find_intersections_of_curves.restype = ctypes.c_int
        cgal_find_intersections_of_curves.argtypes = [
            #array_2d_float,
            array_2d_complex,
            ctypes.c_long,
            #array_2d_float, 
            array_2d_complex,
            ctypes.c_long,
            array_2d_float, ctypes.c_int,
        ]

        new_s_wall = self.s_walls[new_s_wall_index]

        for prev_s_wall in self.s_walls[:new_s_wall_index]:

            if prev_s_wall.label in new_s_wall.parents:
                continue
                    
            # Check if the two S-walls have a common x-range.  
            have_common_x_range = False
            for x_a, x_b in (
                (x_a, x_b) for x_a in new_s_wall.x.T 
                for x_b in prev_s_wall.x.T
            ):
                x_r_range, x_i_range = (
                    find_curve_range_intersection(
                        (x_a.real, x_a.imag),
                        (x_b.real, x_b.imag),
                        cut_at_inflection=False,
                    )
                )
                if ((x_r_range.is_EmptySet is True) or 
                    (x_r_range.is_EmptySet is True) or 
                    x_i_range.is_FiniteSet or
                    x_i_range.is_FiniteSet):
                    continue
                else:
                    have_common_x_range = True
                    break
            if have_common_x_range is False:
                # No common x range, therefore no joint.
                continue

            # Find an intersection on the z-plane.
            try:
                buffer_size = 10
                intersection_search_finished = False
                while not intersection_search_finished:
                    intersections = numpy.empty((buffer_size, 2), 
                                                dtype=numpy.float64)
                    #from plotting import plot_s_walls
                    #plot_s_walls([new_s_wall, prev_s_wall])
                    num_of_intersections = cgal_find_intersections_of_curves(
                        #new_s_wall.z.view(numpy.dtype((numpy.float64,2))),
                        new_s_wall.z,
                        ctypes.c_long(len(new_s_wall.z)), 
                        prev_s_wall.z, 
                        ctypes.c_long(len(prev_s_wall.z)), 
                        intersections, buffer_size
                    )
                    if num_of_intersections == 0:
                        intersection_search_finished = True
                        raise NoIntersection
                    elif num_of_intersections > buffer_size:
                        logging.info('Number of intersections larger than '
                                     'the buffer size; increase its size '
                                     'and run intersection finding again.')
                        buffer_size = num_of_intersections
                    else:
                        intersections.resize((num_of_intersections,2))
                        intersection_search_finished = True

                for ip_x, ip_y in intersections:
                    ip_z = ip_x + 1j*ip_y

                    # t_n: index of new_s_wall.z nearest to ip_z
                    t_n = n_nearest_indices(new_s_wall.z, ip_z, 1)[0]
                    x_n = new_s_wall.x[t_n]

                    # t_p: index of z_seg_p nearest to ip_z
                    t_p = n_nearest_indices(prev_s_wall.z, ip_z, 1)[0]
                    x_p = prev_s_wall.x[t_p]

                    # TODO: need to put the joint into the parent
                    # S-walls?

                    # find the values of x at z = ip_z.
                    ip_xs = find_xs_at_z_0(self.sw_data, ip_z)
                    ip_x_n_0 = n_nearest(ip_xs, x_n[0], 1)[0]
                    ip_x_n_1 = n_nearest(ip_xs, x_n[1], 1)[0]
                    ip_x_p_0 = n_nearest(ip_xs, x_p[0], 1)[0]
                    ip_x_p_1 = n_nearest(ip_xs, x_p[1], 1)[0]

                    a_joint = get_joint(
                        ip_z, ip_x_n_0, ip_x_n_1, ip_x_p_0, ip_x_p_1,
                        new_s_wall.label, 
                        prev_s_wall.label,
                        accuracy=config['accuracy'],
                        xs_at_z=ip_xs,
                        g_data=self.sw_data.g_data,
                    )

                    if(a_joint is None):
                        continue
                    else:
                        new_joints.append(a_joint)

            except NoIntersection:
                pass

        return new_joints


    def get_new_joints_using_interpolation(self, new_s_wall_index, config):
        """
        Find joints between the newly grown segment of the given S-wall
        and the other S-walls by interpolating S-walls with functions and
        finding roots of the pairwise differences of the functions.
        
        This checks joints that are formed by two
        S-walls only, not considering the possibility of a joint of three
        S-walls, which in principle can happen but is unlikely in a numerical
        setup.
        """
        new_joints = []
        if (config['root_system'] in ['A1',]):
            logging.info('There is no joint for the given root system {}.'
                         .format(config['root_system']))
            return new_joints

        new_s_wall = self.s_walls[new_s_wall_index]
        new_tps = new_s_wall.get_turning_points()
        new_z_segs = numpy.split(new_s_wall.z, new_tps, axis=0,)
        new_x_segs = numpy.split(new_s_wall.x, new_tps, axis=0,)

        # NOTE: Here we find only a single joint between two S-walls.
        # If needed, change the part of getting the z-intersection
        # such that it finds multiple z-intersections, or use HitTable.
        for prev_s_wall in self.s_walls[:new_s_wall_index]:
            prev_tps = prev_s_wall.get_turning_points()
            prev_z_segs = numpy.split(prev_s_wall.z, prev_tps, axis=0,)
            prev_x_segs = numpy.split(prev_s_wall.x, prev_tps, axis=0,)

            for i_n in range(len(new_tps)+1):
                z_seg_n = new_z_segs[i_n]
                for i_p in range(len(prev_tps)+1):
                    z_seg_p = prev_z_segs[i_p]

                    # Check if the two segments have a common x-range.  
                    have_common_x_range = False
                    for x_a, x_b in (
                        (x_a, x_b) for x_a in new_x_segs[i_n].T 
                        for x_b in prev_x_segs[i_p].T
                    ):
                        x_r_range, x_i_range = (
                            find_curve_range_intersection(
                                (x_a.real, x_a.imag),
                                (x_b.real, x_b.imag),
                                cut_at_inflection=False,
                            )
                        )
                        if ((x_r_range.is_EmptySet is True) or 
                            (x_r_range.is_EmptySet is True) or 
                            x_i_range.is_FiniteSet or
                            x_i_range.is_FiniteSet):
                            continue
                        else:
                            have_common_x_range = True
                            break
                    if have_common_x_range is False:
                        # No common x range, therefore no joint.
                        continue

                    # Find an intersection on the z-plane.
                    try:
                        ip_x, ip_y = find_intersection_of_segments(
                            (z_seg_n.real, z_seg_n.imag), 
                            (z_seg_p.real, z_seg_p.imag), 
                            config['accuracy'],
                        )
                        ip_z = ip_x + 1j*ip_y

                        # t_n: index of z_seg_n nearest to ip_z
                        t_n = n_nearest_indices(z_seg_n, ip_z, 1)[0]
                        x_n = new_x_segs[i_n][t_n]

                        # t_p: index of z_seg_p nearest to ip_z
                        t_p = n_nearest_indices(z_seg_p, ip_z, 1)[0]
                        x_p = prev_x_segs[i_p][t_p]

                        # TODO: need to put the joint into the parent
                        # S-walls?

                        # find the values of x at z = ip_z.
                        ip_xs = find_xs_at_z_0(self.sw_data, ip_z)
                        ip_x_n_0 = n_nearest(ip_xs, x_n[0], 1)[0]
                        ip_x_n_1 = n_nearest(ip_xs, x_n[1], 1)[0]
                        ip_x_p_0 = n_nearest(ip_xs, x_p[0], 1)[0]
                        ip_x_p_1 = n_nearest(ip_xs, x_p[1], 1)[0]

                        a_joint = get_joint(
                            ip_z, ip_x_n_0, ip_x_n_1, ip_x_p_0, ip_x_p_1,
                            new_s_wall.label, 
                            prev_s_wall.label,
                            accuracy=config['accuracy'],
                            xs_at_z=ip_xs,
                            g_data=self.sw_data.g_data,
                        )

                        if(a_joint is None):
                            continue
                        else:
                            new_joints.append(a_joint)

                    except NoIntersection:
                        pass
        return new_joints
