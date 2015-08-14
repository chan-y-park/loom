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

from geometry import RamificationPoint, SWData, find_xs_at_z_0
from s_wall import SWall, Joint, get_s_wall_seeds, get_joint
from misc import (n_nearest, n_nearest_indices, get_ode, left_right, clock,)
from intersection import (NoIntersection,
                          find_intersection_of_segments,
                          find_curve_range_intersection)
from scipy import interpolate


class SpectralNetwork:
    def __init__(
        self,
        phase=None,
    ):
        self.phase = phase
        self.s_walls = []
        self.joints = []

    def grow(self, config, sw_data):
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
        
        for bp in sw_data.branch_points:
            s_wall_seeds = get_s_wall_seeds(sw_data, self.phase, 
                                                        bp, config)
            for z_0, x_0, M_0 in s_wall_seeds:
                label = 'S-wall #{}'.format(len(self.s_walls))
                self.s_walls.append(
                    SWall(
                        z_0=z_0,
                        x_0=x_0,
                        M_0=M_0,
                        parents=[bp.label],
                        label=label,
                        n_steps=n_steps,
                    )
                )

        logging.info('Setup the ODE integrator...')
        ode = get_ode(sw_data, self.phase, accuracy)

        logging.info('Start growing a new spectral network...')
        ppzs = sw_data.punctures

        bpzs = [bp.z for bp in sw_data.branch_points]
        
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
                self.s_walls[i].grow(ode, bpzs, ppzs, config)
                self.determine_root_types(self.s_walls[i], sw_data)
                new_joints += self.get_new_joints(i, config, sw_data)

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
                if config['mass_limit'] is None or joint.M < config['mass_limit']:
                    self.s_walls.append(
                        SWall(
                            z_0=joint.z,
                            x_0=joint.x,
                            M_0=joint.M,
                            parents=joint.parents,
                            label=label,
                            n_steps=n_steps,
                        )
                    )
            iteration += 1

            logging.info('Iteration #{} finished.'.format(iteration))
#
#        ### Decorate S-walls with trivialization.
#        for k, s_wall in enumerate(self.s_walls):
#            logging.info('Adding trivialization info to S-wall #{}'
#                         .format(k))
#            ### Find weights corresponding to the x's of the endpoint
#            ### of the S-wall.
#            final_z = s_wall.z[-1]
#            sheet_xs = sw.get_sheet_xs_at_z(final_z)
#            sheets = []
#            ### Identify S-wall's x's with sheet #'s. 
#            s_wall_xs = s_wall.x[-1]
#            for x in s_wall_xs:
#                difference = numpy.fromiter(
#                    (abs(x - sheet_x) for sheet_x in sheet_xs),
#                    dtype=float,
#                )
#                i = difference.argsort()[0]
#                sheets.append(i)
#            ### Now sheets = [i_0, i_1], where s_wall.x[0][i_0] corresponds
#            ### to g_data.weights[i_0] and similarly for i_1.
#            s_wall.label += '\nweights = ('
#            for i in sheets:
#                s_wall.label += '{},'.format(sw.g_data.weights[i])
#            s_wall.label += ')'


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


    def get_new_joints(self, new_s_wall_index, config, sw_data):
        try:
            linux_distribution = platform.linux_distribution()[0]
            if linux_distribution != '':
                return self.get_new_joints_using_cgal(
                    new_s_wall_index, config,
                    sw_data, 
                    linux_distribution=linux_distribution,
                )
            else:
                raise OSError
        except OSError:
            logging.warning('CGAL not available; switch from '
                            'get_new_joints_using_cgal() to '
                            'get_new_joints_using_interpolation().')
            return self.get_new_joints_using_interpolation(new_s_wall_index,
                                                           config, sw_data)

    def get_new_joints_using_cgal(self, new_s_wall_index, config,
                                  sw_data, linux_distribution=None):
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
                    M_n = new_s_wall.M[t_n]

                    # t_p: index of z_seg_p nearest to ip_z
                    t_p = n_nearest_indices(prev_s_wall.z, ip_z, 1)[0]
                    M_p = prev_s_wall.M[t_p]

                    # TODO: need to put the joint into the parent
                    # S-walls?

                    # find the values of x at z = ip_z.
                    ip_xs = find_xs_at_z_0(sw_data, ip_z)
                    ip_x_n_0 = n_nearest(ip_xs, x_n[0], 1)[0]
                    ip_x_n_1 = n_nearest(ip_xs, x_n[1], 1)[0]
                    ip_x_p_0 = n_nearest(ip_xs, x_p[0], 1)[0]
                    ip_x_p_1 = n_nearest(ip_xs, x_p[1], 1)[0]

                    a_joint = get_joint(
                        ip_z, 
                        prev_s_wall,
                        new_s_wall,                         
                        t_p, 
                        t_n,
                        sw_data=sw_data.sw_data,
                    )

                    if(a_joint is None):
                        continue
                    else:
                        new_joints.append(a_joint)

            except NoIntersection:
                pass

        return new_joints


    def get_new_joints_using_interpolation(
                                            self, new_s_wall_index, 
                                            config, sw_data
                                        ):
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
        # Use CGAL to find multiple z-intersections.
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
                        ip_xs = find_xs_at_z_0(sw_data, ip_z)
                        ip_x_n_0 = n_nearest(ip_xs, x_n[0], 1)[0]
                        ip_x_n_1 = n_nearest(ip_xs, x_n[1], 1)[0]
                        ip_x_p_0 = n_nearest(ip_xs, x_p[0], 1)[0]
                        ip_x_p_1 = n_nearest(ip_xs, x_p[1], 1)[0]

                        # find mass of parent S-walls: this is approximate,
                        # since we don't interpolate precisely to the joint
                        # TODO: improve by doing precise interpolation
                        M_n = new_s_wall.M[t_n]
                        M_p = prev_s_wall.M[t_p]
                        a_joint = get_joint(
                            ip_z, 
                            prev_s_wall,
                            new_s_wall, 
                            t_p,
                            t_n,
                            sw_data=sw_data
                        )

                        if(a_joint is None):
                            continue
                        else:
                            new_joints.append(a_joint)

                    except NoIntersection:
                        pass
        return new_joints

    def determine_root_types(self, s_wall, sw_data):
        """
        Determine at which points the wall crosses a cut, 
        for instance [55, 107, 231] would mean that 
        it changes root-type 3 times. 
        Hence, s_wall.splittings would have length 3, 
        while s_wall.local_roots would have length 4.
        """
        # Determine the initial root-type
        z_0 = s_wall.z[0]
        x_0 = s_wall.x[0]
        initial_root = get_s_wall_root(z_0, x_0, sw_data,)
        # A list of ordered pairs [...[i, j]...]
        # such that weights[j] - weights[i] = root
        initial_weight_pairs = (
                    sw_data.g_data.ordered_weight_pairs(initial_root,)
                    )
        
        s_wall.local_roots = [initial_root]
        s_wall.local_weight_pairs = [initial_weight_pairs]

        # If the length of the S-wall's coordinates
        # is 2 or less, do not check cuts.
        # Otherwise, the interpolation methood would
        # raise an error.
        if len(s_wall.z) < 3:
            return None

        branch_points = sw_data.branch_points
        bpzs_r = [bp.z.real for bp in branch_points]
        
        # parametrizing the z-coordinate of the k-wall's coordinates
        # as a function of proper time
        traj_t = numpy.array(range(len(s_wall.z)))
        traj_x = numpy.array([w.real for w in s_wall.z])
        
        # Scan over branch cuts, see if path ever crosses one 
        # based on x-coordinates only
        for b_pt_num, x_0 in list(enumerate(bpzs_r)):
            g = interpolate.splrep(traj_t, traj_x - x_0, s=0)
            # now produce a list of integers corresponding to points in the 
            # S-wall's coordinate list that seem to cross branch-cuts
            # based on the z-coordinate's real part.
            # Will get a list [i_0, i_1, ...] of intersections
            intersections = map(int, map(round, interpolate.sproot(g)))
            # removing duplicates
            intersections = list(set(intersections))
            # enforcing imaginary-part of z-coordinate intersection criterion:
            # branch cuts extend vertically
            y_0 = branch_points[b_pt_num].z.imag
            intersections = (
                    [i for i in intersections if s_wall.z[i].imag > y_0 ]
                    )
            # adding the branch-point identifier to each intersection
            intersections = (
                [[branch_points[b_pt_num], i] for i in intersections]
                )
            # dropping intersections of a primary S-wall with the 
            # branch cut emanating from its parent branch-point
            # if such intersections happens at t=0 or t=1
            intersections = (
                    [[bp, i] for bp, i in intersections if (
                    not (bp.label == s_wall.parents[0] and (i == 0 or i==1))
                    )]
                )
            # add the direction to the intersection data: either 'cw' or 'ccw'
            intersections = ([
                        [bp, i, clock(left_right(s_wall.z, i))] 
                        for bp, i in intersections
                    ])

            s_wall.cuts_intersections += intersections
        # Might be worth implementing an algorithm for handling 
        # overlapping branch cuts: e.g. the one with a lower starting point 
        # will be taken to be on the left, or a similar criterion.
        # Still, there will be other sorts of problems, it is necessary
        # to just rotate the z-plane and avoid such situations.

        # now sort intersections according to where they happen in proper 
        # time; recall that the elements of cuts_intersections are organized 
        # as      [..., [branch_point, t, 'ccw'] ,...]
        # where 't' is the integer of proper time at the intersection.
        s_wall.cuts_intersections = sorted(
                                    s_wall.cuts_intersections, 
                                    cmp = lambda k1,k2: cmp(k1[1],k2[1])
                                    )

        logging.debug(
            '\nS-wall {}\nintersects the following'
            'cuts at the points\n{}\n'.format(s_wall.label, intersections))

        # now define the list of splitting points (for convenience) ad the 
        # list of local charges
        s_wall.splittings = [t for bp, t, chi in s_wall.cuts_intersections]

        for k in range(len(s_wall.cuts_intersections)):
            branch_point = s_wall.cuts_intersections[k][0]   # branch-point
            # t = s_wall.cuts_intersections[k][1]           # proper time
            direction = s_wall.cuts_intersections[k][2]     # 'ccw' or 'cw'
            current_root = s_wall.local_roots[-1]
            new_root = sw_data.g_data.weyl_monodromy(
                                    current_root, branch_point, direction)
            new_weight_pairs = sw_data.g_data.ordered_weight_pairs(new_root)
            s_wall.local_roots.append(new_root)
            s_wall.local_weight_pairs.append(new_weight_pairs)

def get_s_wall_root(z, xs, sw_data):
    x_i, x_j = xs

    # The following is a dictionary
    sheets_at_z = sw_data.get_sheets_at_z(z)
    xs_at_z = sheets_at_z.values()
    
    # Sheet matching x_i
    closest_to_x_i = sorted(xs_at_z, key=lambda x: abs(x - x_i))[0]
    i = [k for k, v in sheets_at_z.iteritems() if v == closest_to_x_i][0]

    # Sheet matching x_j
    closest_to_x_j = sorted(xs_at_z, key=lambda x: abs(x - x_j))[0]
    j = [k for k, v in sheets_at_z.iteritems() if v == closest_to_x_j][0]

    return sw_data.g_data.weights[j] - sw_data.g_data.weights[i]
