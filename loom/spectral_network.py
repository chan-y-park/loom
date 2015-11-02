import platform
import os
import numpy
import sympy
import ctypes
import logging
import json
# import pdb

from cmath import exp
from scipy import integrate
from s_wall import SWall, Joint, get_s_wall_seeds, get_joint
from misc import (
    n_nearest_indices, 
)
from intersection import (
    NoIntersection, find_intersection_of_segments,
)

CLIPPING_RADIUS = 30.0

class SpectralNetwork:
    def __init__(
        self,
        phase=None,
        logger_name='loom',
    ):
        self.phase = phase
        self.s_walls = []
        self.joints = []
        self.logger_name = logger_name

    def grow(self, config, sw_data):
        """
        Grow the spectral network by seeding SWall's
        and then calling SWall.grow() for each S-wall.

        Finding joints and growing S-walls from the joints
        is done for a given iteration, therefore it is possible
        that there is a joint from which an S-wall is not grown
        if the depth of the joint is too deep.
        """
        logger = logging.getLogger(self.logger_name)
        
        # Add a signal handler

        accuracy = config['accuracy']
        num_of_iterations = config['num_of_iterations']
        n_steps = config['num_of_steps']
        logger.info('Start growing a new spectral network...')

        logger.info('Seed S-walls at branch points...')
        
        for bp in sw_data.branch_points:
            s_wall_seeds = get_s_wall_seeds(
                sw_data, self.phase, bp, config, self.logger_name,
            )
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
                        logger_name=self.logger_name,
                    )
                )

        logger.debug('Setup the ODE integrator...')
        ode = get_ode(sw_data, self.phase, accuracy)
        punctures = sw_data.regular_punctures + sw_data.irregular_punctures
        ppzs = [p.z for p in punctures]

        bpzs = [bp.z for bp in sw_data.branch_points]
        
        n_finished_s_walls = 0 
        iteration = 1
        while(iteration <= num_of_iterations):
            """
            Iterate until there is no new joint
            or for a specified number of iterations.
            Each S-wall is grown only once.
            """
            new_joints = []     # number of new joints found in each iteration
            for i in range(n_finished_s_walls, len(self.s_walls)):
                s_i = self.s_walls[i]
                logger.info('Growing S-wall #{}...'.format(i))
                s_i.grow(
                    ode, bpzs, ppzs, config, clipping_radius=CLIPPING_RADIUS
                )
                # Cut the grown S-walls at the intersetions with branch cuts
                # and decorate each segment with its root data.
                logger.info('Determining the root type of S-wall #{}...'
                            .format(i))
                s_i.determine_root_types(sw_data)
                new_joints += self.get_new_joints(i, config, sw_data)

            n_finished_s_walls = len(self.s_walls)
            if(len(new_joints) == 0):
                logger.info('No additional joint found: '
                                 'Stop growing this spectral network '
                                 'at iteration #{}.'.format(iteration))
                break
            elif iteration == num_of_iterations:
                # Last iteration finished. Do not form new joints,
                # and do not seed additional S-walls, either.
                break
            else:
                logger.info('Growing S-walls in iteration #{} finished.'
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
                if (
                    config['mass_limit'] is None or 
                    joint.M < config['mass_limit']
                ):
                    self.s_walls.append(
                        SWall(
                            z_0=joint.z,
                            # The numerics of S-walls involves sheets
                            # from the first fundamental cover.
                            x_0=joint.ode_xs,
                            M_0=joint.M,
                            parents=joint.parents,
                            label=label,
                            n_steps=n_steps,
                            logger_name=self.logger_name,
                        )
                    )
            logger.info('Iteration #{} finished.'.format(iteration))
            iteration += 1

    def get_json_data(self):
        """
        Prepare the spectral network data in a JSON-compatible file.
        """
        json_data = {}
        json_data['phase'] = self.phase
        json_data['s_walls'] = [s_wall.get_json_data()
                                for s_wall in self.s_walls]
        json_data['joints'] = [joint.get_json_data()
                               for joint in self.joints]
        return json_data

    def set_from_json_data(self, json_data, sw_data):
        """
        Load the spectral network data from a JSON-compatible file.
        """
        branch_points = sw_data.branch_points

        self.phase = json_data['phase']

        for s_wall_data in json_data['s_walls']:
            an_s_wall = SWall(logger_name=self.logger_name,)
            an_s_wall.set_from_json_data(s_wall_data, branch_points)
            self.s_walls.append(an_s_wall)

        for joint_data in json_data['joints']:
            a_joint = Joint()
            a_joint.set_from_json_data(joint_data)
            self.joints.append(a_joint)

    def get_new_joints(self, new_s_wall_index, config, sw_data):
        logger = logging.getLogger(self.logger_name)
        if (config['root_system'] in ['A1', ]):
            logger.info('There is no joint for the given root system {}.'
                             .format(config['root_system']))
            return [] 
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
            logger.warning('CGAL not available; switch from '
                                'get_new_joints_using_cgal() to '
                                'get_new_joints_using_interpolation().')
            return self.get_new_joints_using_interpolation(new_s_wall_index,
                                                           config, sw_data)

    def get_new_joints_using_cgal(self, new_s_wall_index, config,
                                  sw_data, linux_distribution=None):
        """
        Find new wall-wall intersections using CGAL 2d curve intersection.
        """
        logger = logging.getLogger(self.logger_name)
        new_joints = []
#        if (config['root_system'] in ['A1', ]):
#            logger.info('There is no joint for the given root system {}.'
#                             .format(config['root_system']))
#            return new_joints
        lib_name = 'libcgal_intersection'
        if linux_distribution == 'Ubuntu':
            lib_name += '_ubuntu'
        elif linux_distribution == 'debian':
            # FIXME: Anaconda Python returns 'debian' instead of 'Ubuntu'.
            lib_name += '_ubuntu'
        elif linux_distribution == 'Scientific Linux':
            lib_name += '_het-math2'
        else:
            raise OSError

        logger.info('Using CGAL to find intersections.')

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
            # array_2d_float,
            array_2d_complex,
            ctypes.c_long,
            # array_2d_float, 
            array_2d_complex,
            ctypes.c_long,
            array_2d_float, ctypes.c_int,
        ]

        new_s_wall = self.s_walls[new_s_wall_index]

        for prev_s_wall in self.s_walls[:new_s_wall_index]:

            if prev_s_wall.label in new_s_wall.parents:
                continue
            # Find an intersection on the z-plane.
            try:
                buffer_size = 10
                intersection_search_finished = False
                while not intersection_search_finished:
                    intersections = numpy.empty((buffer_size, 2), 
                                                dtype=numpy.float64)
                    num_of_intersections = cgal_find_intersections_of_curves(
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
                        logger.info('Number of intersections larger than '
                                         'the buffer size; increase its size '
                                         'and run intersection finding again.')
                        buffer_size = num_of_intersections
                    else:
                        intersections.resize((num_of_intersections, 2))
                        intersection_search_finished = True

                for ip_x, ip_y in intersections:
                    ip_z = ip_x + 1j * ip_y

                    # t_n: index of new_s_wall.z nearest to ip_z
                    t_n = n_nearest_indices(new_s_wall.z, ip_z, 1)[0]

                    # t_p: index of z_seg_p nearest to ip_z
                    t_p = n_nearest_indices(prev_s_wall.z, ip_z, 1)[0]

                    # TODO: need to put the joint into the parent
                    # S-walls?

                    logger.debug(
                        'evaluating possible joint at z = {}'.format(ip_z)
                    )
                    a_joint = get_joint(
                        ip_z, 
                        prev_s_wall,
                        new_s_wall,                         
                        t_p, 
                        t_n,
                        sw_data,
                    )

                    if(a_joint is None):
                        continue
                    else:
                        new_joints.append(a_joint)

            except NoIntersection:
                pass

        return new_joints

    def get_new_joints_using_interpolation(
        self, new_s_wall_index, config, sw_data, 
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
        logger = logging.getLogger(self.logger_name)
        new_joints = []
#        if (config['root_system'] in ['A1', ]):
#            logger.info('There is no joint for the given root system {}.'
#                             .format(config['root_system']))
#            return new_joints

        new_s_wall = self.s_walls[new_s_wall_index]
        new_tps = new_s_wall.get_turning_points()
        new_z_segs = numpy.split(new_s_wall.z, new_tps, axis=0,)

        # NOTE: Here we find only a single joint between two S-walls.
        # Use CGAL to find multiple z-intersections.
        for prev_s_wall in self.s_walls[:new_s_wall_index]:
            prev_tps = prev_s_wall.get_turning_points()
            prev_z_segs = numpy.split(prev_s_wall.z, prev_tps, axis=0,)

            for i_n in range(len(new_tps) + 1):
                z_seg_n = new_z_segs[i_n]
                for i_p in range(len(prev_tps) + 1):
                    z_seg_p = prev_z_segs[i_p]
                    # Find an intersection on the z-plane.
                    try:
                        ip_x, ip_y = find_intersection_of_segments(
                            (z_seg_n.real, z_seg_n.imag),
                            (z_seg_p.real, z_seg_p.imag),
                            config['accuracy'],
                        )
                        ip_z = ip_x + 1j * ip_y

                        # t_n: index of z_seg_n nearest to ip_z
                        t_n = n_nearest_indices(new_s_wall.z, ip_z, 1)[0]

                        # t_p: index of z_seg_p nearest to ip_z
                        t_p = n_nearest_indices(prev_s_wall.z, ip_z, 1)[0]

                        # TODO: need to put the joint into the parent
                        # S-walls?

                        # find mass of parent S-walls: this is approximate,
                        # since we don't interpolate precisely to the joint
                        # TODO: improve by doing precise interpolation
                        logger.debug(
                            'evaluating possible joint at z = {}'.format(ip_z)
                        )
                        a_joint = get_joint(
                            ip_z, 
                            prev_s_wall,
                            new_s_wall, 
                            t_p,
                            t_n,
                            sw_data,
                        )

                        if(a_joint is None):
                            continue
                        else:
                            new_joints.append(a_joint)

                    except NoIntersection:
                        pass
        return new_joints


def get_ode(sw, phase, accuracy):
    x, z = sympy.symbols('x z')
    ode_absolute_tolerance = accuracy

    # Even for higher-reps, we always use the 
    # first fundamental representation curve 
    # for evolving the network
    f = sw.ffr_curve.num_eq
    df_dz = f.diff(z)
    df_dx = f.diff(x)
    # F = -(\partial f / \partial z) / (\partial f / \partial x).
    F = sympy.lambdify((z, x), sympy.simplify(-df_dz / df_dx))
    # Do we need to call sw.diff? Or can we just use the values of x_i's?
    v = sympy.lambdify((z, x), sw.diff.num_v)

    def ode_f(t, zx1x2M):
        z_i = zx1x2M[0]
        x1_i = zx1x2M[1]
        x2_i = zx1x2M[2]
        dz_i_dt = exp(phase * 1j) / (v(z_i, x1_i) - v(z_i, x2_i))
        dx1_i_dt = F(z_i, x1_i) * dz_i_dt
        dx2_i_dt = F(z_i, x2_i) * dz_i_dt
        dM_dt = 1
        return [dz_i_dt, dx1_i_dt, dx2_i_dt, dM_dt]

    ode = integrate.ode(ode_f)
    ode.set_integrator(
        'zvode',
        # method='adams',
        atol=ode_absolute_tolerance,
    )

    return ode
