import platform
import os
import numpy
import sympy
import ctypes
import logging
import json
import pdb
import itertools

from cmath import exp
from scipy import integrate
from s_wall import SWall, Joint, get_s_wall_seeds
from misc import (
    n_nearest_indices, is_root,
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
        
        accuracy = config['accuracy']
        num_of_iterations = config['num_of_iterations']
        n_steps = config['num_of_steps']

        # Determine the intersection-finding algorithm
        # according to the availability of CGAL.
        try:
            lib_name = 'libcgal_intersection'

            linux_distribution = platform.linux_distribution()[0]
            if linux_distribution == 'Ubuntu':
                lib_name += '_ubuntu'
            elif linux_distribution == 'debian':
                # FIXME: Anaconda Python returns 'debian' instead of 'Ubuntu'.
                lib_name += '_ubuntu'
            elif linux_distribution == 'Scientific Linux':
                lib_name += '_het-math2'
            else:
                raise OSError

            # Load CGAL shared library.
            libcgal_intersection = numpy.ctypeslib.load_library(
                lib_name, 
                (os.path.dirname(os.path.realpath(__file__)) +
                 '/cgal_intersection/'),
            )
            
            self.get_intersections = (libcgal_intersection.
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

            self.get_intersections.restype = ctypes.c_int
            self.get_intersections.argtypes = [
                # array_2d_float,
                array_2d_complex,
                ctypes.c_long,
                # array_2d_float, 
                array_2d_complex,
                ctypes.c_long,
                array_2d_float, ctypes.c_int,
            ]

            logger.info('Use CGAL to find intersections.')
            self.use_cgal = True

        except OSError:
            logger.warning('CGAL not available; use interpolation '
                           'to find intersections.')
            self.get_intersections = find_intersections_of_curves
            self.use_cgal = False

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
                    ode, bpzs, ppzs, config,
                    clipping_radius=CLIPPING_RADIUS,
                )
                # Cut the grown S-walls at the intersetions with branch cuts
                # and decorate each segment with its root data.
                logger.info('Determining the root type of S-wall #{}...'
                            .format(i))
                s_i.determine_root_types(sw_data)

                logger.info('Finding new joints from S-wall #{}...'.format(i))
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

        # Remove the function pointer that cannot be pickled.
        self.get_intersections = None
        self.use_cgal = None

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
        """
        Find intersections between the new S-wall and 
        previoulsy existing S-walls using 
        either CGAL 2d curve intersection or scipy interpolation
        according to the availability, then form joints
        from the intersection points.
        """

        logger = logging.getLogger(self.logger_name)
        accuracy = config['accuracy']
        new_joints = []

        if (config['root_system'] in ['A1', ]):
            logger.info('There is no joint for the given root system {}.'
                             .format(config['root_system']))
            return [] 

        new_s_wall = self.s_walls[new_s_wall_index]

        for prev_s_wall in self.s_walls[:new_s_wall_index]:

            # First check if the two S-walls are compatible
            # for forming a joint.

            # 1. Check if the new S-wall is a descendant
            # of an existing S-wall. 
            if prev_s_wall.label in new_s_wall.parents:
                continue
            
            # 2. Split the two S-walls into segments 
            # according to the trivialization, then
            # check the compatibility of a pair
            # of segments.
            new_z_segs = new_s_wall.get_z_segs()
            num_new_z_segs = len(new_z_segs)
            prev_z_segs = prev_s_wall.get_z_segs()
            num_prev_z_segs = len(prev_z_segs)
                
            for n_i, p_i in itertools.product(
                range(num_new_z_segs), range(num_prev_z_segs)
            ):
                n_seg_zs = new_z_segs[n_i]
                n_seg_root = new_s_wall.local_roots[n_i]
                p_seg_zs = prev_z_segs[p_i]
                p_seg_root = prev_s_wall.local_roots[p_i]

                if not is_root(p_seg_root + n_seg_root, sw_data.g_data):
                    # The two segments are not compatible for
                    # forming a joint.
                    continue

                # Find an intersection of the two segments on the z-plane.
                if self.use_cgal is True:
                    buffer_size = 10
                    intersection_search_finished = False

                    while not intersection_search_finished:
                        intersections = numpy.empty((buffer_size, 2), 
                                                    dtype=numpy.float64)
                        num_of_intersections = self.get_intersections(
                            n_seg_zs,
                            ctypes.c_long(len(n_seg_zs)), 
                            p_seg_zs,
                            ctypes.c_long(len(p_seg_zs)), 
                            intersections, buffer_size
                        )
                        if num_of_intersections == 0:
                            intersection_search_finished = True
                        elif num_of_intersections > buffer_size:
                            logger.info('Number of intersections larger than '
                                        'the buffer size; increase its size '
                                        'to {} and find intersections again.'
                                        .format(num_of_intersections))
                            buffer_size = num_of_intersections
                        else:
                            intersection_search_finished = True

                        intersections.resize((num_of_intersections, 2))
                else:
                    intersections = self.get_intersections(
                        n_seg_zs, p_seg_zs, accuracy
                    )

                for ip_x, ip_y in intersections:
                    ip_z = ip_x + 1j * ip_y

                    # t_n: index of new_s_wall.z nearest to ip_z
                    t_n = get_nearest_point_index(
                        new_s_wall.z, ip_z, sw_data.branch_points, accuracy,
                    )

                    # t_p: index of z_seg_p nearest to ip_z
                    t_p = get_nearest_point_index(
                        prev_s_wall.z, ip_z, sw_data.branch_points, accuracy,
                    )

                    # TODO: need to put the joint into the parent
                    # S-walls?

                    logger.debug('Joint at z = {}'.format(ip_z))

                    new_joints.append(
                        Joint(
                            z=ip_z, 
                            s_wall_1=prev_s_wall,
                            s_wall_2=new_s_wall,                         
                            t_1=t_p, 
                            t_2=t_n,
                            sw_data=sw_data,
                        )
                    )

        return new_joints

#    def get_new_joints_using_cgal(self, new_s_wall_index, config,
#                                  sw_data, linux_distribution=None):
#        """
#        Find new wall-wall intersections using CGAL 2d curve intersection.
#        """
#        logger = logging.getLogger(self.logger_name)
#        accuracy = config['accuracy']
#        new_joints = []
#
#        lib_name = 'libcgal_intersection'
#        if linux_distribution == 'Ubuntu':
#            lib_name += '_ubuntu'
#        elif linux_distribution == 'debian':
#            # FIXME: Anaconda Python returns 'debian' instead of 'Ubuntu'.
#            lib_name += '_ubuntu'
#        elif linux_distribution == 'Scientific Linux':
#            lib_name += '_het-math2'
#        else:
#            raise OSError
#
#        # Load CGAL shared library.
#        libcgal_intersection = numpy.ctypeslib.load_library(
#            lib_name, 
#            os.path.dirname(os.path.realpath(__file__)) + '/cgal_intersection/'
#        )
#        cgal_find_intersections_of_curves = (libcgal_intersection.
#                                             find_intersections_of_curves)
#        # Prepare types for CGAL library.
#        array_2d_float = numpy.ctypeslib.ndpointer(
#            dtype=numpy.float64,
#            ndim=2,
#            flags=['C_CONTIGUOUS', 'ALIGNED'],
#        )
#        array_2d_complex = numpy.ctypeslib.ndpointer(
#            dtype=numpy.complex128,
#            ndim=1,
#            flags=['C_CONTIGUOUS', 'ALIGNED'],
#        )
#
#        cgal_find_intersections_of_curves.restype = ctypes.c_int
#        cgal_find_intersections_of_curves.argtypes = [
#            # array_2d_float,
#            array_2d_complex,
#            ctypes.c_long,
#            # array_2d_float, 
#            array_2d_complex,
#            ctypes.c_long,
#            array_2d_float, ctypes.c_int,
#        ]
#
#        new_s_wall = self.s_walls[new_s_wall_index]
#
#        for prev_s_wall in self.s_walls[:new_s_wall_index]:
#
#            # First check if the two S-walls are compatible
#            # for forming a joint.
#
#            # 1. Check if the new S-wall is a descendant
#            # of an existing S-wall. 
#            if prev_s_wall.label in new_s_wall.parents:
#                continue
#            
#            # 2. Split the two S-walls into segments 
#            # according to the trivialization, then
#            # check the compatibility of a pair
#            # of segments.
#            new_z_segs = new_s_wall.get_z_segs()
#            num_new_z_segs = len(new_z_segs)
#            prev_z_segs = prev_s_wall.get_z_segs()
#            num_prev_z_segs = len(prev_z_segs)
#                
#            for n_i, p_i in itertools.product(
#                range(num_new_z_segs), range(num_prev_z_segs)
#            ):
#                n_seg_zs = new_z_segs[n_i]
#                n_seg_root = new_s_wall.local_roots[n_i]
#                p_seg_zs = prev_z_segs[p_i]
#                p_seg_root = prev_s_wall.local_roots[p_i]
#
#                if not is_root(p_seg_root + n_seg_root, sw_data.g_data):
#                    # The two segments are not compatible for
#                    # forming a joint.
#                    continue
#                # Find an intersection of the two segments on the z-plane.
#                buffer_size = 10
#                intersection_search_finished = False
#                while not intersection_search_finished:
#                    intersections = numpy.empty((buffer_size, 2), 
#                                                dtype=numpy.float64)
#                    num_of_intersections = cgal_find_intersections_of_curves(
#                        n_seg_zs,
#                        ctypes.c_long(len(n_seg_zs)), 
#                        p_seg_zs,
#                        ctypes.c_long(len(p_seg_zs)), 
#                        intersections, buffer_size
#                    )
#                    if num_of_intersections == 0:
#                        intersection_search_finished = True
#                    elif num_of_intersections > buffer_size:
#                        logger.info('Number of intersections larger than '
#                                    'the buffer size; increase its size '
#                                    'to {} and find intersections again.'
#                                    .format(num_of_intersections))
#                        buffer_size = num_of_intersections
#                    else:
#                        intersection_search_finished = True
#
#                for ip_x, ip_y in intersections:
#                    ip_z = ip_x + 1j * ip_y
#
#                    # t_n: index of new_s_wall.z nearest to ip_z
#                    t_n = get_nearest_point_index(
#                        new_s_wall.z, ip_z, sw_data.branch_points, accuracy,
#                    )
#
#                    # t_p: index of z_seg_p nearest to ip_z
#                    t_p = get_nearest_point_index(
#                        prev_s_wall.z, ip_z, sw_data.branch_points, accuracy,
#                    )
#
#                    # TODO: need to put the joint into the parent
#                    # S-walls?
#
#                    logger.debug('Joint at z = {}'.format(ip_z))
#
#                    new_joints.append(
#                        Joint(
#                            z=ip_z, 
#                            s_wall_1=prev_s_wall,
#                            s_wall_2=new_s_wall,                         
#                            t_1=t_p, 
#                            t_2=t_n,
#                            sw_data=sw_data,
#                        )
#                    )
#
#        return new_joints
#
#    def get_new_joints_using_interpolation(
#        self, new_s_wall_index, config, sw_data, 
#    ):
#        """
#        Find joints between the newly grown segment of the given S-wall
#        and the other S-walls by interpolating S-walls with functions and
#        finding roots of the pairwise differences of the functions.
#
#        This checks joints that are formed by two
#        S-walls only, not considering the possibility of a joint of three
#        S-walls, which in principle can happen but is unlikely in a numerical
#        setup.
#        """
#        logger = logging.getLogger(self.logger_name)
#        new_joints = []
#
#        new_s_wall = self.s_walls[new_s_wall_index]
#        new_tps = new_s_wall.get_turning_points()
#        new_z_segs = numpy.split(new_s_wall.z, new_tps, axis=0,)
#
#        # NOTE: Here we find only a single joint between two S-walls.
#        # Use CGAL to find multiple z-intersections.
#        for prev_s_wall in self.s_walls[:new_s_wall_index]:
#            prev_tps = prev_s_wall.get_turning_points()
#            prev_z_segs = numpy.split(prev_s_wall.z, prev_tps, axis=0,)
#
#            for i_n in range(len(new_tps) + 1):
#                z_seg_n = new_z_segs[i_n]
#                for i_p in range(len(prev_tps) + 1):
#                    z_seg_p = prev_z_segs[i_p]
#                    # Find an intersection on the z-plane.
#                    try:
#                        ip_x, ip_y = find_intersection_of_segments(
#                            (z_seg_n.real, z_seg_n.imag),
#                            (z_seg_p.real, z_seg_p.imag),
#                            config['accuracy'],
#                        )
#                        ip_z = ip_x + 1j * ip_y
#
#                        # t_n: index of z_seg_n nearest to ip_z
#                        t_n = n_nearest_indices(new_s_wall.z, ip_z, 1)[0]
#
#                        # t_p: index of z_seg_p nearest to ip_z
#                        t_p = n_nearest_indices(prev_s_wall.z, ip_z, 1)[0]
#
#                        # TODO: need to put the joint into the parent
#                        # S-walls?
#
#                        # find mass of parent S-walls: this is approximate,
#                        # since we don't interpolate precisely to the joint
#                        # TODO: improve by doing precise interpolation
#                        logger.debug(
#                            'evaluating possible joint at z = {}'.format(ip_z)
#                        )
#                        a_joint = get_joint(
#                            ip_z, 
#                            prev_s_wall,
#                            new_s_wall, 
#                            t_p,
#                            t_n,
#                            sw_data,
#                        )
#
#                        if(a_joint is None):
#                            continue
#                        else:
#                            new_joints.append(a_joint)
#
#                    except NoIntersection:
#                        pass
#        return new_joints


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


def get_nearest_point_index(s_wall_z, p_z, branch_points, accuracy,
                            logger_name='loom',):
    """
    Get the index of the point on the S-wall that is nearest to 
    the given point on the z-plane.

    When the point found is within the accuracy limit from a branch cut,
    look fot the next nearest point and return its index.
    """
    logger = logging.getLogger(logger_name)

    ts = n_nearest_indices(s_wall_z, p_z, 3)
    zs = [s_wall_z[t] for t in ts]

    t = ts[0]

    for bp in branch_points:
        if abs(zs[0].real - bp.z.real) < accuracy:
            logger.info(
                'The nearest point is too close to a branch cut, '
                'find the next nearest point.'
            )
            if (p_z.real - bp.z.real)*(zs[1].real - bp.z.real) > 0:
                t = ts[1]
                break
            elif (p_z.real - bp.z.real)*(zs[2].real - bp.z.real) > 0:
                t = ts[2]
                break
            else:
                logger.warning(
                    'Unable to find the next nearest point '
                    'that is on the same side from the branch cut '
                    'as the reference point.'
                )
                break

    # Final check.
    for bp in branch_points:
        if abs(s_wall_z[t].real - bp.z.real) < accuracy:
            logger.warning('Unable to find the nearest point on the S-wall '
                           'that is outside the accuracy limit '
                           'from branch cuts of {}.'.format(bp.label))
    return t


def find_intersections_of_curves(a_zs, b_zs, accuracy):

    a_tps = get_turning_points(a_zs)
    a_z_segs = numpy.split(a_zs, a_tps, axis=0,)

    b_tps = get_turning_points(b_zs)
    b_z_segs = numpy.split(b_zs, b_tps, axis=0,)

    intersections = []

    for i_a in range(len(a_tps) + 1):
        z_seg_a = a_z_segs[i_a]
        for i_b in range(len(b_tps) + 1):
            z_seg_b = b_z_segs[i_b]

            # Find an intersection on the z-plane.
            try:
                ip_x, ip_y = find_intersection_of_segments(
                    (z_seg_a.real, z_seg_a.imag),
                    (z_seg_b.real, z_seg_b.imag),
                    accuracy,
                )
                intersections.append((ip_x, ip_y))

            except NoIntersection:
                pass

    return intersections


def get_turning_points(zs):
    """
    Return a list of indices of turning points of a curve on the z-plane,
    i.e. dx/dy = 0 or dy/dx = 0, where x = z[t].real and y = z[t].imag.
    """
    tps = []

    if len(zs) < 3:
        return tps

    x_0 = zs[0].real
    y_0 = zs[0].imag

    for t in range(1, len(zs) - 1):
        x_1 = zs[t].real
        y_1 = zs[t].imag
        x_2 = zs[t + 1].real
        y_2 = zs[t + 1].imag
        if (
            (x_1 - x_0) * (x_2 - x_1) < 0 or (y_1 - y_0) * (y_2 - y_1) < 0
        ): 
            tps.append(t)
        x_0 = x_1
        y_0 = y_1
    return tps


