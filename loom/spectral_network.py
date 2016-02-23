import platform
import os
import numpy
import sympy
import ctypes
import logging

import itertools

from cmath import exp
from scipy import integrate
from s_wall import SWall, Joint, get_s_wall_seeds
from misc import (
    n_nearest_indices, get_turning_points, get_splits_with_overlap,
    get_descendant_roots, sort_roots,
)
from intersection import (
    NoIntersection, find_intersection_of_segments,
)


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
        self.data_attributes = ['phase', 's_walls', 'joints']

    def set_z_rotation(self, z_rotation):
        for s_wall in self.s_walls:
            s_wall.set_z_rotation(z_rotation)
        for joint in self.joints:
            joint.set_z_rotation(z_rotation)

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
            
            get_intersections = (libcgal_intersection.
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

            get_intersections.restype = ctypes.c_int
            get_intersections.argtypes = [
                # array_2d_float,
                array_2d_complex,
                ctypes.c_long,
                # array_2d_float, 
                array_2d_complex,
                ctypes.c_long,
                array_2d_float, ctypes.c_int,
            ]

            logger.info('Use CGAL to find intersections.')
            use_cgal = True

        except OSError:
            logger.warning('CGAL not available; use interpolation '
                           'to find intersections.')
            get_intersections = find_intersections_of_curves
            use_cgal = False

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
                        # TODO: change this to take bp itself.
                        parents=[bp.label],
                        parent_roots=[root for root in bp.positive_roots],
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
            logger.info('Start iteration #{}...'.format(iteration))
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
                )
                # Cut the grown S-walls at the intersetions with branch cuts
                # and decorate each segment with its root data.
                logger.info('Determining the root type of S-wall #{}...'
                            .format(i))
                s_i.determine_root_types(
                    sw_data, cutoff_radius=config['size_of_small_step'],
                )

                logger.info('Finding new joints from S-wall #{}...'.format(i))
                new_joints += self.get_new_joints(i, config, sw_data,
                                                  get_intersections,
                                                  use_cgal,)

            n_finished_s_walls = len(self.s_walls)
            if(len(new_joints) == 0):
                logger.info(
                    'No additional joint found: '
                    'Stop growing this spectral network '
                    'at iteration #{}.'.format(iteration)
                )
                break
            elif iteration == num_of_iterations:
                # Last iteration finished. Do not form new joints,
                # and do not seed additional S-walls, either.
                break
            else:
                logger.info(
                    'Growing S-walls in iteration #{} finished.'
                    .format(iteration)
                )

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
                            # XXX: parent_roots != (roots of parents)
                            parent_roots=joint.roots,
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
        branch_loci = sw_data.branch_points + sw_data.irregular_singularities

        self.phase = json_data['phase']

        for s_wall_data in json_data['s_walls']:
            an_s_wall = SWall(logger_name=self.logger_name,)
            an_s_wall.set_from_json_data(s_wall_data, branch_loci)
            self.s_walls.append(an_s_wall)

        for joint_data in json_data['joints']:
            a_joint = Joint()
            a_joint.set_from_json_data(joint_data)
            self.joints.append(a_joint)

    def get_new_joints(self, new_s_wall_index, config, sw_data,
                       get_intersections, use_cgal,):
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
            logger.info(
                'There is no joint for the given root system {}.'
                .format(config['root_system'])
            )
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
            n_z_splits = new_s_wall.get_splits(endpoints=True)
            num_n_z_segs = len(new_s_wall.local_roots)
            p_z_splits = prev_s_wall.get_splits(endpoints=True)
            num_p_z_segs = len(prev_s_wall.local_roots)
                
            for n_z_seg_i, p_z_seg_i in itertools.product(
                range(num_n_z_segs), range(num_p_z_segs)
            ):
                n_z_i = n_z_splits[n_z_seg_i]
                n_z_f = n_z_splits[n_z_seg_i + 1]
                p_z_i = p_z_splits[p_z_seg_i]
                p_z_f = p_z_splits[p_z_seg_i + 1]

#                n_seg_root = new_s_wall.local_roots[n_z_seg_i]
#                p_seg_root = prev_s_wall.local_roots[p_z_seg_i]
                descendant_roots = get_descendant_roots(
                    (prev_s_wall.multiple_local_roots[p_z_seg_i] +
                     new_s_wall.multiple_local_roots[n_z_seg_i]),
                    sw_data.g_data,
                )

#                if not is_root(p_seg_root + n_seg_root, sw_data.g_data):
                if len(descendant_roots) == 0:
                    # The two segments are not compatible for
                    # forming a joint.
                    continue

                # Find an intersection of the two segments on the z-plane.
                if use_cgal is True:
                    buffer_size = 10
                    intersection_search_finished = False

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
                        prev_s_wall.parents == new_s_wall.parents 
                        and abs(ip_z - prev_s_wall.z[0]) < accuracy
                        and abs(ip_z - new_s_wall.z[0]) < accuracy
                    ):
                        continue
                    elif (
                        n_nearest_indices(prev_s_wall.z, ip_z, 1)[0] == 0
                        or n_nearest_indices(new_s_wall.z, ip_z, 1)[0] == 0
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

                    logger.debug('Intersection at z = {}'.format(ip_z))

                    # TODO: check if the following descendant-roots
                    # finding is necessary, note that we calculate
                    # descendant roots above.
                    descendant_roots = get_descendant_roots(
                        (prev_s_wall.get_roots_at_t(t_p) +
                         new_s_wall.get_roots_at_t(t_n)),
                        sw_data.g_data,
                    )

                    joint_data = get_joint_data(
                        descendant_roots, ip_z, sw_data
                    )
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

                    joint_M = prev_s_wall.M[t_p] + new_s_wall.M[t_n]
                    joint_parents = [prev_s_wall.label, new_s_wall.label]
                    for roots, ode_xs in joint_data_groups:
                        new_joints.append(
                            Joint(
                                z=ip_z,
                                M=joint_M,
                                ode_xs=ode_xs,
                                parents=joint_parents,
                                roots=sort_roots(roots, sw_data.g_data),
                            )
                        )
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
    v = sympy.lambdify((z, x), sw.diff.num_v)

    def ode_f(t, z_x1_x2_M):
        z_i = z_x1_x2_M[0]
        x1_i = z_x1_x2_M[1]
        x2_i = z_x1_x2_M[2]
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

    t_0 = n_nearest_indices(s_wall_z, p_z, 1)[0]

    t = t_0

    t_max = len(s_wall_z) - 1

    t_before = t_0

    t_after = t_0

    for bp in branch_points:
        if bp.z.real - p_z.real == 0.0 and bp.z.imag <= p_z.imag:
            raise Exception(
                'Intersection point lies exactly above {}. '
                'Try perturbing the phase.'
                .format(bp.label)
            )

        elif (
            abs(s_wall_z[t].real - bp.z.real) < accuracy 
            and s_wall_z[t].imag - bp.z.imag > 0
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
                        (p_z.real - bp.z.real) * (z_m.real - bp.z.real) > 0
                        and abs(s_wall_z[t_m].real - bp.z.real) > accuracy
                        and t_m <= t_before
                    ):
                        t_before = t_m

                        break

                t_p += 1
                if t_p <= t_max:
                    z_p = s_wall_z[t_p]
                    if (
                        (p_z.real - bp.z.real) * (z_p.real - bp.z.real) > 0
                        and abs(s_wall_z[t_p].real - bp.z.real) > accuracy
                        and t_p >= t_after
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
        raise Exception(
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


def get_joint_data(descendant_roots, z, sw_data):
    """
    From a set of descendant roots find a set of joint data,
        [(root, ode_xs), ...].
    """
    joint_data = []
    for root in descendant_roots:
        # FIXME: The following, including self.ode_xs, can be removed
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

    return joint_data
