import numpy
import logging
import signal
import multiprocessing
import json
import pdb

from cmath import exp
from itertools import combinations

from geometry import RamificationPoint, SWData
from s_wall import SWall, Joint, get_s_wall_seeds, get_joint
from misc import (n_nearest, n_nearest_indices, find_xs_at_z_0, get_ode)
from intersection import (HitTable, NoIntersection,
                          find_intersection_of_segments, 
                          find_curve_range_intersection)
#
from plotting import plot_segments
#


class SpectralNetwork:
    def __init__(
        self,
        phase=None,
        ramification_points=[],
        config=None,
    ):
        self.phase = phase
        self.ramification_points = ramification_points
        self.hit_table = None 
        #self.hit_table = HitTable(config['size_of_bin'])

        self.s_walls = []
        self.joints = []


    def grow(self, sw, config):
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
        for rp in self.ramification_points:
            s_wall_seeds = get_s_wall_seeds(sw, self.phase, rp, config)
            for z_0, x_0 in s_wall_seeds:
                label = 'S-wall #{}'.format(len(self.s_walls))
                self.s_walls.append(
                    SWall(
                        z_0=z_0,
                        x_0=x_0,
                        parents=[rp.label],
                        label=label,
                        n_steps=n_steps,
                    )
                )

        logging.info('Setup the ODE integrator...')
        ode = get_ode(sw, self.phase, accuracy)

        logging.info('Start growing a new spectral network...')
        ppzs = sw.punctures

        rpzs = []
        for rp in self.ramification_points:
            rpzs.append(rp.z)

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
                self.s_walls[i].grow(ode, rpzs, ppzs, config,)
                new_joints += self.get_new_joints(i, sw, config)
                #new_joints += self.get_new_joints_with_hit_table(i, sw, config)

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
                self.joints.append(joint)
                label = 'S-wall #{}'.format(len(self.s_walls))
                self.s_walls.append(
                    SWall(
                        z_0=joint.z,
                        x_0=joint.x, 
                        parents=joint.parents,
                        label=label,
                        n_steps=n_steps,
                    )
                )
            iteration += 1

            logging.info('Iteration #{} finished.'.format(iteration))


#    def get_data(self):
#        """
#        Returns a dictionary of spectral network data that is used 
#        when comminicating with parallel children processes
#        that generates spectral networks.
#        """
#        data = {
#            'phase': self.phase,
#            'ramification_points' : self.ramification_points,
#            's_walls': self.s_walls,
#            'joints': self.joints,
#            'hit_table': self.hit_table,
#        }
#        return data


    def save_json_data(self, file_object, **kwargs):
        """
        Save the spectral network data in a JSON-compatible file.
        """
        json_data = {}

        json_data['phase'] = self.phase
        json_data['ramification_points'] = [
            rp.get_json_data() for rp in self.ramification_points
        ]
        json_data['s_walls'] = [s_wall.get_json_data()
                                for s_wall in self.s_walls]
        json_data['joints'] = [joint.get_json_data()
                               for joint in self.joints]
        if self.hit_table is None:
            json_data['hit_table'] = None 
        else:
            json_data['hit_table'] = self.hit_table.get_json_data()
        json.dump(json_data, file_object, **kwargs)
    

    def set_from_json_data(self, file_object, **kwargs):
        """
        Load the spectral network data from a JSON-compatible file.
        """
        json_data = json.load(file_object, **kwargs)

        self.phase = json_data['phase']

        for rp_data in json_data['ramification_points']:
            rp = RamificationPoint()
            rp.set_json_data(rp_data)
            self.ramification_points.append(rp)

        for s_wall_data in json_data['s_walls']:
            an_s_wall = SWall()
            an_s_wall.set_json_data(s_wall_data)
            self.s_walls.append(an_s_wall)

        for joint_data in json_data['joints']:
            a_joint = Joint()
            a_joint.set_json_data(joint_data)
            self.joints.append(a_joint)

        if json_data['hit_table'] is not None:
            self.hit_table = HitTable()
            self.hit_table.load_json_data(json_data['hit_table'])


    def get_new_joints(self, new_s_wall_index, sw, config):
        """
        Find joints between the newly grown segment of the given S-wall
        and the other S-walls. This checks joints that are formed by two
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

                    #plot_segments(
                    #    [(z_seg_n.real, z_seg_n.imag), 
                    #     (z_seg_p.real, z_seg_p.imag)], 
                    #    [(rp.z.real, rp.z.imag) 
                    #     for rp in self.ramification_points], 
                    #)

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
                        ip_xs = find_xs_at_z_0(sw.curve.num_eq, ip_z)
                        ip_x_n_0 = n_nearest(ip_xs, x_n[0], 1)[0]
                        ip_x_n_1 = n_nearest(ip_xs, x_n[1], 1)[0]
                        ip_x_p_0 = n_nearest(ip_xs, x_p[0], 1)[0]
                        ip_x_p_1 = n_nearest(ip_xs, x_p[1], 1)[0]

                        a_joint = get_joint(
                            ip_z, ip_x_n_0, ip_x_n_1, ip_x_p_0, ip_x_p_1,
                            new_s_wall.label, 
                            prev_s_wall.label,
                            config['accuracy'],
                            root_system=config['root_system'],
                            label='joint #{}'.format(len(self.joints)),
                        )

                        if(a_joint is None):
                            continue
                        else:
                            new_joints.append(a_joint)

                    except NoIntersection:
                        pass
        return new_joints


    def get_new_joints_with_hit_table(self, new_s_wall_index, sw, config):
        """
        Find joints between the newly grown segment of the given S-wall
        and the other S-walls. This checks joints that are formed by two
        S-walls only, not considering the possibility of a joint of three
        S-walls, which in principle can happen but is unlikely in a numerical
        setup.
        """
        if self.hit_table is None:
            self.hit_table = HitTable(config['size_of_bin'])
        new_joints = []
        if (config['root_system'] in ['A1',]):
            logging.info('There is no joint for the given root system {}.'
                         .format(config['root_system']))
            return new_joints

        i_c = new_s_wall_index
        new_s_wall = self.s_walls[i_c]
        new_curve = numpy.array([new_s_wall.z.real, new_s_wall.z.imag]).T
        new_bin_keys = self.hit_table.fill(i_c, new_curve)

        # first get the intersections on the z-plane
        for bin_key in new_bin_keys:
            if len(self.hit_table[bin_key]) == 1:
                # only one S-wall in the bin, skip the rest.
                continue
            for i_a, i_b in combinations(self.hit_table[bin_key], 2):
                # Don't check self-intersections.
                if i_a == i_c:
                    i_d = i_b   # i_b != i_c
                elif i_b == i_c:
                    i_d = i_a   # i_a != i_c
                else:
                    # Both segments are not in the newly added curve.
                    # Don't check the intersection.
                    continue

                for t_c_i, t_c_f in self.hit_table[bin_key][i_c]:
                    for t_d_i, t_d_f in self.hit_table[bin_key][i_d]:
                        # Check if the two segments have a common x-range.  
                        x_c = [x_i for x_i 
                               in self.s_walls[i_c].x[t_c_i:t_c_f+1].T]
                        x_d = [x_i for x_i
                               in self.s_walls[i_d].x[t_d_i:t_d_f+1].T]
                        have_common_x_range = False
                        for x_a, x_b in (
                            (x_c_i, x_d_i) for x_c_i in x_c for x_d_i in x_d
                        ):
                            x_r_range, x_i_range = (
                                find_curve_range_intersection(
                                    (x_a.real, x_a.imag),
                                    (x_b.real, x_b.imag),
                                )
                            )
                            if (x_r_range.is_EmptySet and 
                                x_r_range.is_EmptySet and 
                                x_i_range.is_FiniteSet and
                                x_i_range.is_FiniteSet):
                                continue
                            else:
                                have_common_x_range = True
                                break
                        if have_common_x_range is False:
                            # No common x range, therefore no joint.
                            continue

                        # Find an intersection on the z-plane.
                        seg_c_z = self.s_walls[i_c].z[t_c_i:t_c_f+1]
                        seg_d_z = self.s_walls[i_d].z[t_d_i:t_d_f+1]
                        try:
                            ip_x, ip_y = find_intersection_of_segments(
                                (seg_c_z.real, seg_c_z.imag),
                                (seg_d_z.real, seg_d_z.imag),
                                config['accuracy'],
                                self.hit_table.get_bin_location(bin_key),
                                self.hit_table.get_bin_size(),
                            )
                            ip_z = ip_x + 1j*ip_y

                            # segment_c
                            #seg_c_z = self.s_walls[i_c].z[t_c_i:t_c_f+1]
                            # index of z of segment_c nearest to ip_z
                            ip_t_c = (
                                n_nearest_indices(seg_c_z, ip_z, 1)[0] +
                                t_c_i
                            )
                            #z_c = s_wall_c.z[ip_t_c]
                            c_x = self.s_walls[i_c].x[ip_t_c]

                            # segment_d
                            #seg_d_z = s_wall_d.z[t_d_i:t_d_f+1]
                            # index of z of segment_d nearest to ip_z
                            ip_t_d = (
                                n_nearest_indices(seg_d_z, ip_z, 1)[0] +
                                t_d_i
                            )
                            #z_d = s_wall_d.z[ip_t_d]
                            d_x = self.s_walls[i_d].x[ip_t_d]

                            # TODO: need to put the joint into the parent
                            # S-walls?

                            # find the values of x at z = ip_z.
                            ip_xs = find_xs_at_z_0(sw.curve.num_eq, ip_z)
                            c_ip_x0 = n_nearest(ip_xs, c_x[0], 1)[0]
                            c_ip_x1 = n_nearest(ip_xs, c_x[1], 1)[0]
                            d_ip_x0 = n_nearest(ip_xs, d_x[0], 1)[0]
                            d_ip_x1 = n_nearest(ip_xs, d_x[1], 1)[0]
                            
                            a_joint_label = ('joint ' +
                                             '#{}'.format(len(self.joints)))
                            a_joint = get_joint(
                                ip_z, c_ip_x0, c_ip_x1, d_ip_x0, d_ip_x1,
                                self.s_walls[i_c].label, 
                                self.s_walls[i_d].label,
                                config['accuracy'],
                                root_system=config['root_system'],
                                label=a_joint_label,
                            )

                            if(a_joint is None):
                                continue
                            else:
                                new_joints.append(a_joint)
                            
                        except NoIntersection:
                            pass

        return new_joints
