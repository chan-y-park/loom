import numpy
import sympy
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
                          find_intersection_of_segments)
from plotting import SpectralNetworkPlot, plot_segments

x, z = sympy.symbols('x z')


class SpectralNetwork:
    def __init__(
        self,
        phase=None,
        ramification_points=[],
        config=None,
    ):
        self.phase = phase
        self.ramification_points = ramification_points
        self.hit_table = HitTable(config['size_of_bin'])

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
                        n_steps=config['num_of_steps']
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
                self.s_walls[i].grow(ode, rpzs, ppzs, config,)
                new_joints += get_new_joints(self, i, sw, config)

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
                    SWall(joint.z, joint.x1, joint.x2, joint.parents, label)
                )
            iteration += 1

            logging.info('Iteration #{} finished.'.format(iteration))


    def get_data(self):
        """
        Returns a dictionary of spectral network data that is used 
        when comminicating with parallel children processes
        that generates spectral networks.
        """
        data = {
            'phase': self.phase,
            'ramification_points' : self.ramification_points,
            's_walls': self.s_walls,
            'joints': self.joints,
            'hit_table': self.hit_table,
        }
        return data


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
            spectral_network.s_walls.append(an_s_wall)

        for joint_data in json_data['joints']:
            a_joint = Joint()
            a_joint.set_json_data(joint_data)
            self.append(a_joint)

        spectral_network.hit_table.load_json_data(json_data['hit_table'])


def get_new_joints(spectral_network, new_s_wall_index, sw, config):
    """
    Find joints between the newly grown segment of the given S-wall
    and the other S-walls. This checks joints that are formed by two
    S-walls only, not considering the possibility of a joint of three
    S-walls, which in principle can happen but is unlikely in a numerical
    setup.
    """
    s_walls = spectral_network.s_walls

    new_joints = []
    if (config['root_system'] in ['A2',]):
        # There is no joint for the given root system.
        return new_joints

    i_c = new_s_wall_index
    new_s_wall = s_walls[i_c]
    new_curve = numpy.array([new_s_wall.z.real, new_s_wall.z.imag]).T
    new_bin_keys = self.hit_table.fill(i_c, new_curve)

    # first get the intersections on the z-plane
    for bin_key in new_bin_keys:
        if len(self.hit_table[bin_key]) == 1:
            # only one S-wall in the bin, skip the rest.
            continue
        for i_a, i_b in combinations(self.hit_table[bin_key], 2):
            # don't check self-intersections.
            if i_a == i_c:
                i_d = i_b   # i_b != i_c
            elif i_b == i_c:
                i_d = i_a   # i_a != i_c
            else:
                # both segments are not in the newly added curve.
                # don't check the intersection.
                continue

            s_wall_c = s_walls[i_c]
            s_wall_d = s_walls[i_d]

            for t_c_i, t_c_f in self.hit_table[bin_key][i_c]:
                for t_d_i, t_d_f in self.hit_table[bin_key][i_d]:
                    segment_c = s_wall_c.z[t_c_i:t_c_f+1]
                    segment_d = s_wall_c.z[t_d_i:t_d_f+1]
                    try:
                        # find an intersection on the z-plane.
                        ip_x, ip_y = find_intersection_of_segments(
                            (segment_c.real, segment_c.imag),
                            (segment_d.real, segment_d.imag),
                            self.hit_table.get_bin_location(bin_key),
                            self.hit_table.get_bin_size(),
                            config['accuracy']
                        )
                        ip_z = ip_x + 1j*ip_y

                        # segment_c
                        seg_c_z = s_wall_c.z[t_c_i:t_c_f+1]
                        # index of z of segment_c nearest to ip_z
                        ip_t_c = (n_nearest_indices(seg_c_z, ip_z, 1)[0] +
                                  t_c_i)
                        #z_c = s_wall_c.z[ip_t_c]
                        c_x = s_wall_c.x[ip_t_c]

                        # segment_d
                        seg_d_z = s_wall_d.z[t_d_i:t_d_f+1]
                        # index of z of segment_d nearest to ip_z
                        ip_t_d = (n_nearest_indices(seg_d_z, ip_z, 1)[0] +
                                  t_d_i)
                        #z_d = s_wall_d.z[ip_t_d]
                        d_x = s_wall_d.x[ip_t_d]

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
                            s_wall_c.label, 
                            s_wall_d.label,
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
