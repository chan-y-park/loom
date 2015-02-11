import numpy
import sympy
import logging
import signal
import multiprocessing
import time
import os
import time
import zipfile, zlib
import glob
import json
import pdb

from cmath import exp
from itertools import combinations

from scipy import integrate

from curve import RamificationPoint, SWCurve, SWDiff
from s_wall import SWall, Joint, get_s_wall_seeds, get_joint
from misc import (n_nearest, n_nearest_indices, find_xs_at_z_0)
from intersection import (HitTable, NoIntersection,
                          find_intersection_of_segments)
from plotting import SpectralNetworkPlot, plot_segments
#from data_io import save_spectral_network_data, load_spectral_network_data

x, z = sympy.symbols('x z')


class SpectralNetwork:
    def __init__(
        self,
        #sw_curve,
        #sw_diff,
        phase=None,
        ramification_points=[],
        config=None,
    ):
        self.phase = phase
        self.ramification_points = ramification_points
        self.hit_table = HitTable(config['size_of_bin'])

        self.ode = None
        self.s_walls = []
        self.joints = []


    def seed(self, sw_curve, sw_diff, config):
        # Initialize S-walls at ramification points
        for rp in self.ramification_points:
            s_wall_seeds = get_s_wall_seeds(
                sw_curve, sw_diff, self.phase, rp, config
            )
            for z_0, x1_0, x2_0 in s_wall_seeds:
                label = 'S-wall #{}'.format(len(self.s_walls))
                self.s_walls.append(SWall(z_0, x1_0, x2_0, [rp.label], label))
    
    def get_data(self):
        data = {
            'phase': self.phase,
            'ramification_points' : self.ramification_points,
            's_walls': self.s_walls,
            'joints': self.joints,
            'hit_table': self.hit_table,
        }
        return data

    def set_ode(self, sw_curve, sw_diff, config):
        ode_absolute_tolerance = config['accuracy']

        f = sw_curve.num_eq
        df_dz = f.diff(z)
        df_dx = f.diff(x)
        #F = -(\partial f/\partial z)/(\partial f/\partial x)
        F = sympy.lambdify((z, x), -df_dz/df_dx)
        v = sympy.lambdify((z, x), sw_diff.num_v)

        def ode_f(t, zx1x2):
            z_i = zx1x2[0]
            x1_i = zx1x2[1]
            x2_i = zx1x2[2]
            dz_i_dt = exp(self.phase*1j)/(v(z_i, x1_i) - v(z_i, x2_i))
            dx1_i_dt = F(z_i, x1_i) * dz_i_dt
            dx2_i_dt = F(z_i, x2_i) * dz_i_dt
            return [dz_i_dt, dx1_i_dt, dx2_i_dt]

        self.ode = integrate.ode(ode_f)
        self.ode.set_integrator(
            'zvode',
            #method='adams',
            atol=ode_absolute_tolerance,
        )

    def get_new_joints(self, new_s_wall_index, sw_curve, sw_diff, config):
        """
        Find joints between the newly grown segment of the given S-wall
        and the other S-walls. This checks joints that are formed by two
        S-walls only, not considering the possibility of a joint of three
        S-walls, which in principle can happen but is unlikely in a numerical
        setup.
        """
        new_joints = []
        if (config['root_system'] in ['A2',]):
            # There is no joint for the given root system.
            return new_joints

        i_c = new_s_wall_index
        new_curve = self.s_walls[i_c].get_zxzys()
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
                logging.debug('i_c, i_d = {}, {}'.format(i_c, i_d))

                for t_c_i, t_c_f in self.hit_table[bin_key][i_c]:
                    for t_d_i, t_d_f in self.hit_table[bin_key][i_d]:
                        segment_c = self.s_walls[i_c].get_zxzys(t_c_i, t_c_f+1)
                        segment_d = self.s_walls[i_d].get_zxzys(t_d_i, t_d_f+1)

                        try:
                            # find an intersection on the z-plane.
                            ip_x, ip_y = find_intersection_of_segments(
                                segment_c,
                                segment_d,
                                self.hit_table.get_bin_location(bin_key),
                                self.hit_table.get_bin_size(),
                                config['accuracy']
                            )
                            ip_z = ip_x + 1j*ip_y

                            # segment_c
                            zs_c = self.s_walls[i_c].get_zs(t_c_i, t_c_f+1)
                            ip_t_c = (n_nearest_indices(zs_c, ip_z, 1)[0] +
                                      t_c_i)
                            z_c, x1_c, x2_c = self.s_walls[i_c].data[ip_t_c]

                            # segment_d
                            zs_d = self.s_walls[i_d].get_zs(t_d_i, t_d_f+1)
                            ip_t_d = (n_nearest_indices(zs_d, ip_z, 1)[0] +
                                      t_d_i)
                            z_d, x1_d, x2_d = self.s_walls[i_d].data[ip_t_d]

                            # TODO: need to put the joint into the parent
                            # S-walls?

                            # find the values of x at z = ip_z.
                            ip_xs = find_xs_at_z_0(sw_curve.num_eq, ip_z)
                            ip_x1_c = n_nearest(ip_xs, x1_c, 1)[0]
                            ip_x2_c = n_nearest(ip_xs, x2_c, 1)[0]
                            ip_x1_d = n_nearest(ip_xs, x1_d, 1)[0]
                            ip_x2_d = n_nearest(ip_xs, x2_d, 1)[0]
                            
                            a_joint_label = ('joint ' +
                                             '#{}'.format(len(self.joints)))
                            a_joint = get_joint(
                                ip_z, ip_x1_c, ip_x2_c, ip_x1_d, ip_x2_d,
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

    def grow(self, sw_curve, sw_diff, config):
        rpzs = []
        ppzs = []

        for rp in self.ramification_points:
            if rp.is_puncture is False:
                rpzs.append(rp.z)
            elif rp.is_puncture is True:
                ppzs.append(rp.z)

        n_finished_s_walls = 0 
        iteration = 0
        while(iteration < config['num_of_iterations']):
            """
            Iterate until there is no new joint.
            Each S-wall is grown only once.
            """
            new_joints = []     # number of new joints found in each iteration
            # first grow seeded S-walls.
            for i in range(n_finished_s_walls, len(self.s_walls)):
                self.s_walls[i].grow(
                    self.ode, rpzs, ppzs, 
                    z_range_limits=config['z_range_limits'],
                    num_of_steps=config['num_of_steps'],
                    size_of_small_step=config['size_of_small_step'],
                    size_of_large_step=config['size_of_large_step'],
                    size_of_neighborhood=config['size_of_neighborhood'],
                )
                new_joints += self.get_new_joints(i, sw_curve, sw_diff,
                                                  config)

            n_finished_s_walls = len(self.s_walls)
            if(len(new_joints) == 0):
                break

            # then seed an S-wall for each new joint.
            for joint in new_joints:
                joint_is_new = True
                # check if this is not in the previous joint list
                for old_joint in self.joints:
                    if joint.is_equal_to(old_joint, config['accuracy']):
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


def init_process():
    """
    Initializer of each process that generates a spectral network.
    Take care of a keyboard interrupt.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def parallel_get_spectral_network(
    sw_curve, 
    sw_diff,
    phase,
    ramification_points,
    config,
    data_save_dir,
    shared_n_started_spectral_networks,
    shared_n_finished_spectral_networks
):
    shared_n_started_spectral_networks.value += 1
    job_id = shared_n_started_spectral_networks.value
    logging.info('Start generating spectral network #{}: theta = {}.'.format(
        job_id, phase
    ))

    spectral_network = SpectralNetwork(phase, ramification_points, config) 

    spectral_network.seed(sw_curve, sw_diff, config)

    spectral_network.set_ode(sw_curve, sw_diff, config)

    spectral_network.grow(sw_curve, sw_diff, config)

    shared_n_finished_spectral_networks.value += 1
    logging.info('Finished generating spectral network #{}.'.format(
        shared_n_finished_spectral_networks.value
    ))

    spectral_network_data = spectral_network.get_data()

    # Save spectral network data to a file
    data_file_name = os.path.join(data_save_dir,
                                  'data_{}.json'.format(job_id))
    logging.info('Saving data to {}.'.format(data_file_name))
    with open(data_file_name, 'wb') as fp:
        save_spectral_network_data(
            spectral_network_data, 
            fp, 
        )

    #return spectral_network
    return spectral_network_data


def save_spectral_network_data(data, file_object, **kwargs):
    json_data = {}
    json_data['phase'] = data['phase']
    json_data['ramification_points'] = [rp.get_json_data()
                                        for rp in data['ramification_points']]
    json_data['s_walls'] = [s_wall.get_json_data()
                            for s_wall in data['s_walls']]
    json_data['joints'] = [joint.get_json_data()
                           for joint in data['joints']]
    json_data['hit_table'] = data['hit_table'].get_json_data()
    json.dump(json_data, file_object, **kwargs)
    

def generate_spectral_network(opts, config):
    sw_curve = SWCurve(config)
    sw_diff = SWDiff(config)
    spectral_network_data_list = []
    file_list = []
    ramification_points = sw_curve.get_ramification_points()
    start_time = time.time()
    logging.info('start cpu time: %s', start_time)

    logging.debug('\nList of ramification points')
    for rp in ramification_points:
        logging.debug('%s', rp)
    logging.debug('\n')

    # Prepare to save spectral network data to files.
    timestamp = str(int(time.time()))
    data_save_dir = os.path.join(
        config['root_dir'], 
        config['data_dir'], 
        timestamp
    )
    os.makedirs(data_save_dir)
    # Save configuration to a file.
    config_file_name = os.path.join(data_save_dir, 'config.ini')
    with open(config_file_name, 'wb') as fp:
        config.parser.write(fp)
        file_list.append(config_file_name)

    if(opts['phase'] is not None):
        # Generate a single spectral network.
        spectral_network = SpectralNetwork(
            phase=opts['phase'], 
            ramification_points=ramification_points,
            config=config
        ) 

        spectral_network.seed(sw_curve, sw_diff, config)

        spectral_network.set_ode(sw_curve, sw_diff, config)

        spectral_network.grow(sw_curve, sw_diff, config)

        spectral_network_data_list.append(spectral_network.get_data())

        # Save spectral network data to a file
        data_file_name = os.path.join(data_save_dir, 'data_0.json')
        logging.info('Saving data to {}.'.format(data_file_name))
        with open(data_file_name, 'wb') as fp:
            save_spectral_network_data(
                spectral_network_data_list[0], 
                fp, 
            )
            file_list.append(data_file_name)

    elif(config['phase_range'] is not None):
        # Generate multiple spectral networks.
        theta_i, theta_f, theta_n = config['phase_range']
        phases = [theta_i + i * (theta_f - theta_i) / theta_n
                  for  i in range(theta_n)]

        manager = multiprocessing.Manager()
        shared_n_started_spectral_networks = manager.Value('i', 0)
        shared_n_finished_spectral_networks = manager.Value('i', 0)
        if(config['n_processes'] == 0):
            n_processes = multiprocessing.cpu_count()
        elif(config['n_processes'] < 0):
            n_processes = multiprocessing.cpu_count()
            if(n_processes > config['n_processes']):
                n_processes -= config['n_processes']
            else:
                logging.warning('The number of CPUs is smaller than '
                                '{}.'.format(-config['n_processes']))
                logging.warning('Set n_processes to 1.')
                n_processes = 1
        else:
            n_processes = config['n_processes']
        logging.info('Number of processes in the pool: {}'.format(n_processes))
        pool =  multiprocessing.Pool(n_processes, init_process)
        try:
            results = [
                pool.apply_async(
                    parallel_get_spectral_network, 
                    args=(
                        sw_curve,
                        sw_diff,
                        phase,
                        ramification_points,
                        config,
                        data_save_dir,
                        shared_n_started_spectral_networks,
                        shared_n_finished_spectral_networks,
                    )
                ) for phase in phases
            ]
            pool.close()

            for result in results:
                spectral_network_data_list.append(result.get())
                logging.info('job progress: {}/{} finished.'.format(
                    shared_n_finished_spectral_networks.value,
                    theta_n
                ))

            file_list += glob.glob(os.path.join(data_save_dir, 'data_*.json'))

        except KeyboardInterrupt:
            logging.warning('Caught ^C; terminates processes...')
            pool.terminate()
            pool.join()
        # End of generating multiple spectral networks

    end_time = time.time()
    logging.info('end cpu time: %.8f', end_time)
    logging.info('elapsed cpu time: %.8f', end_time - start_time)

    # Make a compressed data file.
    zipped_file_name = data_save_dir + '.zip'
    logging.info('Save compressed data to {}.'.format(zipped_file_name))
    with zipfile.ZipFile(zipped_file_name, 'w', zipfile.ZIP_DEFLATED) as fp:
        for a_file in file_list:
            fp.write(a_file, os.path.relpath(a_file, data_save_dir))

    # Plot spectral networks.
    if(opts['show-plot'] is True):
        spectral_network_plot = SpectralNetworkPlot(
            config,
            #plot_data_points=True,
            #plot_joints=True,
            #plot_bins=True,
            #plot_segments=True,
        )

        if (len(spectral_network_data_list) > 0):
            for data in spectral_network_data_list:
                spectral_network_plot.set_data(data)

        spectral_network_plot.show()

    return spectral_network_data_list


def load_spectral_network_data(file_object, sw_curve, sw_diff, config,
                               **kwargs):
    json_data = json.load(file_object, **kwargs)
    spectral_network = SpectralNetwork(
        phase=json_data['phase'],
        config=config,
    )

    for rp_data in json_data['ramification_points']:
        rp = RamificationPoint()
        rp.set_json_data(rp_data)
        spectral_network.ramification_points.append(rp)

    for joint_data in json_data['joints']:
        a_joint = Joint()
        a_joint.set_json_data(joint_data)
        spectral_network.joints.append(a_joint)

    for s_wall_data in json_data['s_walls']:
        an_s_wall = SWall()
        an_s_wall.set_json_data(s_wall_data)
        spectral_network.s_walls.append(an_s_wall)

    spectral_network.hit_table.load_json_data(json_data['hit_table'])

    # NOTE: Loaded data indicates parents as labels,
    # but they should be changed to the corresponding instances.

    return spectral_network

def load_spectral_network(data_dir, config):
    sw_curve = SWCurve(config)
    sw_diff = SWDiff(config)
    spectral_network_data_list = []

    data_file_list = glob.glob(os.path.join(data_dir, 'data_*.json'))
    for data_file in data_file_list:
        with open(data_file, 'r') as fp:
            spectral_network = load_spectral_network_data(
                fp, sw_curve, sw_diff, config
            )
            spectral_network_data_list.append(spectral_network.get_data())

    # Make plots from the loaded data
    spectral_network_plot = SpectralNetworkPlot(
        config,
        #sw_curve,
        #sw_diff,
    )
    for data in spectral_network_data_list:
        spectral_network_plot.set_data(data)
    spectral_network_plot.show()

    return spectral_network_data_list

