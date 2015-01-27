import numpy
import sympy
import logging
import pdb
import signal
import multiprocessing

from cmath import exp, pi, phase
from itertools import combinations

from scipy import integrate

from curve import RamificationPoint, SWCurve, SWDiff
from s_wall import SWall, Joint
from misc import (cpow, gather, remove_duplicate, n_nearest, n_nearest_indices,
                  unravel, find_xs_at_z_0, LocalDiffError, GetSWallSeedsError)
from intersection import (HitTable, NoIntersection,
                          find_intersection_of_segments)
from plotting import SpectralNetworkPlot

x, z = sympy.symbols('x z')


class SpectralNetwork:
    def __init__(self, sw_curve, sw_diff, theta, config_data):
        self.sw_curve = sw_curve
        self.sw_diff = sw_diff
        self.theta = theta
        self.config_data = config_data
        self.hit_table = HitTable(config_data.size_of_bin)

        self.s_walls = []
        self.joints = []

        self.set_ode()
        # Initialize S-walls at ramification points
        for rp in self.sw_curve.ramification_points:
            s_wall_seeds = get_s_wall_seeds(self.sw_curve, self.sw_diff,
                                            self.theta, rp, self.config_data)
            for z_0, x1_0, x2_0 in s_wall_seeds:
                label = 'S-wall #{}'.format(len(self.s_walls))
                self.s_walls.append(SWall(z_0, x1_0, x2_0, rp, label))

    def set_ode(self):
        ode_absolute_tolerance = self.config_data.accuracy

        f = self.sw_curve.num_eq
        df_dz = f.diff(z)
        df_dx = f.diff(x)
        # F = -(\partial f/\partial z)/(\partial f/\partial x)
        F = sympy.lambdify((z, x), -df_dz/df_dx)
        v = sympy.lambdify((z, x), self.sw_diff.num_v)

        def ode_f(t, zx1x2):
            z_i = zx1x2[0]
            x1_i = zx1x2[1]
            x2_i = zx1x2[2]
            dz_i_dt = exp(self.theta*1j)/(v(z_i, x1_i) - v(z_i, x2_i))
            dx1_i_dt = F(z_i, x1_i) * dz_i_dt
            dx2_i_dt = F(z_i, x2_i) * dz_i_dt
            return [dz_i_dt, dx1_i_dt, dx2_i_dt]

        self.ode = integrate.ode(ode_f)
        self.ode.set_integrator(
            'zvode',
            # method='adams',
            atol=ode_absolute_tolerance,
            # min_step = ode_min_step,
            # max_step = ode_max_step,
        )

    def get_new_joints(self, current_s_wall_index, previous_length):
        """
        Find joints between the newly grown segment of the given S-wall
        and the other S-walls. This checks joints that are formed by two
        S-walls only, not considering the possibility of a joint of three
        S-walls, which in principle can happen but is unlikely in a numerical
        setup.
        """
        i_c = current_s_wall_index
        t_0 = previous_length
        new_joints = []
        new_curve = self.s_walls[i_c].get_zxzys(t_0-1)
        new_bin_keys = self.hit_table.fill(i_c, new_curve, t_0-1)

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
                for t_c_i, t_c_f in self.hit_table[bin_key][i_c]:
                    if t_0 > t_c_i:
                        # this segment is not a new one.
                        continue
                    for t_d_i, t_d_f in self.hit_table[bin_key][i_d]:
                        segment_c = self.s_walls[i_c].get_zxzys(t_c_i, t_c_f+1)
                        segment_d = self.s_walls[i_d].get_zxzys(t_d_i, t_d_f+1)

                        try:
                            # find an intersection on the z-plane.
                            ip_x, ip_y = find_intersection_of_segments(
                                segment_c,
                                segment_d,
                                self.hit_table.get_bin_location(bin_key),
                                self.hit_table.get_bin_size()
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
                            ip_xs = find_xs_at_z_0(self.sw_curve.num_eq, ip_z)
                            ip_x1_c = n_nearest(ip_xs, x1_c, 1)[0]
                            ip_x2_c = n_nearest(ip_xs, x2_c, 1)[0]
                            ip_x1_d = n_nearest(ip_xs, x1_d, 1)[0]
                            ip_x2_d = n_nearest(ip_xs, x2_d, 1)[0]
                            
                            a_joint_label = ('joint ' +
                                             '#{}'.format(len(self.joints)))
                            a_joint = get_joint(ip_z,
                                                ip_x1_c, ip_x2_c,
                                                ip_x1_d, ip_x2_d,
                                                self.s_walls[i_c],
                                                self.s_walls[i_d],
                                                self.config_data.accuracy,
                                                label=a_joint_label)

                            if(a_joint is None):
                                continue
                            else:
                                self.joints.append(a_joint)
                                new_joints.append(a_joint)
                            
                        except NoIntersection:
                            pass

        return new_joints

    def grow(self):
        rpzs = []
        ppzs = []

        for rp in self.sw_curve.ramification_points:
            if rp.is_puncture is False:
                rpzs.append(rp.z)
            elif rp.is_puncture is True:
                ppzs.append(rp.z)
        for iteration in range(self.config_data.num_of_iterations):
            new_joints = []
            for i in range(len(self.s_walls)):
                if self.s_walls[i].out_of_range(
                    self.config_data.z_range_limits
                ):
                    continue
                previous_length = self.s_walls[i].grow(
                    self.ode, rpzs, ppzs, 
                    z_range_limits=self.config_data.z_range_limits,
                    num_of_steps=self.config_data.num_of_steps,
                    size_of_small_step=self.config_data.size_of_small_step,
                    size_of_large_step=self.config_data.size_of_large_step,
                    size_of_neighborhood=self.config_data.size_of_neighborhood,
                )
                new_joints += self.get_new_joints(i, previous_length)

            for joint in new_joints:
                label = 'S-wall #{}'.format(len(self.s_walls))
                self.s_walls.append(
                    SWall(joint.z, joint.x1, joint.x2, joint.parents, label)
                )
                

def get_local_sw_diff(sw_curve, sw_diff, ramification_point):
    rp = ramification_point
    local_curve = sw_curve.num_eq.series(x, rp.x, rp.i+1)
    local_curve = local_curve.series(z, rp.z, 2).removeO()
    # curve_at_rp = a(z - rp.z) + b(x - rp.x)^(rp.i)
    # translate z such that rp.z = 0
    a = local_curve.subs(z, z+rp.z).coeff(z).coeff(x, 0)
    # translate x such that rp.x = 0
    b = local_curve.subs(x, x+rp.x).coeff(x**rp.i).coeff(z, 0)
    local_x = rp.x + (-(a/b)*(z - rp.z))**sympy.Rational(1, rp.i)
    # substitute x with x(z)
    local_diff = sw_diff.num_v.subs(x, local_x)
    # series expansion in z at rp.z
    local_diff = local_diff.series(z, rp.z, 1).subs(z, z+rp.z)
    # get the coefficient and the exponent of the leading term
    (diff_c, diff_e) = local_diff.leadterm(z)
    if diff_e == 0:
        # remove the constant term from the local_diff
        local_diff -= local_diff.subs(z, 0)
        (diff_c, diff_e) = local_diff.leadterm(z)

    return (complex(diff_c), diff_e)


def get_s_wall_seeds(sw_curve, sw_diff, theta, ramification_point,
                     config_data,):
    rp = ramification_point
    delta = config_data.accuracy
    dt = config_data.size_of_small_step

    ###
    # 1. Find the first-order approximations of the starting points 
    # of S-walls around a given branch point, which is of the form
    # \Delta z_i = c_i / (\lambda_0)^(rp.i/(rp.i+1))
    # at a branch point from a ramification point, and
    # \Delta z_i = c_i / (\lambda_0)^(rp.i) exp(rp.i theta I)
    # at a branch point from a massless regular puncture
    ###

    # 1.1 find the coefficient and the exponent of the leading term
    # of the SW differential at the ramification point.
    lambda_0, diff_e = get_local_sw_diff(sw_curve, sw_diff, rp)

    # 1.2 find c_i, a phase factor for each S-wall.
    omega_1 = exp(2*pi*1j/rp.i)
    omega = [omega_1**k for k in range(rp.i)]

    beta_1 = exp(rp.i*2*pi*1j/(rp.i+1))
    beta = [beta_1**k for k in range(rp.i+1)]

    cs = []
    if diff_e == sympy.Rational(1, rp.i):
        # the branch point is a ramification point
        # go over pairs of omegas that differ by \omega_1^i
        for i in range(1, rp.i):
            new_locs = []
            # and go over all the omegas
            for j in range(rp.i): 
                if j+i < rp.i:
                    new_loc = 1/(omega[j]-omega[j+i])
                else:
                    new_loc = 1/(omega[j]-omega[j+i-rp.i])
                new_loc = cpow(new_loc, rp.i, rp.i+1)
                new_locs += [new_loc*beta_i for beta_i in beta]
            new_locs = remove_duplicate(new_locs, 
                                        lambda l1, l2: abs(l1 - l2) < delta)
            cs += [(c*exp(theta*1j*sympy.Rational(rp.i, rp.i+1)) /
                    cpow(lambda_0, rp.i, rp.i+1)) for c in new_locs]

    elif diff_e == -1 + sympy.Rational(1, rp.i):
        # the branch point is a massless regular puncture
        # go over pairs of omegas that differ by \omega_1^i
        for i in range(1, rp.i):
            cs.append(exp(rp.i*theta*1j) /
                      cpow(((omega[0]-omega[i])*lambda_0), rp.i))

    else:
        logging.error('unknown form of sw_diff at rp ({}, {}): '
                      'diff_e  = {}'.format(rp.z, rp.x, diff_e))
        raise GetSWallSeedsError(diff_e)
        
    cs = gather(cs, lambda c1, c2: abs(c1 - c2) < delta)
    logging.debug('list of c = %s, # = %d', cs, len(cs))

    # 2. Now calculate \Delta z_i for each S-wall and 
    # find the two points on the curve that are projected onto it.  
    seeds = []
    for c in cs:
        cv = c[0]   # value of c
        cm = c[1]   # multiplicity of c
        # resize to the size of the small step 
        Delta_z = cv/abs(cv)*dt
        z_0 = rp.z + Delta_z
        xs_at_z_0 = find_xs_at_z_0(sw_curve.num_eq, z_0, rp.x, rp.i)
        dev_phases = [pi for i in range(len(xs_at_z_0)**2)] 
        for i in range(len(xs_at_z_0)):
            diffx = sw_diff.num_v.subs(z, z_0) 
            v_i = complex(diffx.subs(x, xs_at_z_0[i]))
            for j in range(len(xs_at_z_0)):
                if i == j:
                    continue
                else:
                    v_j = complex(diffx.subs(x, xs_at_z_0[j])) 
                    delta_z = exp(1j*theta)/(v_i - v_j)*dt
                    # flattened index
                    fij = i*len(xs_at_z_0) + j
                    dev_phases[fij] = phase((delta_z/Delta_z))
        min_dev_indices = n_nearest_indices(dev_phases, 0.0, cm)
        for k in min_dev_indices:
            i, j = unravel(k, len(xs_at_z_0))
            seeds.append((z_0, xs_at_z_0[i], xs_at_z_0[j]))
        
    return seeds


def get_joint(z, x1_i, x2_i, x1_j, x2_j, parent_i, parent_j, accuracy, 
              spectral_network_type='A', label=None):
    """
    return a joint if formed, otherwise return None.
    """
    # TODO: implementation of joint rule changes according to the spectral
    # network type.
    if(abs(x1_i - x2_j) < accuracy and abs(x1_j - x2_i) < accuracy):
        return None
    elif(abs(x2_i - x1_j) < accuracy):
        return Joint(z, x1_i, x2_j, [parent_i, parent_j], label)
    elif(abs(x2_j - x1_i) < accuracy):
        return Joint(z, x1_j, x2_i, [parent_j, parent_i], label)


def init_process():
    """
    Initializer of each process that generates a spectral network.
    Take care of a keyboard interrupt.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def parallel_get_spectral_network(sw_curve, sw_diff, phase, config_data,
                                  shared_n_started_spectral_networks,
                                  shared_n_finished_spectral_networks,):
    shared_n_started_spectral_networks.value += 1
    logging.info('Start generating spectral network #{}.'.format(
        shared_n_started_spectral_networks.value
    ))

    spectral_network = SpectralNetwork(sw_curve, sw_diff, phase, config_data) 
    spectral_network.grow()

    shared_n_finished_spectral_networks.value += 1
    logging.info('Finished generating spectral network #{}.'.format(
        shared_n_finished_spectral_networks.value
    ))

    return spectral_network

def generate_spectral_network(opts, config_data):
    sw_curve = SWCurve(config_data)
    sw_curve.find_ramification_points()
    sw_diff = SWDiff(config_data)

    logging.info('\nList of ramification points')
    for rp in sw_curve.ramification_points:
        logging.info('%s', rp)
    logging.info('\n')

    if(opts['phase'] is not None):
        spectral_network = SpectralNetwork(sw_curve, sw_diff, opts['phase'], 
                                           config_data) 

        spectral_network.grow()

        if(opts['show-plot'] is True):
            spectral_network_plot = SpectralNetworkPlot(
                spectral_network,
                # plot_data_points=True,
            )
            spectral_network_plot.show()

    elif(config_data.phase_range is not None):
        theta_i, theta_f, theta_n = config_data.phase_range
        phases = [theta_i + i * (theta_f - theta_i) / theta_n
                  for  i in range(theta_n)]
        spectral_networks = []

        manager = multiprocessing.Manager()
        shared_n_started_spectral_networks = manager.Value('i', 0)
        shared_n_finished_spectral_networks = manager.Value('i', 0)
        if(config_data.n_processes == 0):
            n_processes = multiprocessing.cpu_count()
        logging.info('Number of processes: {}'.format(n_processes))
        pool =  multiprocessing.Pool(n_processes, init_process)
        try:
            results = [
                pool.apply_async(
                    parallel_get_spectral_network, 
                    args=(
                        sw_curve, sw_diff, phase, config_data,
                        shared_n_started_spectral_networks,
                        shared_n_finished_spectral_networks,
                    )
                ) for phase in phases
            ]
            pool.close()

            for result in results:
                spectral_networks.append(result.get())
                logging.info('job progress: {}/{} finished.'.format(
                    shared_n_finished_spectral_networks.value,
                    theta_n
                ))

        except KeyboardInterrupt:
            logging.warning('Caught ^C; terminates processes...')
            pool.terminate()
            pool.join()

        pdb.set_trace()

