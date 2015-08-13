import os
import signal
import multiprocessing
import logging

from spectral_network import SpectralNetwork

def init_process():
    """
    Initializer of each child process that generates a spectral network.
    Take care of a keyboard interrupt.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def a_child_process(
    sw,
    phase,
    config,
    shared_n_started_spectral_networks,
    shared_n_finished_spectral_networks
):
    theta_i, theta_f, theta_n = config['phase_range']

    shared_n_started_spectral_networks.value += 1
    job_id = shared_n_started_spectral_networks.value
    logging.info('Start generating spectral network #{}/{}: theta = {}.'
                 .format(job_id, theta_n, phase)
    )

    spectral_network = SpectralNetwork(
        phase, 
    ) 

    spectral_network.grow(config, sw)

    shared_n_finished_spectral_networks.value += 1
    logging.info('Finished generating spectral network #{}/{}.'
                 .format(shared_n_finished_spectral_networks.value, theta_n)
    )

    return spectral_network


def parallel_get_spectral_network(
    sw, 
    config,
):
    spectral_network_list = []
    n_processes = config['n_processes']

    theta_i, theta_f, theta_n = config['phase_range']
    phases = [(float(theta_i) + i * float(theta_f - theta_i) / (theta_n-1))
              for i in range(theta_n)]

    manager = multiprocessing.Manager()
    shared_n_started_spectral_networks = manager.Value('i', 0)
    shared_n_finished_spectral_networks = manager.Value('i', 0)

    n_cpu = multiprocessing.cpu_count()
    if (n_processes == 0):
        ### Use all the CPUs.
        n_processes = n_cpu
    elif (n_processes < 0):
        ### Leave |n_processes| CPUs.
        if(n_cpu > -n_processes):
            n_processes = n_cpu - (-n_processes)
        else:
            logging.warning('The number of CPUs is smaller than {}.'
                            .format(-config['n_processes']))
            logging.warning('Set n_processes to 1.')
            n_processes = 1
    elif (n_cpu < n_processes):
            logging.warning('The number of CPUs is smaller than {}.'
                            .format(config['n_processes']))
            logging.warning('Set n_processes to {}.'.format(n_cpu))
            n_processes = n_cpu

    ### Use n_processes CPUs.
    multiprocessing.freeze_support()
    pool =  multiprocessing.Pool(n_processes, init_process)
    logging.info('Number of processes in the pool: {}'.format(n_processes))

    try:
        results = [
            pool.apply_async(
                a_child_process,
                args=(
                    sw,
                    phase,
                    config,
                    shared_n_started_spectral_networks,
                    shared_n_finished_spectral_networks,
                )
            ) for phase in phases
        ]
        pool.close()

        for result in results:
            spectral_network_list.append(result.get())

    except KeyboardInterrupt:
        logging.warning('Caught ^C; terminates processes...')
        pool.terminate()
        pool.join()

    return spectral_network_list
