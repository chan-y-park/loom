import os
import signal
import multiprocessing
import logging

def init_process():
    """
    Initializer of each child process that generates a spectral network.
    Take care of a keyboard interrupt.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def a_child_process(
    sw,
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

    spectral_network.grow(sw, config)

    shared_n_finished_spectral_networks.value += 1
    logging.info('Finished generating spectral network #{}.'.format(
        shared_n_finished_spectral_networks.value
    ))

    spectral_network_data = spectral_network.get_data()

    # Save spectral network data to a file
    data_file_name = os.path.join(data_save_dir, 'data_{}.json'.format(job_id))
    logging.info('Saving data to {}.'.format(data_file_name))
    with open(data_file_name, 'wb') as fp:
        save_spectral_network_data(
            spectral_network_data, 
            fp, 
        )

    logging.info('job progress: {}/{} finished.'
                 .format(shared_n_finished_spectral_networks.value, theta_n))

    return spectral_network_data


def parallel_get_spectral_network(
    sw, 
    ramification_points, 
    config,
    data_save_dir,
):
    spectral_network_data_list = []
    n_processes = config['n_processes']

    theta_i, theta_f, theta_n = config['phase_range']
    phases = [(theta_i + i * (theta_f - theta_i) / theta_n)
              for  i in range(theta_n)]

    manager = multiprocessing.Manager()
    shared_n_started_spectral_networks = manager.Value('i', 0)
    shared_n_finished_spectral_networks = manager.Value('i', 0)

    if(n_processes == 0):
        n_processes = multiprocessing.cpu_count()
    elif(n_processes < 0):
        n_cpu = multiprocessing.cpu_count()
        if(n_cpu > -n_processes):
            n_processes = n_cpu - (-n_processes)
        else:
            logging.warning('The number of CPUs is smaller than {}.'
                            .format(-config['n_processes']))
            logging.warning('Set n_processes to 1.')
            n_processes = 1

    pool =  multiprocessing.Pool(n_processes, init_process)
    logging.info('Number of processes in the pool: {}'.format(n_processes))

    try:
        results = [
            pool.apply_async(
                a_child_process, 
                args=(
                    sw,
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

    except KeyboardInterrupt:
        logging.warning('Caught ^C; terminates processes...')
        pool.terminate()
        pool.join()

    return spectral_network_data_list
