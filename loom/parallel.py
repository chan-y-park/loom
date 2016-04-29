import signal
import multiprocessing
import logging
import os

from spectral_network import SpectralNetwork

#from spectral_network import SpectralNetwork


def child_sigint_handler(signum, frame):
    # print("SIGINT catched by a child of loom.parallel.")
    raise KeyboardInterrupt("SIGINT catched by a child of loom.parallel.")


def child_sigterm_handler(signum, frame):
    print("SIGTERM catched by a child of loom.parallel.")
    raise SystemExit("SIGTERM catched by a child of loom.parallel.")


def init_process():
    """
    Initializer of each child process that generates a spectral network.
    Take care of a keyboard interrupt and a SIGTERM.
    """
    signal.signal(signal.SIGINT, child_sigint_handler)
    # signal.signal(signal.SIGTERM, child_sigterm_handler)


def a_child_process(
    config=None,
    sw_data=None,
    spectral_network=None,
    job_id=None,
    n_jobs=None,
    logger_name='loom',
    cache_file_path=None,
    #shared_n_started_spectral_networks,
    #shared_n_finished_spectral_networks,
    additional_n_steps=0,
    new_mass_limit=None,
    additional_iterations=0,
):
    logger = logging.getLogger(logger_name)

    logger.info('Start generating spectral network #{}/{}: phase = {}.'
                .format(job_id, n_jobs, spectral_network.phase))
   
    try:
        spectral_network.grow(
            config, sw_data,
            additional_iterations=additional_iterations,
            additional_n_steps=additional_n_steps,
            new_mass_limit=new_mass_limit,
            cache_file_path=cache_file_path,
        )
    except Exception as e:
        error_msg = (
            'A child process calculating phase = {} caught an exception: {}'
            .format(spectral_network.phase, e)
        )
        logger.warning(error_msg)
        spectral_network.errors.append = ('Unknown', error_msg)

    logger.info('Finished generating spectral network #{}/{}.'
                .format(job_id, n_jobs))

    if cache_file_path is None:
        return spectral_network
    else:
#        cache_file_path = os.path.join(
#            cache_dir,
#            'data_{}.json'.format(
#                str(job_id).zfill(len(str(n_jobs - 1)))
#            )
#        )
#        logger.info('Saving cache data to {}.'.format(cache_file_path))
#        spectral_network.save(cache_file_path)
        return cache_file_path


def parallel_get_spectral_network(
    config=None,
    sw_data=None, 
    spectral_networks=None,
    n_processes=0,
    additional_n_steps=0,
    new_mass_limit=None,
    additional_iterations=0,
    logger_name='loom',
    cache_dir=None,
    data_file_prefix='data',
):
    logger = logging.getLogger(logger_name)
    phase = config['phase']
    if spectral_networks is None and isinstance(phase, list) is True:
        theta_i, theta_f, theta_n = phase
        phases = [
            (float(theta_i) + i * float(theta_f - theta_i) / (theta_n - 1))
            for i in range(theta_n)
        ]
        spectral_networks = [
            SpectralNetwork(
                phase=sn_phase, 
                logger_name=logger_name,
            ) 
            for sn_phase in phases
        ]

    n_cpu = multiprocessing.cpu_count()
    if (n_processes == 0):
        # Use all the CPUs.
        n_processes = n_cpu
    elif (n_processes < 0):
        # Leave |n_processes| CPUs.
        if(n_cpu > -n_processes):
            n_processes = n_cpu - (-n_processes)
        else:
            logger.warning('The number of CPUs is smaller than {}.'
                           .format(-n_processes))
            logger.warning('Set n_processes to 1.')
            n_processes = 1
    elif (n_cpu < n_processes):
            logger.warning('The number of CPUs is smaller than {}.'
                           .format(n_processes))
            logger.warning('Set n_processes to {}.'.format(n_cpu))
            n_processes = n_cpu

    # Use n_processes CPUs.
    n_jobs = len(spectral_networks)
    if n_jobs < n_processes:
        n_processes = n_jobs
    multiprocessing.freeze_support()
    pool = multiprocessing.Pool(n_processes, init_process)
    logger.info('Number of processes in the pool: {}'.format(n_processes))


    try:
        results = [
            pool.apply_async(
                a_child_process,
                kwds=dict(
                    config=config,
                    sw_data=sw_data,
                    spectral_network=sn,
                    job_id=i+1,
                    n_jobs=n_jobs,
                    logger_name=logger_name,
                    additional_n_steps=additional_n_steps,
                    new_mass_limit=new_mass_limit,
                    additional_iterations=additional_iterations,
                    cache_file_path=get_cache_file_path(
                        cache_dir, data_file_prefix, i,
                    )
                )
            ) for i, sn in enumerate(spectral_networks)
        ]
        pool.close()
        pool.join()

        new_spectral_networks = []
        for result in results:
            if cache_dir is None:
                new_sn = result.get()
            else:
                # FIXME: Do we really need the following?
                new_sn_file_path = result.get() 
                new_sn = SpectralNetwork(logger_name=logger_name)
                new_sn.load(new_sn_file_path, sw_data)
            new_spectral_networks.append(new_sn)

    except (KeyboardInterrupt, SystemExit) as e:
        logger.warning('loom.parallel caught {}: {}; terminates processes...'
                       .format(type(e), e.args,))
        pool.terminate()
        pool.join()
        raise e
    except Exception as e:
        logger.warning('loom.parallel caught {}.'.format(e))
        raise e

    return new_spectral_networks


def get_cache_file_path(cache_dir, data_file_prefix, i):
    if cache_dir is not None and data_file_prefix is not None:
        cache_file_path=os.path.join(
            cache_dir,
            '{}_{}.json'.format(data_file_prefix, i)
        )
        return cache_file_path
    else:
        return None
