import signal
import multiprocessing
import logging
import os

from spectral_network import SpectralNetwork
from spectral_network import SolitonTree
from geometry import BranchPoint


def child_sigint_handler(signum, frame):
    raise KeyboardInterrupt("SIGINT catched by a child of loom.parallel.")


def child_sigterm_handler(signum, frame):
    raise SystemExit("SIGTERM catched by a child of loom.parallel.")


def init_process():
    """
    Initializer of each child process that generates a spectral network.
    Take care of a keyboard interrupt and a SIGTERM.
    """
    signal.signal(signal.SIGINT, child_sigint_handler)


def get_spectral_network_child_process(
    config=None,
    sw_data=None,
    spectral_network=None,
    job_id=None,
    n_jobs=None,
    logger_name='loom',
    cache_file_path=None,
    additional_n_steps=0,
    new_mass_limit=None,
    additional_iterations=0,
    z_plane_rotation=None,
    two_way_streets_only=False,
    search_radius=None,
    task='grow',
    method=None,
    downsample=False,
    downsample_ratio=None,
):
    logger = logging.getLogger(logger_name)

    if task == 'grow':
        msg = 'growing'
    elif task == 'trivialize':
        msg = 'trivializing'
    elif task == 'find_two_way_streets':
        msg = 'finding two-way streets of'

    logger.info('Start {} spectral network #{}/{}: phase = {}.'
                .format(msg, job_id, n_jobs - 1, spectral_network.phase))
    try:
        if task == 'grow':
            spectral_network.grow(
                config, sw_data,
                additional_iterations=additional_iterations,
                additional_n_steps=additional_n_steps,
                new_mass_limit=new_mass_limit,
                cache_file_path=cache_file_path,
                method=method,
                downsample=downsample,
                downsample_ratio=downsample_ratio,
            )
        elif task == 'trivialize':
            spectral_network.trivialize(
                config, sw_data,
                cache_file_path=cache_file_path,
                z_plane_rotation=z_plane_rotation,
                two_way_streets_only=two_way_streets_only,
            )
        elif task == 'find_two_way_streets':
            spectral_network.find_two_way_streets(
                config=config,
                sw_data=sw_data,
                search_radius=search_radius,
                cache_file_path=cache_file_path,
            )
    except Exception as e:
        error_msg = (
            'A child process {} spectral network #{} @ phase = {} '
            'caught an exception: {}'
            .format(msg, job_id, spectral_network.phase, e)
        )
        logger.warning(error_msg)
        spectral_network.errors.append = ('Unknown', error_msg)

    logger.info('Finished {} spectral network #{}/{}.'
                .format(msg, job_id, n_jobs - 1))

    if cache_file_path is None:
        return spectral_network
    else:
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
    z_plane_rotation=None,
    two_way_streets_only=False,
    search_radius=None,
    data_file_prefix='data',
    task='grow',
    method=None,
    downsample=False,
    downsample_ratio=None,
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

    n_jobs = len(spectral_networks)
    n_processes = get_n_processes(n_processes, n_jobs, logger_name)

    multiprocessing.freeze_support()
    pool = multiprocessing.Pool(n_processes, init_process)
    logger.info('Number of processes in the pool: {}'.format(n_processes))

    try:
        results = [
            pool.apply_async(
                get_spectral_network_child_process,
                kwds=dict(
                    config=config,
                    sw_data=sw_data,
                    spectral_network=sn,
                    job_id=i,
                    n_jobs=n_jobs,
                    logger_name=logger_name,
                    additional_n_steps=additional_n_steps,
                    new_mass_limit=new_mass_limit,
                    additional_iterations=additional_iterations,
                    cache_file_path=get_cache_file_path(
                        cache_dir, data_file_prefix, i,
                    ),
                    z_plane_rotation=z_plane_rotation,
                    search_radius=search_radius,
                    two_way_streets_only=two_way_streets_only,
                    task=task,
                    method=method,
                    downsample=downsample,
                    downsample_ratio=downsample_ratio,
                )
            ) for i, sn in enumerate(spectral_networks)
        ]
        pool.close()
        pool.join()

        new_spectral_networks = []
        for result in results:
            new_sn = SpectralNetwork(logger_name=logger_name)
            if cache_dir is None:
                new_sn.set_from_json_data(
                    result.get().get_json_data(),
                    sw_data,
                )
                # NOTE: Need to replace all the branch points
                # in a spectral network with those in sw_data,
                # because a child process gets a copy of them.
#               new_sn = result.get()
#                for s_wall in new_sn.s_walls:
#                    if len(s_wall.parents) == 1:
#                        if isinstance(s_wall.parents[0], BranchPoint):
#                            bp_label = s_wall.parents[0].label
#                            s_wall.parents = [
#                                bp for bp in sw_data.branch_points
#                                if bp.label == bp_label
#                            ]
            else:
                # XXX: Do we really need the following?
                # It's probably good to use cache files
                # just in case the internal queue size
                # of the pool is too small.
                new_sn_file_path = result.get()
                #new_sn = SpectralNetwork(logger_name=logger_name)
                new_sn.load(new_sn_file_path, sw_data)
            new_spectral_networks.append(new_sn)
    except Exception as e:
        logger.warning(
            'parallel_get_spectral_network() caught {}: {}; '
            .format(type(e), e.args,)
        )
        if e == KeyboardInterrupt or e == SystemExit:
            logger.warning('Terminates children processes...')
            pool.terminate()
            pool.join()
        raise e
#    except (KeyboardInterrupt, SystemExit) as e:
#        logger.warning(
#            'parallel_get_spectral_network() caught {}: {}; '
#            'terminates children processes...'
#            .format(type(e), e.args,)
#        )
#        pool.terminate()
#        pool.join()
#        raise e
#    except Exception as e:
#        logger.warning(
#            'parallel_get_spectral_network() caught {}.'
#            .format(e)
#        )
#        raise e

    return new_spectral_networks


def get_improved_soliton_tree_child_process(
    config=None,
    sw_data=None,
    search_radius=None,
    tree=None,
    job_id=None,
    n_jobs=None,
    cache_file_path=None,
    label=None,
    logger_name='loom',
):
    logger = logging.getLogger(logger_name)
    logger.info('Improve soliton tree #{}/{}: phase = {}.'
                .format(job_id, n_jobs - 1, tree.phase))
    try:
        improved_tree = tree.get_improved_tree(
            config=config,
            sw_data=sw_data,
            search_radius=search_radius,
            label=label,
            logger_name=logger_name,
            cache_file_path=cache_file_path,
        )
    except Exception as e:
        logger.warning(
            'A child process improving soliton tree #{} @ phase = {} '
            'caught an exception: {}'
            .format(job_id, tree.phase, e)
        )

    logger.info('Finished improving soliton tree #{}/{}.'
                .format(job_id, n_jobs - 1))

    if cache_file_path is None:
        return improved_tree
    else:
        return cache_file_path
        

def parallel_get_improved_soliton_tree(
    config=None,
    sw_data=None,
    trees=[],
    search_radius=None,
    n_processes=0,
    cache_dir=None,
    logger_name='loom',
):
    logger = logging.getLogger(logger_name)

    n_jobs = len(trees)
    n_processes = get_n_processes(n_processes, n_jobs, logger_name)

    multiprocessing.freeze_support()
    pool = multiprocessing.Pool(n_processes, init_process)
    logger.info('Number of processes in the pool: {}'.format(n_processes))

    try:
        results = [
            pool.apply_async(
                get_improved_soliton_tree_child_process,
                kwds=dict(
                    config=config,
                    sw_data=sw_data,
                    search_radius=search_radius,
                    tree=tree,
                    job_id=i,
                    n_jobs=n_jobs,
                    cache_file_path=get_cache_file_path(
                        cache_dir, 'tree', i,
                    ),
                    logger_name=logger_name,
                ),
            ) for i, tree in enumerate(trees)
        ]
        pool.close()
        pool.join()

        improved_trees = []

        for result in results:
            if cache_dir is None:
                improved_tree = result.get()
                # NOTE: Need to replace all the branch points
                # in a soliton tree with those in sw_data
                # because a child process gets a copy of each.
                improved_tree.root_branch_point = sw_data[
                    improved_tree.root_branch_point.label
                ]
                for street in improved_tree.streets:
                    if (len(street.parents) == 1):
                        parent = street.parents[0]
                        if isinstance(parent, BranchPoint):
                            street.parents = [sw_data[parent.label]]
            else:
                # XXX: Do we really need the following?
                # It's probably good to use cache files
                # just in case the internal queue size
                # of the pool is too small.
                improved_tree_file_path = result.get()
                #improved_tree = SolitonTree(logger_name=logger_name)
                improved_tree = SolitonTree()
                improved_tree.load(improved_tree_file_path, sw_data)

            improved_trees.append(improved_tree)
    except Exception as e:
        logger.warning(
            'parallel_get_improved_soliton_tree() caught {}: {}; '
            .format(type(e), e.args,)
        )
        if e == KeyboardInterrupt or e == SystemExit:
            logger.warning('Terminates children processes...')
            pool.terminate()
            pool.join()
        raise e

    return improved_trees

def get_n_processes(
    n_processes=None,
    n_jobs=None,
    logger_name='loom',
):
    logger = logging.getLogger(logger_name)
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

    if n_jobs < n_processes:
        n_processes = n_jobs

    return n_processes


def get_cache_file_path(cache_dir, data_file_prefix, i):
    if cache_dir is not None and data_file_prefix is not None:
        cache_file_path = os.path.join(
            cache_dir,
            '{}_{}.json'.format(data_file_prefix, i)
        )
        return cache_file_path
    else:
        return None
