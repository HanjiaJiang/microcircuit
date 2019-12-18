import os
import sys
import pickle
import multiprocessing as mp
import time
import microcircuit.network as network
import microcircuit.tools as tools
from microcircuit.create_params import params_single

if __name__ == "__main__":
    # simulation settings
    run_sim = True
    on_server = False

    # analysis settings
    plot_half_len = 200.0   # ms
    fr_interval = [1000.0, 2000.0]
    sf_interval = 20  # ms

    # check for: parameter scan or single-run
    try:
        pickle_path = sys.argv[1]    # path to pickle file
        # os.system('cp Snakefile run_network.py ' + os.path.dirname(pickle_path))
        # os.mkdir(os.path.join(os.path.dirname(pickle_path), 'microcircuit'))
        # os.system('cp microcircuit/*.py ' + os.path.join(os.path.dirname(pickle_path), 'microcircuit'))
    except IndexError:  # single-run if no path input
        print('No argv[1]; single-run.')
        cwd = os.getcwd()
        data_path = os.path.join(cwd, 'data')
        if not os.path.isdir(data_path):
            os.mkdir(data_path)
        pickle_path = os.path.join(cwd, 'para_dict.pickle')
        params_single(pickle_path)
        os.system('cp run_network.py ' + data_path)
        os.mkdir(os.path.join(data_path, 'microcircuit'))
        os.system('cp microcircuit/*.py ' + os.path.join(data_path, 'microcircuit'))
    # assign parameters
    with open(pickle_path, 'rb') as handle:
        para_dict = pickle.load(handle)

    # cpu number / cluster or not
    cpu_ratio = 0.5
    if mp.cpu_count() > 10:
        cpu_ratio = 1
        on_server = True
    else:
        # exec(tools.set2txt())
        pass
    para_dict['sim_dict']['local_num_threads'] = \
        int(mp.cpu_count() * cpu_ratio)

    # run simulation
    net = network.Network(para_dict['sim_dict'], para_dict['net_dict'],
                          para_dict['stim_dict'], para_dict['special_dict'])
    net.setup()
    if run_sim:
        net.simulate()

    # analysis
    tools.plot_raster(
        para_dict['sim_dict']['data_path'], 'spike_detector',
        para_dict['stim_dict']['th_start'][0] - plot_half_len,
        para_dict['stim_dict']['th_start'][0] + plot_half_len)
    mean_fr_cache, std_fr = \
        tools.fire_rate(para_dict['sim_dict']['data_path'], 'spike_detector',
                        fr_interval[0], fr_interval[1])
    # tools.response(para_dict['sim_dict']['data_path'], 'spike_detector',
    #                para_dict['stim_dict']['th_start'][0], sf_interval)
    tools.boxplot(para_dict['net_dict'], para_dict['sim_dict']['data_path'])

    t0 = time.time()
    tools.ai_score(para_dict['sim_dict']['data_path'], 'spike_detector',
        fr_interval[0], fr_interval[1])
    print('ai analysis time = {}'.format(time.time() - t0))

    # delete .gdf files to save space
    if os.path.isdir(para_dict['sim_dict']['data_path']):
        for item in os.listdir(para_dict['sim_dict']['data_path']):
            if item.endswith('.gdf'):
                os.remove(os.path.join(para_dict['sim_dict']['data_path'], item))

    if not on_server:
        # exec(tools.end2txt())
        pass



