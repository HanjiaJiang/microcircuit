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
    run_analysis = True
    run_ai = False
    run_response = True

    # timing, in ms
    plot_half_len = 100.0
    start = 2000.0
    segment = 2000.0
    n_segment = 2
    total_length = segment*n_segment
    interval = [start, start + total_length]

    # thalamic input
    th_starts = list(range(int(start), int(start + total_length), int(segment)))
    th_rate = 300.0

    # check for: parameter scan or single-run
    try:
        pickle_path = sys.argv[1]    # path to pickle file
    except IndexError:  # single-run if no path input
        print('No argv[1]; single-run.')
        cwd = os.getcwd()

        # create pickle file
        pickle_path = os.path.join(cwd, 'para_dict.pickle')
        params_single(pickle_path, th_starts=th_starts, th_rate=th_rate)

        # handle data path
        data_path = os.path.join(cwd, 'data')
        if not os.path.isdir(data_path):
            os.mkdir(data_path)
        os.system('cp run_network.py ' + data_path)

        # handle microcircuit files
        m_path = os.path.join(data_path, 'microcircuit')
        if not os.path.isdir(m_path):
            os.mkdir(m_path)
        os.system('cp microcircuit/*.py ' + os.path.join(data_path, 'microcircuit'))

    # assign parameters
    with open(pickle_path, 'rb') as handle:
        para_dict = pickle.load(handle)
    handle.close()

    # cpu number / on server or not
    cpu_ratio = 0.5
    if mp.cpu_count() > 10:
        cpu_ratio = 1
        on_server = True
    exec(tools.set2txt(para_dict['sim_dict']['data_path']))

    # set dictionary
    para_dict['sim_dict']['local_num_threads'] = \
        int(mp.cpu_count() * cpu_ratio)
    para_dict['sim_dict']['t_sim'] = start + total_length

    # run simulation
    net = network.Network(para_dict['sim_dict'], para_dict['net_dict'],
                          para_dict['stim_dict'], para_dict['special_dict'])
    net.setup()
    if run_sim:
        net.simulate()

    # analysis
    if run_analysis:
        spikes = tools.Spikes(para_dict['sim_dict']['data_path'], 'spike_detector')
        spikes.load_data_all()
        mean_fr, std_fr = \
            tools.fire_rate(para_dict['sim_dict']['data_path'], 'spike_detector',
                            interval[0], interval[1])
        if run_ai:
            t0 = time.time()
            tools.ai_score(para_dict['sim_dict']['data_path'], 'spike_detector',
                           interval[0], interval[1], seg_len=segment)
            print('ai analysis time = {}'.format(time.time() - t0))
        if run_response:
            t1 = time.time()
            tools.response(spikes, para_dict['sim_dict']['data_path'],
                           para_dict['stim_dict']['th_start'][0],
                           window=20.0,
                           n_stim=n_segment,
                           interval=segment)
            print('response analysis time = {}'.format(time.time() - t1))
        tools.plot_raster(
            para_dict['sim_dict']['data_path'], 'spike_detector',
            para_dict['stim_dict']['th_start'][-1] - plot_half_len,
            para_dict['stim_dict']['th_start'][-1] + plot_half_len)
        tools.fr_boxplot(para_dict['net_dict'], para_dict['sim_dict']['data_path'])

    # delete .gdf files to save space
    if on_server and os.path.isdir(para_dict['sim_dict']['data_path']):
        for item in os.listdir(para_dict['sim_dict']['data_path']):
            if item.endswith('.gdf'):
                os.remove(os.path.join(para_dict['sim_dict']['data_path'], item))

    exec(tools.end2txt())
