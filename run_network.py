import os
import sys
import pickle
import multiprocessing as mp
import time
import microcircuit.network as network
import microcircuit.tools as tools
from microcircuit.create_params import params_single, set_thalamic

if __name__ == "__main__":
    # simulation settings
    run_sim = True
    on_server = False
    run_analysis = True
    print_to_file = True

    # set ai segments
    n_seg_ai = 0
    start_ai = 2000.0
    seg_ai = 10000.0
    len_ai = seg_ai*n_seg_ai

    # set thalamic input
    n_stim = 20
    th_rate = 120.0 # Bruno, Simons, 2002: 1.4 spikes/20-ms deflection
    interval_stim = 2000.0
    ana_win = 40.0
    start_stim = start_ai + len_ai
    len_stim = interval_stim*n_stim
    stims = list(range(int(start_stim + interval_stim/2), int(start_stim + len_stim), int(interval_stim)))
    print('stims = {}'.format(stims))
    set_thalamic(stims, th_rate)

    # set others
    plot_half_len = 100.0
    if n_stim == 0:
        plot_center = start_ai
    else:
        plot_center = stims[-1]

    # check for: parameter scan or single-run
    try:
        pickle_path = sys.argv[1]    # path to pickle file
    except IndexError:  # single-run if no path input
        print('No argv[1]; single-run.')
        cwd = os.getcwd()

        # create pickle file
        pickle_path = os.path.join(cwd, 'para_dict.pickle')
        params_single(pickle_path)

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

    # cpu number / on server or not
    cpu_ratio = 0.5
    if mp.cpu_count() > 10:
        cpu_ratio = 1
        on_server = True

    # set print to file
    if print_to_file:
        exec(tools.set2txt(para_dict['sim_dict']['data_path']))

    # set simulation condition
    para_dict['sim_dict']['local_num_threads'] = \
        int(mp.cpu_count() * cpu_ratio)
    para_dict['sim_dict']['t_sim'] = start_ai + len_ai + len_stim

    # run simulation
    net = network.Network(para_dict['sim_dict'], para_dict['net_dict'],
                          para_dict['stim_dict'], para_dict['special_dict'])
    net.setup()
    if run_sim:
        net.simulate()

    # analysis
    if run_analysis:
        spikes = tools.Spikes(para_dict['sim_dict']['data_path'], 'spike_detector')
        mean_fr, std_fr = \
            tools.fire_rate(spikes, start_ai, start_ai + len_ai)
        if n_seg_ai > 0:
            t0 = time.time()
            tools.ai_score(spikes, start_ai, start_ai + len_ai, seg_len=seg_ai)
            print('ai analysis time = {}'.format(time.time() - t0))
        if n_stim > 0:
            t1 = time.time()
            tools.response(spikes, start_stim,
                           para_dict['stim_dict']['th_start'],
                           window=ana_win,
                           exportplot=True)
            print('response analysis time = {}'.format(time.time() - t1))
        tools.plot_raster(spikes, plot_center - plot_half_len, plot_center + plot_half_len)
        tools.fr_boxplot(para_dict['net_dict'], para_dict['sim_dict']['data_path'])

    # delete .gdf files to save space
    if on_server and os.path.isdir(para_dict['sim_dict']['data_path']):
        for item in os.listdir(para_dict['sim_dict']['data_path']):
            if item.endswith('.gdf'):
                os.remove(os.path.join(para_dict['sim_dict']['data_path'], item))

    if print_to_file:
        exec(tools.end2txt())
