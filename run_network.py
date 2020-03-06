import os
import sys
import shutil
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
    run_ai = True
    run_response = True

    # timing, in ms
    plot_half_len = 100.0
    analysis_start = 1000.0
    analysis_segment = 1000.0
    n_segment = 1
    analysis_total_length = analysis_segment*n_segment
    analysis_interval = [analysis_start, analysis_start + analysis_total_length]

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
    handle.close()

    # cpu number / cluster or not
    cpu_ratio = 0.5
    if mp.cpu_count() > 10:
        cpu_ratio = 1
        on_server = True
    # else:
    #     exec(tools.set2txt(para_dict['sim_dict']['data_path']))
    #     pass
    exec(tools.set2txt(para_dict['sim_dict']['data_path']))
    para_dict['sim_dict']['local_num_threads'] = \
        int(mp.cpu_count() * cpu_ratio)

    para_dict['sim_dict']['t_sim'] = analysis_start + analysis_total_length

    # run simulation
    net = network.Network(para_dict['sim_dict'], para_dict['net_dict'],
                          para_dict['stim_dict'], para_dict['special_dict'])
    net.setup()
    if run_sim:
        net.simulate()

    # analysis
    if run_analysis:
        tools.plot_raster(
            para_dict['sim_dict']['data_path'], 'spike_detector',
            para_dict['stim_dict']['th_start'][0] - plot_half_len,
            para_dict['stim_dict']['th_start'][0] + plot_half_len)
        mean_fr, std_fr = \
            tools.fire_rate(para_dict['sim_dict']['data_path'], 'spike_detector',
                            analysis_interval[0], analysis_interval[1])
        tools.fr_boxplot(para_dict['net_dict'], para_dict['sim_dict']['data_path'])
        if run_ai:
            t0 = time.time()
            tools.ai_score(para_dict['sim_dict']['data_path'], 'spike_detector',
                           analysis_interval[0], analysis_interval[1], seg_len=analysis_segment)
            print('ai analysis time = {}'.format(time.time() - t0))
        if run_response:
            tools.response(para_dict['sim_dict']['data_path'], 'spike_detector', para_dict['stim_dict']['th_start'][0], window=20.0, n_stim=n_segment)

    # delete .gdf files to save space
    if on_server and os.path.isdir(para_dict['sim_dict']['data_path']):
        for item in os.listdir(para_dict['sim_dict']['data_path']):
            if item.endswith('.gdf'):
                os.remove(os.path.join(para_dict['sim_dict']['data_path'], item))

    exec(tools.end2txt())
    # if not on_server:
    #     exec(tools.end2txt())
    #     pass
