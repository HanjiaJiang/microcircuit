import os
import sys
import pickle
import multiprocessing as mp
import time
import microcircuit.network as network
import microcircuit.analysis as analysis
import microcircuit.create_params as create

if __name__ == "__main__":
    # simulation settings
    run_sim = True
    on_server = False
    run_analysis = True
    print_to_file = False

    #  settings
    do_ai = False
    do_response = False
    do_selectivity = False

    # set ai segments
    n_seg_ai = 1
    start_ai = 2000.0
    seg_ai = 2000.0
    len_ai = seg_ai*n_seg_ai

    # set thalamic input
    n_stim = 0
    # Bruno, Simons, 2002: 1.4 spikes/20-ms deflection
    # Landisman, Connors, 2007, Cerebral Cortex: VPM >300 spikes/s in burst
    th_rate = 120.0
    interval_stim = 1000.0
    ana_win = 40.0
    orient = False
    duration = 10.0
    start_stim = start_ai + len_ai
    len_stim = interval_stim*n_stim
    stims = list(range(int(start_stim + interval_stim/2), int(start_stim + len_stim), int(interval_stim)))

    # set others
    plot_half_len = 100.0
    if n_stim == 0:
        plot_center = start_ai
    else:
        plot_center = stims[-1]

    # initiate ScanParams
    scanparams = create.ScanParams()

    # check for: parameter scan or single-run
    try:
        pickle_path = sys.argv[1]    # path to pickle file
    except IndexError:  # single-run if no path input
        print('No scanning input; do single simulation')
        cwd = os.getcwd()

        # handle data path and copy files
        dpath = os.path.join(cwd, 'data')
        os.system('mkdir -p ' + dpath)
        os.system('cp run_network.py ' + dpath)
        os.system('cp config.yml ' + dpath)
        os.system('cp -r microcircuit/ ' + dpath)

        # create pickle file
        pickle_path = os.path.join(dpath, 'para_dict.pickle')
        scanparams.do_single(pickle_path)

    # assign parameters
    with open(pickle_path, 'rb') as handle:
        para_dict = pickle.load(handle)
    data_path = para_dict['sim_dict']['data_path']

    # cpu number / on server or not
    cpu_ratio = 0.5
    if mp.cpu_count() > 10:
        cpu_ratio = 1
        on_server = True

    # set print to file
    if print_to_file:
        exec(analysis.set2txt(data_path))

    # set simulation condition
    create.set_thalamic(para_dict, stims, th_rate, orient=orient, duration=duration)
    para_dict['sim_dict']['local_num_threads'] = int(mp.cpu_count() * cpu_ratio)
    para_dict['sim_dict']['t_sim'] = start_ai + len_ai + len_stim
    print('start_ai = {}, len_ai = {}'.format(start_ai, len_ai))
    print('thalamic_input = {}'.format(para_dict['stim_dict']['thalamic_input']))
    print('stims = {}'.format(para_dict['stim_dict']['th_start']))

    # initialize and run
    net = network.Network(para_dict['sim_dict'], para_dict['net_dict'],
                          para_dict['stim_dict'], para_dict['special_dict'])
    net.setup()
    if run_sim:
        # print parameters
        create.print_all(para_dict)
        net.simulate()

    # analysis
    if run_analysis:
        spikes = analysis.Spikes(data_path, 'spike_detector')
        mean_fr, std_fr = \
            analysis.fire_rate(spikes, start_ai, start_ai + len_ai)
        if do_ai and n_seg_ai > 0:
            t0 = time.time()
            analysis.ai_score(spikes, start_ai, start_ai + len_ai, bw=10.0, seg_len=seg_ai)
            print('ai_score() running time = {}'.format(time.time() - t0))
        if n_stim > 0:
            t1 = time.time()
            if do_response:
                analysis.response(spikes, start_stim, stims, window=ana_win)
            t2 = time.time()
            if do_selectivity:
                analysis.selectivity(spikes, para_dict['stim_dict']['th_start'], duration=duration, raw=True)
                analysis.selectivity(spikes, para_dict['stim_dict']['th_start'], duration=duration, raw=False)
            t3 = time.time()
            print('response() runing time = {}'.format(t2 - t1))
            print('selectivity() runing time = {}'.format(t3 - t2))
            # move hist-ltc.png to root folder
            # if on_server:
            #     os.system('mv {} {}'.format(os.path.join(data_path, 'hist-ltc.png'), 'hist-ltc_{}.png'.format(data_path.split('/')[-2])))

        if not on_server:
            analysis.plot_raster(spikes, plot_center - plot_half_len, plot_center + plot_half_len)
            analysis.fr_plot(spikes)

        spikes.verify_print(data_path)

    # delete .gdf files to save space
    if on_server and os.path.isdir(data_path):
        for item in os.listdir(data_path):
            if item.endswith('.gdf'):
                os.remove(os.path.join(data_path, item))

    if print_to_file:
        exec(analysis.end2txt())
