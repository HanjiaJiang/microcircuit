import os
import sys
import pickle
import numpy as np
import multiprocessing as mp
import time
import microcircuit.network as network
import microcircuit.analysis as analysis
import microcircuit.create_params as create

if __name__ == "__main__":
    # simulation settings
    run_sim = True
    run_analysis = True
    print_to_file = False

    #  settings
    do_ai = True
    do_response = False
    do_selectivity = False
    weight_record = True
    test = True

    # set ai segments
    n_seg_ai, start_ai, seg_ai = 1, 0., 200.
    len_ai = seg_ai*n_seg_ai
    t_sim = start_ai + len_ai

    # set background input
    # indgs = [1000,1500,800,1000]
    indgs = [750,1500,500,250]

    # set thalamic input:
    # Bruno, Simons, 2002: 1.4 spikes/20-ms deflection
    # Landisman, Connors, 2007, Cerebral Cortex: VPM >300 spikes/s in burst
    n_stim, th_rate, stim_intrv = 0, 200., 1000.
    duration, ana_win, orient = 10., 40., False
    start_stim, len_stim = t_sim, stim_intrv*n_stim
    stims = list(range(int(start_stim + stim_intrv/2), int(start_stim + len_stim), int(stim_intrv)))
    t_sim += len_stim

    # set paradox effect input
    paradox_type = 'dc'
    n_paradox, paradox_start, = 0, t_sim
    paradox_duration, paradox_intrv = 600., 1000.
    paradox_pops = [1, 5, 8, 11]
    paradox_offsets = [0., 20., 40., 60., 80., 100., 120., 140., 160., 180.]
    # paradox_offsets = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90.]
    paradox_freq, paradox_ac_amp = 10., 0.1 # ac

    # set dc_extra
    dc_extra_targets, dc_extra_amps = [], []

    # set others
    plot_half_len = 100.0
    plot_center = 100.
    if len(stims) > 0:
        plot_center = stims[0]
    # plot_center = (paradox_start+n_paradox*(len(paradox_offsets)-1)*paradox_duration*2 if n_stim == 0 else stims[0])

    # initiate ScanParams
    scanparams = create.ScanParams(indgs)
    scanparams.vip2som(True)
    scanparams.set_g(8.)
    scanparams.set_bg(4.)
    scanparams.set_stp(2)
    if weight_record:
        scanparams.net_dict['rec_dev'].append('weight_recorder')
    # scanparams.net_dict['K_ext'] = np.array([750, 1500, 500, 1000,
    #                                          750, 1500, 500,
    #                                          1500, 1500, 0,
    #                                          1500, 1500, 0])

    # get pickle, scan or single
    cwd = os.getcwd()
    try:
        # load pickle file
        pickle_path = sys.argv[1]
        scanparams.load_pickle(pickle_path)
        lvls_str, lvls = scanparams.read_levels(pickle_path)
        # set parameters
        scanparams.set_stp(sys.argv[3])
        # to be improved
        # if int(sys.argv[3]) == 0:
        #     scanparams.set_indgs([1000,1500,750,1000])
        # elif int(sys.argv[3]) == 2:
        #     scanparams.set_indgs([750,1500,500,250])
        scanparams.vip2som(sys.argv[4])
        scanparams.set_epsp(sys.argv[5])
        scanparams.set_ipsp(sys.argv[6])
        scanparams.set_g(lvls[0])
        scanparams.set_bg(lvls[1])
        # scanparams.set_som(lvls[2])
        scanparams.save_pickle(pickle_path)
    except IndexError:  # single-run if no path input
        print('No scanning input; do single simulation')

        # handle data path and copy files
        dpath = os.path.join(cwd, 'data')
        os.system('mkdir -p ' + dpath)
        os.system('cp run_network.py ' + dpath)
        os.system('cp config.yml ' + dpath)
        os.system('cp -r microcircuit/ ' + dpath)

        # create pickle file
        pickle_path = os.path.join(dpath, 'para_dict.pickle')
        scanparams.save_pickle(pickle_path)

    # get parameters from pickle
    with open(pickle_path, 'rb') as handle:
        para_dict = pickle.load(handle)
    data_path = para_dict['sim_dict']['data_path']

    # cpu number / on server or not
    on_server = (False if mp.cpu_count() <= 10 else True)
    cpu_ratio = (0.5 if mp.cpu_count() <= 10 else 1.)

    # set print to file
    if print_to_file:
        exec(analysis.set2txt(data_path))

    # set simulation condition
    create.set_thalamic(para_dict, stims, th_rate, orient=orient, duration=duration)
    t_sim += create.set_paradox(para_dict, paradox_type, n_paradox, paradox_pops,
        paradox_offsets, paradox_start, paradox_duration, paradox_intrv,
        paradox_ac_amp, paradox_freq)
    para_dict['sim_dict']['local_num_threads'] = int(mp.cpu_count() * cpu_ratio)
    para_dict['sim_dict']['t_sim'] = t_sim
    for target, amp in zip(dc_extra_targets, dc_extra_amps):
        para_dict['net_dict']['dc_extra'][target] = amp
    print('stims = {}'.format(para_dict['stim_dict']['th_start']))

    # delete existing files
    if run_sim == True and os.path.isdir(data_path):
        os.system('rm {}/*.gdf {}/*.csv'.format(data_path, data_path))

    # initialize and run
    net = network.Network(para_dict['sim_dict'], para_dict['net_dict'],
                          para_dict['stim_dict'], test=test)
    net.setup()
    if run_sim:
        # print parameters
        create.print_all(para_dict)
        net.simulate()

    # analysis
    if run_analysis:
        spikes = analysis.Spikes(data_path, para_dict['net_dict']['rec_dev'])
        mean_fr, std_fr = \
            analysis.fire_rate(spikes, start_ai, start_ai + len_ai)
        if do_ai and n_seg_ai > 0:
            t0 = time.time()
            analysis.gs_analysis(spikes, start_ai, start_ai + len_ai, bw=10, seg_len=seg_ai)
            print('gs_analysis() running time = {}'.format(time.time() - t0))
        if n_stim > 0:
            t1 = time.time()
            if do_response:
                analysis.response(spikes, start_stim, stims, window=ana_win, interpol=True, bw=1.)
            t2 = time.time()
            if do_selectivity:
                analysis.selectivity(spikes, para_dict['stim_dict']['th_start'], duration=duration, raw=True)
                analysis.selectivity(spikes, para_dict['stim_dict']['th_start'], duration=duration, raw=False)
            t3 = time.time()
            print('response() runing time = {}'.format(t2 - t1))
            print('selectivity() runing time = {}'.format(t3 - t2))
        analysis.paradox_calc(spikes, para_dict['stim_dict']['paradox'])

        # if not on_server:
        analysis.plot_raster(spikes, plot_center - plot_half_len, plot_center + plot_half_len)
        analysis.fr_plot(spikes)
        if weight_record:
            spikes.plot_weight()
            spikes.compare_musig(start_ai, start_ai + len_ai)

        spikes.verify_print(data_path)

    # delete .gdf files to save space
    if on_server and os.path.isdir(data_path):
        os.system('rm {}/*.gdf {}/*.csv'.format(data_path, data_path))
        if n_paradox > 0:
            affix = cwd.replace('/', '-') + '-' + data_path.replace('/', '-')
            os.chdir(data_path)
            os.system('for f in *.png; do mv -- \"$f\" \"${{f%}}{}\"; done'.format('.' + affix + '.png'))
            os.chdir(cwd)

    if print_to_file:
        exec(analysis.end2txt())

    # copy exception
    xcpt_path = os.path.join(data_path, 'ai_xcpt.dat')
    xcpt_cp_path = xcpt_path.replace('/', '_')
    # print(xcpt_cp_path)
    if os.path.isfile(xcpt_path):
        os.system('mkdir -p ../exception/')
        os.system('cp {} ../exception/{}'.format(xcpt_path, xcpt_cp_path))
