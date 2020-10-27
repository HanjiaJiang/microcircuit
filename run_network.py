import os
import sys
import pickle
import numpy as np
import multiprocessing as mp
import microcircuit.network as network
import microcircuit.analysis as analysis
import microcircuit.create_params as create

if __name__ == "__main__":
    # running settings
    run_sim = False
    run_analysis = True
    print_to_file = False

    #  model settings
    do_ai = True
    do_response = False
    do_selectivity = False
    do_weight, testmode_weight, weight_seg_width = False, True, 100.

    # ai segments
    n_seg_ai, start_ai, seg_ai = 1, 2000., 5000.
    len_ai = seg_ai*n_seg_ai
    t_sim = start_ai + len_ai

    # background input
    # indgs = [1000,1500,750,1000]
    indgs = [750,1500,500,1000]

    # thalamic input
    # Bruno, Simons, 2002: 1.4 spikes/20-ms deflection
    # Landisman, Connors, 2007, Cerebral Cortex: VPM >300 spikes/s in burst
    n_stim, th_rate, stim_intrv = 0, 120., 1000.
    duration, ana_win, orient = 10., 40., False
    start_stim, len_stim = t_sim, stim_intrv*n_stim
    stims = list(range(int(start_stim + stim_intrv/2), int(start_stim + len_stim), int(stim_intrv)))
    conn_probs_th = None
    # conn_probs_th = np.array([0., 0., 0., 0., 0.4, 0.4, 0., 0., 0., 0.0, 0., 0., 0.0])
    t_sim += len_stim

    # paradox effect
    paradox_type = 'dc'
    n_paradox, paradox_start, = 0, t_sim
    paradox_duration, paradox_intrv = 600., 1000.
    paradox_pops = [3] #[1, 5, 8, 11]
    paradox_offsets = [0., 20., 40., 60., 80., 100., 120., 140., 160., 180.]
    # paradox_offsets = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90.]
    paradox_freq, paradox_ac_amp = 10., 0.1 # ac

    # dc_extra of all populations
    dc_extra_targets, dc_extra_amps = [], []

    # others
    plot_center, plot_half_len = start_ai, 100.
    if len(stims) > 0:
        plot_center = stims[0]

    # initiate ScanParams
    scanparams = create.ScanParams(indgs)
    scanparams.set_g(8.)
    scanparams.set_bg(4.5)
    scanparams.set_stp(2)
    if do_weight:
        scanparams.net_dict['rec_dev'].append('weight_recorder')
    # scanparams.net_dict['mean_delay_exc'] = 0.2
    # scanparams.net_dict['mean_delay_inh'] = 0.2
    # scanparams.net_dict['poisson_delay'] = 0.2
    # scanparams.stim_dict['delay_th'] = np.full(13, 0.2)
    # scanparams.stim_dict['delay_th_sd'] = np.full(13, 0.1)
    # scanparams.net_dict['K_ext'] = np.array([750, 1500, 500, 1000,
    #                                          750, 1500, 500,
    #                                          1250, 1250, 0,
    #                                          750, 1500, 500])

    # get pickle, scan or single
    cwd = os.getcwd()
    try:
        # scan, load pickle file
        pickle_path = sys.argv[1]
        scanparams.load_pickle(pickle_path)
        lvls_str, lvls = scanparams.read_levels(pickle_path)
        # set constant parameters
        scanparams.set_indgs(indgs) # use the defined, if not scanned
        # set scanned parameters
        scanparams.set_stp(sys.argv[3])
        scanparams.set_vip(sys.argv[4])
        scanparams.set_epsp(sys.argv[5])
        scanparams.set_ipsp(sys.argv[6])
        scanparams.set_exc(lvls[0])
        scanparams.set_pv(lvls[1])
        scanparams.set_som(lvls[2])
        scanparams.save_pickle(pickle_path)
    # single
    except IndexError:
        print('No scanning input; do single simulation')
        # handle data path and copy files
        single_path = os.path.join(cwd, 'data')
        os.system('mkdir -p ' + single_path)
        os.system('cp run_network.py ' + single_path)
        os.system('cp config.yml ' + single_path)
        os.system('cp -r microcircuit/ ' + single_path)
        # create pickle file
        pickle_path = os.path.join(single_path, 'para_dict.pickle')
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

    # set other parameters
    create.set_thalamic(para_dict, stims, th_rate, orient=orient,
        duration=duration, conn_probs=conn_probs_th)
    t_sim += create.set_paradox(para_dict, paradox_type, n_paradox, paradox_pops,
        paradox_offsets, paradox_start, paradox_duration, paradox_intrv,
        paradox_ac_amp, paradox_freq)
    for target, amp in zip(dc_extra_targets, dc_extra_amps):
        para_dict['net_dict']['dc_extra'][target] = amp
    print('stims = {}'.format(para_dict['stim_dict']['th_start']))

    # set simulation
    para_dict['sim_dict']['local_num_threads'] = int(mp.cpu_count() * cpu_ratio)
    para_dict['sim_dict']['t_sim'] = t_sim

    # delete existing files
    if run_sim == True and os.path.isdir(data_path):
        os.system('rm {}/*.gdf {}/*.csv'.format(data_path, data_path))

    # initialize and run
    net = network.Network(para_dict['sim_dict'], para_dict['net_dict'],
                          para_dict['stim_dict'], test=testmode_weight)
    net.setup()
    if run_sim:
        # print parameters
        create.print_all(para_dict)
        net.simulate()

    # analysis
    if run_analysis:
        spikes = analysis.Spikes(data_path, para_dict['net_dict'])
        mean_fr, std_fr = \
            analysis.fire_rate(spikes, start_ai, start_ai + len_ai)
        if do_ai and n_seg_ai > 0:
            analysis.gs_analysis(spikes, start_ai, start_ai + len_ai, bw=10, seg_len=seg_ai)
        if n_stim > 0:
            if do_response:
                analysis.response(spikes, start_stim, stims, window=ana_win, interpol=True, bw=1.)
            if do_selectivity:
                analysis.selectivity(spikes, para_dict['stim_dict']['th_start'], duration=duration, raw=True)
                analysis.selectivity(spikes, para_dict['stim_dict']['th_start'], duration=duration, raw=False)
        analysis.paradox_calc(spikes, para_dict['stim_dict']['paradox'])
        analysis.plot_raster(spikes, plot_center - plot_half_len, plot_center + plot_half_len)
        analysis.fr_plot(spikes)
        if do_weight:
            spikes.stationary_musig(start_ai, start_ai + len_ai, sw=weight_seg_width, verify=testmode_weight)
        spikes.verify_print(data_path)

    # delete recording files, move .png files
    if on_server and os.path.isdir(data_path):
        os.system('rm {}/*.gdf'.format(data_path))
        if testmode_weight is True:
            os.system('rm {}/*.csv'.format(data_path))
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
    if os.path.isfile(xcpt_path):
        os.system('mkdir -p ../exception/')
        os.system('cp {} ../exception/{}'.format(xcpt_path, xcpt_cp_path))
