import os
import sys
import pickle
import numpy as np
import multiprocessing as mp
import microcircuit.network as network
import microcircuit.analysis as analysis
import create_params as create

class RunNetwork:
    def __init__(self, run_sim=True, run_analysis=True, stp=2,
                do_ai=True, do_response=True, do_sel=False, indgs=None):
        self.run_sim = run_sim
        self.run_analysis = run_analysis
        self.do_ai = do_ai
        self.do_response = do_response
        self.do_selectivity = do_sel
        self.stp = stp
        if indgs is None:
            self.indgs = [750,1500,500,1250] if stp == 2 else [1000,1500,750,1000]
        else:
            self.indgs = indgs
        self.t_sim = 0.

    def set_wr(self, enable=False, testmode=True, seg_w=100.):
        self.wr = {
            'enabled': enable, 'testmode': testmode, 'seg_width': seg_w
        }

    def set_ai(self, n=1, start=2000., length=2000.):
        self.ai = {
            'start': start,
            'end': start + n*length,
            'n_seg': n,
            'len_seg': length
        }
        self.t_sim += self.ai['end']

    def set_th(self, n=0, rate=200., seg_w=1000., duration=10.,
                conn_probs=None, ana_win=40.,):
        len_stim = seg_w*n
        starts = list(range(int(self.t_sim + seg_w/2), \
            int(self.t_sim + len_stim), int(seg_w)))
        starts = np.array(starts).astype(float)
        self.thalamus_start = self.t_sim    # start of thalamus section
        self.thalamus_ana_win = ana_win
        self.thalamus = {
            'th_start': starts, # times of stimuli
            'th_rate': rate,
            'th_duration': duration,
            'conn_probs_th': conn_probs,
        }
        self.t_sim += seg_w*n

    def set_perturb(self, stype='poisson', n=0, duration=600., intrv=1000.,
                    targets=[0], levels=None):
        if levels is None:
            levels = np.arange(0., 401., 50.)
        self.perturbs = {
            'type': stype,
            'n_repeat': n,
            'start': self.t_sim,
            'duration': duration,
            'interval': intrv,
            'targets': targets,
            'levels': levels
        }
        self.t_sim += n*len(levels)*(duration+intrv)

    def set_dc_extra(self, targets=[], amps=[]):
        self.dc_extra = {'targets': targets, 'amplitudes': amps}

    def set_raster_plot(self, center=None, half=None):
        self.raster_half = 100.
        if isinstance(center, float) and isinstance(half, float):
            self.raster_center = center
            self.raster_half = half
        else:
            self.raster_center = run.ai['start']
            if len(run.thalamus['th_start']) > 0:
                self.raster_center = self.thalamus['th_start'][0]
            if run.perturbs['n_repeat'] > 0:
                self.raster_center = self.perturbs['start']

if __name__ == "__main__":
    # initiate with running and model settings
    run = RunNetwork(run_sim=True, run_analysis=True, stp=2,
                    do_ai=True, do_response=True, indgs=[750, 1500, 500, 1250])

    # ai segments
    run.set_ai(n=1, start=2000., length=5000.)

    # set thalamic input
    # Bruno, Simons, 2002: 1.4 spikes/20-ms deflection
    # Landisman, Connors, 2007, Cerebral Cortex: VPM >300 spikes/s in burst
    run.set_th(n=0, rate=100.)

    # perturb effect
    # perturb_levels = np.arange(0., 201., 100.).tolist()
    run.set_perturb(n=0, stype='poisson', duration=600., intrv=1000.,
                    targets=[2], levels=np.arange(0., 201., 100.))

    # weight recorder
    run.set_wr(enable=False, testmode=True, seg_w=100.)

    # dc_extra
    run.set_dc_extra(targets=[], amps=[])

    # others
    run.set_raster_plot(center=None, half=None)

    # initiate ScanParams
    scanparams = create.ScanParams(indgs=run.indgs, conn='0715', stp=run.stp, g=-8., bg=4.5)
    # scanparams.renew_conn(extrapolate=False)
    # scanparams.set_lognormal(False)
    if run.wr['enabled']:
        scanparams.net_dict['rec_dev'].append('weight_recorder')

    # get pickle, scan or single
    cwd = os.getcwd()
    single_path = None
    try:
        # scan, load pickle file
        arg1, arg2, arg3 = sys.argv[1], sys.argv[2], sys.argv[3]
        arg4, arg5, arg6 = sys.argv[4], sys.argv[5], sys.argv[6]
        pickle_path = arg1
        scanparams.load_pickle(pickle_path)
        lvls_str, lvls = scanparams.read_levels(pickle_path)
        # set constant parameters
        scanparams.set_indgs(run.indgs) # use the defined, if not scanned
        # set scanned parameters
        run.perturbs['targets'] = [int(arg3)]
        # scanparams.set_stp(sys.argv[3])
        scanparams.set_vip2som(arg4)
        scanparams.set_epsp(arg5)
        scanparams.set_ipsp(arg6)
        scanparams.set_g(lvls[0])
        scanparams.set_bg(lvls[1])
        scanparams.set_vip(lvls[2])
        # scanparams.set_exc(lvls[0])
        # scanparams.set_pv(lvls[1])
        # scanparams.set_som(lvls[2])
        scanparams.save_pickle(pickle_path)
    # single
    except IndexError:
        print('No scanning input; do single simulation')
        # handle data path and copy files
        # single_path = sys.argv[1]
        single_path = os.path.join(cwd, 'data/')
        os.system('mkdir -p ' + single_path)
        os.system('cp -r *.py *.yml microcircuit/ ' + single_path)
        # create pickle file
        pickle_path = os.path.join(single_path, 'para_dict.pickle')
        scanparams.save_pickle(pickle_path)

    # get parameters from pickle
    with open(pickle_path, 'rb') as handle:
        para_dict = pickle.load(handle)
    if single_path is not None:
        para_dict['sim_dict']['data_path'] = single_path
    data_path = para_dict['sim_dict']['data_path']

    # cpu number / on server or not
    on_server = (False if mp.cpu_count() <= 10 else True)
    cpu_ratio = (0.5 if mp.cpu_count() <= 10 else 1.)

    # set other parameters
    scanparams.set_thalamic(para_dict, run.thalamus)
    scanparams.set_perturb(para_dict, run.perturbs)
    for target, amp in zip(run.dc_extra['targets'], run.dc_extra['amplitudes']):
        para_dict['net_dict']['dc_extra'][target] = amp
    print('th_start = {}'.format(para_dict['stim_dict']['th_start']))

    # set simulation
    para_dict['sim_dict']['local_num_threads'] = int(mp.cpu_count() * cpu_ratio)
    para_dict['sim_dict']['t_sim'] = run.t_sim

    # delete existing files
    if os.path.isdir(data_path):
        os.system('rm {}'.format(os.path.join(data_path, '*.png')))
        if run.run_sim == True:
            os.system('rm {} {}'. \
                format(os.path.join(data_path, '*.gdf'), \
                os.path.join(data_path, '*.csv')))

    # initialize and run
    net = network.Network(para_dict['sim_dict'], para_dict['net_dict'],
                          para_dict['stim_dict'], test=run.wr['testmode'])
    net.setup()
    if run.run_sim:
        # print parameters
        scanparams.print_all(para_dict)
        net.simulate()

    # analysis
    if run.run_analysis:
        spikes = analysis.Spikes(data_path, para_dict['net_dict'])
        mean_fr, std_fr = \
            analysis.fire_rate(spikes, run.ai['start'], run.ai['end'])
        if run.do_ai and run.ai['n_seg'] > 0:
            analysis.gs_analysis(spikes, run.ai['start'], run.ai['end'], bw=10, seg_len=run.ai['len_seg'])
        if len(run.thalamus['th_start']) > 0:
            if run.do_response:
                analysis.response(spikes, run.thalamus_start, run.thalamus['th_start'],
                    window=run.thalamus_ana_win, bw=1., ptitle='no STPs', figsize=(8,5),
                    ylims=(-0.04, 0.35))
            if run.do_selectivity:
                analysis.selectivity(spikes, para_dict['stim_dict']['th_start'], duration=run.thalamus['duration'], raw=True)
                analysis.selectivity(spikes, para_dict['stim_dict']['th_start'], duration=run.thalamus['duration'], raw=False)
        analysis.perturb_calc(spikes, para_dict['stim_dict']['perturbs'])
        analysis.plot_raster(spikes, run.raster_center - run.raster_half, run.raster_center + run.raster_half)
        analysis.fr_plot(spikes)
        if run.wr['enabled']:
            spikes.stationary_musig(run.ai['start'], run.ai['end'], sw=run.wr['seg_width'], verify=run.wr['testmode'])
        spikes.verify_print(data_path)

    # delete recording files, copy .png files
    if on_server and os.path.isdir(data_path):
        os.system('rm {}'.format(os.path.join(data_path, '*.gdf')))
        if run.wr['testmode'] is True:
            os.system('rm -f {}'.format(os.path.join(data_path, '*.csv')))
        if run.perturbs['n_repeat'] > 0:
            affix = cwd.replace('/', '-') + '-' + data_path.replace('/', '-')
            os.chdir(data_path)
            os.system('for f in *.png; do cp -- \"$f\" \"../${{f%}}{}\"; done'.format('.' + affix + '.png'))
            os.chdir(cwd)
        # os.system('cp {} ./'.format(os.path.join(data_path, '*')))

    # copy exception
    xcpt_path = os.path.join(data_path, 'ai_xcpt.dat')
    xcpt_cp_path = xcpt_path.replace('/', '_')
    if os.path.isfile(xcpt_path):
        os.system('mkdir -p ../exception/')
        os.system('cp {} ../exception/{}'.format(xcpt_path, xcpt_cp_path))
