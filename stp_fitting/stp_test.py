import os
import sys
import copy
import nest
import pickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from microcircuit.functions import verify_collect, verify_print
matplotlib.rcParams['font.size'] = 20.0
np.set_printoptions(precision=2, linewidth=500, suppress=True)

class ConnTest:
    def __init__(self, syn_dict,
                    pre_subtype,
                    post_subtype,
                    pprs=None,
                    peaks=None,
                    spk_n=10,
                    spk_isi=50.0,
                    verify=False):
        self.syn_dict = syn_dict
        pre_subtype, post_subtype = self.set_subtype(pre_subtype, post_subtype)
        self.pre_subtype = pre_subtype
        self.post_subtype = post_subtype
        self.spk_n = spk_n
        self.spk_isi = spk_isi
        self.verify = verify
        self.bisyn_delay = 4.0
        self.set_params()
        self.set_labels()
        self.set_neurons()
        self.set_spkgen()
        self.set_mm()
        self.set_data(pprs, peaks)

    # change 'L23_Exc' to 'Exc'
    def set_subtype(self, pretype, posttype):
        subtypes = ['Exc', 'PV', 'SOM', 'VIP']
        for subtype in subtypes:
            if subtype in pretype:
                pretype = subtype
            if subtype in posttype:
                posttype = subtype
        return pretype, posttype

    def set_params(self):
        # neuron params
        self.neuron_model = 'iaf_psc_exp'
        self.subtypes = ['Exc', 'PV', 'SOM', 'VIP']
        self.ctsp = {
                # Reset membrane potential of the neurons (in mV).
                'E_L': {'default': 0.0, 'Exc': 0.0, 'PV': 0.0, 'SOM': 0.0, 'VIP': 0.0}, # Neske, Patrick, Connors, 2015 (in vitro)
                # 'E_L': {'default': -67.0, 'Exc': -63.3, 'PV': -66.8, 'SOM': -61.6, 'VIP': -65.7}, # Neske, Patrick, Connors, 2015 (in vitro)
                # Threshold potential of the neurons (in mV).
                'V_th': {'default': 40.0, 'Exc': 40.0, 'PV': 40.0, 'SOM': 40.0, 'VIP': 40.0}, # Gentet, Petersen, 2012 (in vivo)
                # 'V_th': {'default': -40.0, 'Exc': -41.0, 'PV': -40.5, 'SOM': -40.3, 'VIP': -41.2}, # Gentet, Petersen, 2012 (in vivo)
                # Membrane capacitance (in pF).
                'C_m': {'default': 200.0, 'Exc': 322.0, 'PV': 86.2, 'SOM': 134.0, 'VIP': 86.5}, # Neske, Patrick, Connors, 2015 (in vitro)
                # Membrane time constant (in ms).
                'tau_m': {'default': 10.0, 'Exc': 13.0, 'PV': 3.6, 'SOM': 11.8, 'VIP': 10.9}, # Neske, Patrick, Connors, 2015 (in vitro)
                # Time constant of postsynaptic excitatory currents (in ms).
                'tau_syn_ex': 2.0,
                # Time constant of postsynaptic inhibitory currents (in ms).
                'tau_syn_in': 4.0
            }

    def set_labels(self):
        self.color_labels = {
            'Exc': (68/255,119/255,170/255),
            'PV': (238/255,102/255,119/255),
            'SOM': (34/255,136/255,51/255),
            'VIP': (204/255,187/255,68/255)
        }

    def set_neurons(self, w=100.0):
        # set neurons
        self.pre_pop = nest.Create(self.neuron_model, self.spk_n)
        nest.SetStatus(self.pre_pop, {'E_L': 0.0, 'V_reset': 0.0, 'V_th': 5.0, 'V_m': 0.0})
        self.post_pop = nest.Create(self.neuron_model, self.spk_n)
        self.post_static = nest.Create(self.neuron_model)
        self.post_Uis1 = nest.Create(self.neuron_model)
        params_dict = {
                'E_L': self.ctsp['E_L'][self.post_subtype],
                'V_th': self.ctsp['V_th'][self.post_subtype],
                'C_m': self.ctsp['C_m'][self.post_subtype],
                'tau_m': self.ctsp['tau_m'][self.post_subtype],
                'tau_syn_ex': self.ctsp['tau_syn_ex'],
                'tau_syn_in': self.ctsp['tau_syn_in'],
                'V_m': self.ctsp['E_L'][self.post_subtype],
            }
        nest.SetStatus(self.post_pop, params_dict)
        nest.SetStatus(self.post_static, params_dict)
        nest.SetStatus(self.post_Uis1, params_dict)
        # syn_dict handling
        syn_dict = self.syn_dict
        if self.pre_subtype == 'Exc':
            tau_psc = params_dict['tau_syn_ex']
        else:
            w *= -1
            tau_psc = params_dict['tau_syn_in']
        syn_dict['tau_psc'] = tau_psc
        # static version
        syn_dict_stat = {'model': 'static_synapse'}
        syn_dict_stat['weight'] = w
        # U = 1 version
        syn_dict_Uis1 = copy.deepcopy(syn_dict)
        syn_dict_Uis1['weight'] = w
        syn_dict_Uis1['U'] = 1.0
        # tested version; adjust w for the release probability U
        syn_dict['weight'] = w
        if syn_dict['U'] != 0.0:
            syn_dict['weight'] /= syn_dict['U']
        # connect
        nest.Connect(self.pre_pop, self.post_pop, conn_spec={'rule': 'one_to_one'}, syn_spec=syn_dict)
        nest.Connect([self.pre_pop[-1]], self.post_static, syn_spec=syn_dict_stat)
        nest.Connect([self.pre_pop[-1]], self.post_Uis1,  syn_spec=syn_dict_Uis1)

    def set_spkgen(self, spk_w=1000.0):
        isi = self.spk_isi
        n = self.spk_n
        self.spks = nest.Create('spike_generator', n)
        block = isi*n
        self.spk_ts = []
        # spk_ts = np.array([])
        for i in range(n):
            self.spk_ts.append(np.arange(block, block+(i+1)*isi, isi))
            # spk_ts = np.concatenate((spk_ts, np.arange((2*i+1)*block, (2*i+1)*block+(i+1)*isi, isi)))
            nest.SetStatus([self.spks[i]], {'spike_times': self.spk_ts[i]})
        nest.Connect(self.spks, self.pre_pop, conn_spec={'rule': 'one_to_one'}, syn_spec={'weight': spk_w})

    def set_mm(self, resol=0.1):
        # for post-synaptic neurons
        self.mm = nest.Create('multimeter')
        nest.SetStatus(self.mm, {"withtime": True, "record_from": ["V_m"], 'interval': resol})
        nest.Connect(self.mm, self.pre_pop)
        nest.Connect(self.mm, self.post_pop)
        #
        self.mm_extra = nest.Create('multimeter')
        nest.SetStatus(self.mm_extra, {"withtime": True, "record_from": ["V_m"], 'interval': resol})
        nest.Connect(self.mm_extra, self.post_static)
        nest.Connect(self.mm_extra, self.post_Uis1)


    def reshape_mm(self, Vms, ts, cell_n, resolution=0.1):
        print('len(Vms)={}'.format(len(Vms)))
        freq = int(1.0/resolution)
        Vms = np.reshape(Vms, (int(len(Vms) / (cell_n*freq)), cell_n, freq))
        ts = np.reshape(ts, (int(len(ts) / (cell_n*freq)), cell_n, freq))
        Vms_new = []
        ts_new = []
        for i in range(cell_n):
            tmp_ts = tmp_Vms = None
            for j, k in zip(ts, Vms):
                if tmp_ts is None and tmp_Vms is None:
                    tmp_ts = j[i]
                    tmp_Vms = k[i]
                else:
                    tmp_ts = np.concatenate((tmp_ts, j[i]))
                    tmp_Vms = np.concatenate((tmp_Vms, k[i]))
            ts_new.append(tmp_ts)
            Vms_new.append(tmp_Vms)
        ts_new = np.array(ts_new)
        Vms_new = np.array(Vms_new)
        return Vms_new, ts_new

    def set_data(self, exp_pprs, exp_peaks):
        self.exp_data = {
            'PPRs': exp_pprs,
            'peaks': exp_peaks,
        }
        if isinstance(exp_peaks, np.ndarray):
            self.exp_data['peaks_norm'] = exp_peaks/exp_peaks[0]
        else:
            self.exp_data['peaks_norm'] = np.nan

        self.result = {
            # PSPs subtracted from lagging waveform
            'PSPs': np.full(self.spk_n, np.nan),
            # Paired-pulse ratios
            'PPRs': np.full(self.spk_n, np.nan),
            # original peaks, calculated from baseline
            'peaks': np.full(self.spk_n, np.nan),
            'peaks_norm': np.full(self.spk_n, np.nan),
            # fitness
            'fitness': np.nan
        }

    def run_sim(self, time):
        nest.Simulate(time)

    def plot_analysis(self, vs, ts, vs_extra=None):
        # plot for viewing
        len_pre = len(self.pre_pop)
        len_post = len(self.post_pop)
        fig, axes = plt.subplots(2, 1, figsize=(16, 12), constrained_layout=True)
        # presynaptic
        axes[0].plot(ts[len_pre-1], vs[len_pre-1], color=self.color_labels[self.pre_subtype], label='presynaptic')
        # postsynaptic
        axes[1].plot(ts[len_pre+len_post-1], vs[len_pre+len_post-1], color=self.color_labels[self.post_subtype],label='postsynaptic')
        axes[1].plot(ts[0], vs_extra[0], color='black', label='postsynaptic_static')
        axes[1].plot(ts[0], vs_extra[1], color='grey', label='postsynaptic_U=1')
        axes[1].set_ylim(self.ctsp['E_L'][self.post_subtype]-2.0, self.ctsp['E_L'][self.post_subtype]+2.0)
        axes[0].legend()
        axes[1].legend()
        plt.savefig('stp_test:plot_analysis().png')
        plt.show()
        plt.close()

    def run_analysis(self, png_name):
        # data for principle STP evaluation
        dmm = nest.GetStatus(self.mm)[0]
        Vms, ts = self.reshape_mm(dmm['events']['V_m'], dmm['events']['times'], len(self.pre_pop)+len(self.post_pop))
        Vms_post = Vms[len(self.pre_pop):len(self.pre_pop)+len(self.post_pop)]
        ts_post = ts[len(self.pre_pop):len(self.pre_pop)+len(self.post_pop)]

        # data for static and U=1 synapse
        dmm_extra = nest.GetStatus(self.mm_extra)[0]
        Vms_extra = self.reshape_mm(dmm_extra['events']['V_m'], dmm_extra['events']['times'], 2)[0]

        # plot
        if self.verify:
            self.plot_analysis(Vms, ts, Vms_extra)

        # calculate results
        self.calc_psp(Vms_post, ts_post)
        self.calc_peaks(Vms_post[-1], ts_post[-1])
        self.calc_fitness(png_name)

        # print experimental data and results
        print('exp_data:')
        for key, value in self.exp_data.items():
            print('{} = {}'.format(key, value))
        print()
        print('result:')
        for key, value in self.result.items():
            print('{} = {}'.format(key, value))

    # peaks: (de)polarizations
    def calc_peaks(self, vs, ts):
        # get parameters
        spk_n = self.spk_n
        isi = self.spk_isi
        spk_ts = self.spk_ts[-1]
        # plot
        if self.verify:
            fig = plt.figure(figsize=(16, 12))
        # calculation
        baseline = vs[ts==spk_ts[0]]
        for i in range(spk_n):
            # data with spk_n stimulations (the last one)
            t_start = spk_ts[i]
            # data of ith pulse
            vs_ith = vs[(ts>=t_start+self.bisyn_delay)&(ts<t_start+isi+self.bisyn_delay)]
            # determine the peak by Exc/Inh
            if self.pre_subtype == 'Exc':
                v_peak = np.max(vs_ith)
            else:
                v_peak = np.min(vs_ith)
            # print('baseline,peak={},{}'.format(baseline, v_peak))
            peak = float(np.abs(v_peak - baseline))
            self.result['peaks'][i] = peak
            self.result['peaks_norm'][i] = self.result['peaks'][i]/self.result['peaks'][0]
            if self.pre_subtype == 'Exc':
                plt.text(t_start+self.bisyn_delay, peak, '{:.4f}'.format(peak))
            else:
                plt.text(t_start+self.bisyn_delay, -peak, '{:.4f}'.format(peak))

        if self.verify:
            plt.plot(ts, vs, color='b')
            if self.post_subtype == 'Exc':
                plt.scatter(spk_ts+self.bisyn_delay, self.result['peaks'], color='green', label='peak amplitudes')
            else:
                plt.scatter(spk_ts+self.bisyn_delay, self.result['peaks'], color='green', label='peak amplitudes')
            plt.savefig('stp_test:calc_peaks()')


    def calc_psp(self, vs_raw, ts_raw):
        # get parameters
        isi = self.spk_isi
        spk_ts = self.spk_ts

        # plot raw data
        if self.verify:
            fig = plt.figure(figsize=(16, 12))
            for i in range(self.spk_n):
                if i == 0:
                    plt.plot(ts_raw[i], vs_raw[i], color='grey', label='raw')
                else:
                    plt.plot(ts_raw[i], vs_raw[i], color='grey')

        # to numpy
        # print(np.array([vs_raw[0]]))
        # print(np.diff(np.array(vs_raw), axis=0))
        vs = np.concatenate((np.array([vs_raw[0]]), np.diff(vs_raw, axis=0)), axis=0)
        ts = np.array(ts_raw)
        if vs.shape != ts.shape:
            print('calc_psp(): error in data shape!')
            return

        # calculation
        for i in range(self.spk_n):
            vs_calc = vs[i][(ts[i]>=spk_ts[i][i])&(ts[i]<spk_ts[i][i]+isi+self.bisyn_delay)]
            ts_tmp = ts[i][(ts[i]>=spk_ts[i][i])&(ts[i]<spk_ts[i][i]+isi+self.bisyn_delay)]
            if self.pre_subtype == 'Exc':
                psp = np.max(vs_calc) - vs_calc[0] # vs_calc[0] is baseline
            else:
                psp = vs_calc[0] - np.min(vs_calc)
            if self.verify:
                if i == 0:
                    plt.plot(ts_tmp, vs_calc, color='r', label='PSPs')
                else:
                    plt.plot(ts_tmp, vs_calc, color='r')
                plt.text(self.spk_ts[-1][i]+self.bisyn_delay, psp, '{:.4f}'.format(psp))
            # save result
            self.result['PSPs'][i] = psp
            # paired-pulse ratios
            self.result['PPRs'][i] = psp/self.result['PSPs'][0]

        # verify plot
        if self.verify:
            plt.scatter(self.spk_ts[-1]+self.bisyn_delay, self.result['PSPs'], color='green', label='PSP amplitudes')
            plt.legend()
            plt.savefig('stp_test:calc_psp()')
            plt.show()
            plt.close()


    def calc_fitness(self, png_name):
        exp_pprs = self.exp_data['PPRs']
        exp_peaks_norm = self.exp_data['peaks_norm']
        SumSqErr = 0.0
        fig, axs = plt.subplots(1, 1, figsize=(8, 6), constrained_layout=True)
        plt.xlabel('pulse')
        # when using PPRs
        if isinstance(exp_pprs, np.ndarray) and exp_pprs.ndim == 1:
            cnt = 0
            for i, ppr in enumerate(exp_pprs):
                # print('ppr={}'.format(ppr))
                if i==0 or np.isnan(ppr) or ppr <= 0:
                    continue
                SumSqErr += (self.result['PPRs'][i] - ppr)**2
                cnt += 1
                # print('fitness={:.2f}'.format(fitness))
            self.result['fitness'] = np.sqrt(SumSqErr/cnt)  # RMSE
            plt.ylabel('paired-pulse ratio')
            axs.plot(exp_pprs, marker='.', linestyle='solid', color='b', label='exp.')
            axs.plot(self.result['PPRs'][:len(exp_pprs)], marker='.', linestyle='solid', color='r', label='sim.')

        # when using peaks (depolarization/polarization)
        elif isinstance(exp_peaks_norm, np.ndarray) and exp_peaks_norm.ndim == 1:
            cnt = 0
            for i, peak in enumerate(exp_peaks_norm):
                if np.isnan(peak):
                    continue
                SumSqErr += (self.result['peaks_norm'][i] - peak)**2
                cnt += 1
            self.result['fitness'] = np.sqrt(SumSqErr/cnt)  # RMSE
            plt.ylabel('normalized (de)polarization')
            axs.plot(exp_peaks_norm, marker='.', linestyle='solid', color='b', label='exp.')
            axs.plot(self.result['peaks_norm'][:len(exp_peaks_norm)], marker='.', linestyle='solid', color='r', label='sim.')
        plt.ylim(bottom=0.0)
        plt.legend()
        plt.savefig(png_name)

if __name__ == '__main__':
    # neuron parameters
    pre_subtype = 'Exc'
    post_subtype = 'PV'

    # stimulation parameters
    spk_n = 10
    spk_isi = 10.0

    # experimental data
    pprs = None
    peaks = np.array([0.88, 0.83, 0.68, 0.57, 0.54, 0.46, 0.41, 0.38, 0.34, 0.34])

    # tsodyks Parameters
    U = 0.5
    F = 0.0
    D = 50.0

    # scanning input
    try:
        pickle_path = sys.argv[1]
        single_verify = False
        with open(pickle_path, 'rb') as h:
            para_dict = pickle.load(h)
            pre_subtype = para_dict['pre_subtype']
            post_subtype = para_dict['post_subtype']
            spk_n = para_dict['spk_n']
            spk_isi = para_dict['spk_isi']
            pprs = para_dict['pprs']
            peaks = para_dict['peaks']
            U = para_dict['U']
            F = para_dict['F']
            D = para_dict['D']
    except IndexError:
        single_verify = True
        pickle_path = 'test:calc_fitness()'
        para_dict = [pre_subtype, post_subtype, spk_n, spk_isi, pprs, peaks, U, F, D]
        print('No scanning input. Run single simulation.')

    verify_collect('{}'.format(para_dict), pickle_path)

    # STP dictionary
    if D == 0.0:
        D = 0.01
    syn_dict = {
        'model': 'tsodyks_synapse',
        'U': U,
        'tau_fac': F,
        'tau_rec': D
    }

    # initiate and run
    conntest = ConnTest(syn_dict,
                        pre_subtype,
                        post_subtype,
                        pprs=pprs,
                        peaks=peaks,
                        spk_n=spk_n,
                        spk_isi=spk_isi,
                        verify=single_verify)
    conntest.run_sim(spk_n*spk_isi*3)
    conntest.run_analysis('stp_' + pickle_path.replace('.pickle', '.png'))
    verify_print()

    # write results to file
    try:
        dat_path = sys.argv[2]
    except IndexError:
        dat_path = 'stp_test.dat'
    with open(dat_path, 'w') as f:
        f.write(pre_subtype + '\n')
        f.write(post_subtype + '\n')
        f.write(str(spk_n) + '\n')
        f.write(str(spk_isi) + '\n')
        f.write('{}\n'.format(pprs))
        f.write('{}\n'.format(peaks))
        f.write('{:.2f}\n'.format(U))
        f.write(str(F) + '\n')
        f.write(str(D) + '\n')
        f.write(str(conntest.result['fitness']))
        f.close()
    if single_verify:
        with open(dat_path.replace('.dat', '.txt'),'w') as f:
            f.write('exp_data:\n')
            for key, value in conntest.exp_data.items():
                f.write('{} = {}\n'.format(key, value))
            f.write('\n')
            f.write('result:\n')
            for key, value in conntest.result.items():
                f.write('{} = {}\n'.format(key, value))
            f.close()
