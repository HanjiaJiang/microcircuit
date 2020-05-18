import sys
import copy
import nest
import pickle
import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(precision=2, linewidth=500, suppress=True)

class ConnTest:
    def __init__(self, syn_dict,
                    pre_subtype,
                    post_subtype,
                    pprs=None,
                    peaks=None,
                    spk_n=10,
                    spk_isi=50.0):
        self.syn_dict = syn_dict
        pre_subtype, post_subtype = self.set_subtype(pre_subtype, post_subtype)
        self.pre_subtype = pre_subtype
        self.post_subtype = post_subtype
        self.spk_n = spk_n
        self.spk_isi = spk_isi
        self.set_params()
        self.set_labels()
        self.set_neurons()
        self.set_spkgen()
        self.set_mm()
        self.set_data(pprs, peaks)

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
                'E_L': {'default': -67.0, 'Exc': -63.3, 'PV': -66.8, 'SOM': -61.6, 'VIP': -65.7}, # Neske, Patrick, Connors, 2015 (in vitro)
                # Threshold potential of the neurons (in mV).
                'V_th': {'default': -40.0, 'Exc': -41.0, 'PV': -40.5, 'SOM': -40.3, 'VIP': -41.2}, # Gentet, Petersen, 2012 (in vivo)
                # Membrane capacitance (in pF).
                'C_m': {'default': 200.0, 'Exc': 322.0, 'PV': 86.2, 'SOM': 134.0, 'VIP': 86.5}, # Neske, Patrick, Connors, 2015 (in vitro)
                # Membrane time constant (in ms).
                'tau_m': {'default': 10.0, 'Exc': 13.0, 'PV': 3.6, 'SOM': 11.8, 'VIP': 10.9}, # Neske, Patrick, Connors, 2015 (in vitro)
                # Time constant of postsynaptic excitatory currents (in ms).
                'tau_syn_ex': 2.0, #2.1, # Allen Institue,
                # Time constant of postsynaptic inhibitory currents (in ms).
                'tau_syn_in': 4.0, #3.2, # Allen Institue,
            }

    def set_labels(self):
        self.color_labels = {
            'Exc': (68/255,119/255,170/255),
            'PV': (238/255,102/255,119/255),
            'SOM': (34/255,136/255,51/255),
            'VIP': (204/255,187/255,68/255)
        }

    def set_neurons(self, w=100.0):
        self.pre_neuron = nest.Create(self.neuron_model)
        self.post_neuron = nest.Create(self.neuron_model)
        params_dict = {
                'E_L': self.ctsp['E_L'][self.post_subtype],
                'V_th': self.ctsp['V_th'][self.post_subtype],
                'C_m': self.ctsp['C_m'][self.post_subtype],
                'tau_m': self.ctsp['tau_m'][self.post_subtype],
                'tau_syn_ex': self.ctsp['tau_syn_ex'],
                'tau_syn_in': self.ctsp['tau_syn_in'],
                'V_m': self.ctsp['E_L'][self.post_subtype],
            }
        nest.SetStatus(self.post_neuron, params_dict)
        syn_dict = self.syn_dict
        if self.pre_subtype == 'Exc':
            syn_dict['weight'] = w
            syn_dict['tau_psc'] = params_dict['tau_syn_ex']
        else:
            syn_dict['weight'] = -w
            syn_dict['tau_psc'] = params_dict['tau_syn_in']
        nest.Connect(self.pre_neuron, self.post_neuron, syn_spec=syn_dict)

    def set_spkgen(self, spk_w=3000.0):
        self.spks = nest.Create('spike_generator')
        isi = self.spk_isi
        n = self.spk_n
        block = isi*n
        self.spk_ts = []
        spk_ts = np.array([])
        for i in range(n):
            self.spk_ts.append(np.arange((2*i+1)*block, (2*i+1)*block+(i+1)*isi, isi))
            spk_ts = np.concatenate((spk_ts, np.arange((2*i+1)*block, (2*i+1)*block+(i+1)*isi, isi)))
        # print(self.spk_ts)
        nest.SetStatus(self.spks, {'spike_times': spk_ts})
        nest.Connect(self.spks, self.pre_neuron, syn_spec={'weight': spk_w})

    def set_mm(self, resol=0.1):
        self.mm = nest.Create('multimeter')
        nest.SetStatus(self.mm, {"withtime": True, "record_from": ["V_m"], 'interval': resol})
        nest.Connect(self.mm, self.pre_neuron)
        nest.Connect(self.mm, self.post_neuron)

    def reshape_mm(self, Vms, ts, cell_n, resolution=0.1):
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
            self.exp_data['peaks_norm'] = (exp_peaks - np.mean(exp_peaks))/np.std(exp_peaks)
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

    def run_analysis(self):
        # get data
        dmm = nest.GetStatus(self.mm)[0]
        Vms, ts = self.reshape_mm(dmm['events']['V_m'], dmm['events']['times'], 2)

        # plot for viewing
        # fig, axes = plt.subplots(2, 1, figsize=(16, 16))
        # axes[0].plot(ts[0], Vms[0], color=self.color_labels[self.pre_subtype])
        # axes[1].plot(ts[1], Vms[1], color=self.color_labels[self.post_subtype])
        # plt.show()
        # plt.savefig('stp_test.png')

        # calculate results
        self.calc_psp(Vms[1], ts[1])
        self.calc_fitness()

        # print experimental data and results
        print('exp_data:')
        for key, value in self.exp_data.items():
            print('{} = {}'.format(key, value))
        print()
        print('result:')
        for key, value in self.result.items():
            print('{} = {}'.format(key, value))

    def calc_psp(self, Vms, ts):
        # transpose data for calculation
        data = np.array([ts, Vms]).T

        # get parameters
        spk_n = self.spk_n
        isi = self.spk_isi
        spk_ts = self.spk_ts

        # verify
        # plt.plot(data[:, 0], data[:, 1], color='grey')

        # peaks = depolarizations/polarizations
        for i in range(spk_n):
            t_start = spk_ts[-1][i]
            baseline = data[ts==t_start][0, 1]
            # data of ith pulse of the last round
            data_ith = data[(ts>=t_start)&(ts<t_start+isi)]
            peak = np.max(data_ith[:, 1]) - baseline
            self.result['peaks'][i] = peak

        # normalized peaks
        peak_mean = np.mean(self.result['peaks'])
        peak_std = np.std(self.result['peaks'])
        for i in range(spk_n):
            self.result['peaks_norm'][i] = (self.result['peaks'][i] - peak_mean)/peak_std

        # verify
        # plt.scatter(spk_ts[-1], self.result['peaks'], color='g')

        # PSPs: PSP(i) is obtained by subtracting the 1st to (i-1)th pulses
        vms_base = np.zeros(len(Vms))
        for i in range(spk_n):
            t_start = spk_ts[i][0]
            # data of ith round (one round is 1~10 pulses)
            data_ith = data[(ts>=t_start) & (ts<t_start+spk_n*isi)]
            # verify
            # plt.plot(data_ith[:, 0], data_ith[:, 1], color='b')
            # copy data to keep the original
            data_calc = copy.deepcopy(data_ith)
            # subtract baseline (1st to (i-1)th pulses)
            data_calc[:, 1] -= vms_base[:len(data_calc)]
            # truncate the data of ith pulse
            data_calc = data_calc[(data_calc[:, 0]>=spk_ts[i][i]) & (data_calc[:, 0]<spk_ts[i][i]+isi)]
            # verify
            # plt.plot(data_calc[:, 0], data_calc[:, 1], color='r')
            # psp = amplitude of the obtained pulse data
            psp = np.max(data_calc[:, 1]) - np.min(data_calc[:, 1])
            # save pulse data of this round for subtraction in the next round
            vms_base = data_ith[:, 1]
            # save result
            self.result['PSPs'][i] = psp
            # paired-pulse ratios
            self.result['PPRs'][i] = psp/self.result['PSPs'][0]

        # verify
        # plt.show()


    def calc_fitness(self):
        exp_pprs = self.exp_data['PPRs']
        exp_peaks_norm = self.exp_data['peaks_norm']
        fitness = 0.0
        # when using PPRs
        if isinstance(exp_pprs, np.ndarray) and exp_pprs.ndim == 1:
            for i, ppr in enumerate(exp_pprs):
                # print('ppr={}'.format(ppr))
                if np.isnan(ppr):
                    continue
                fitness += (self.result['PPRs'][i] - ppr)**2
                # print('fitness={:.2f}'.format(fitness))
            self.result['fitness'] = fitness
        # when using peaks (depolarization/polarization)
        elif isinstance(exp_peaks_norm, np.ndarray) and exp_peaks_norm.ndim == 1:
            for i, peak in enumerate(exp_peaks_norm):
                if np.isnan(peak):
                    continue
                fitness += (self.result['peaks_norm'][i] - peak)**2
            self.result['fitness'] = fitness

if __name__ == '__main__':
    # simulation settings
    on_server = True

    # neuron parameters
    pre_subtype = 'Exc'
    post_subtype = 'PV'

    # stimulation parameters
    spk_n = 10
    spk_isi = 50.0

    # experimental data
    pprs = np.full(5, np.nan)
    pprs[0] = 1.0
    pprs[4] = 1.0
    peaks = np.arange(0.7, 0.2, -0.05)

    # tsodyks Parameters
    U = 0.75
    tau_fac = 0.0
    tau_rec = 100.0

    # scanning input
    try:
        pickle_path = sys.argv[1]
        with open(pickle_path, 'rb') as h:
            para_dict = pickle.load(h)
            pre_subtype = para_dict['pre_subtype']
            post_subtype = para_dict['post_subtype']
            spk_n = para_dict['spk_n']
            spk_isi = para_dict['spk_isi']
            pprs = para_dict['pprs']
            peaks = para_dict['peaks']
            U = para_dict['U']
            tau_fac = para_dict['tau_fac']
            tau_rec = para_dict['tau_rec']
    except IndexError:
        print('No scanning input. Run single simulation.')
        pass

    # STP dictionary
    syn_dict = {
        'model': 'tsodyks_synapse',
        'U': U,
        'tau_fac': tau_fac,
        'tau_rec': tau_rec,
    }

    # initiate and run
    conntest = ConnTest(syn_dict,
                        pre_subtype,
                        post_subtype,
                        pprs=pprs,
                        peaks=peaks,
                        spk_n=spk_n,
                        spk_isi=spk_isi)
    conntest.run_sim(spk_n*spk_isi*(spk_n*2+1))
    conntest.run_analysis()

    # write results to file
    try:
        dat_path = sys.argv[2]
    except IndexError:
        dat_path = 'test.dat'
    with open(dat_path, 'w') as f:
        f.write(str(conntest.result['fitness']))
        f.close()
    if not on_server:
        with open(dat_path.replace('.dat', '.txt'),'w') as f:
            f.write('exp_data:\n')
            for key, value in conntest.exp_data.items():
                f.write('{} = {}\n'.format(key, value))
            f.write('\n')
            f.write('result:\n')
            for key, value in conntest.result.items():
                f.write('{} = {}\n'.format(key, value))
            f.close()
