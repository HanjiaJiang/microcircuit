import os
from microcircuit.functions import *
from microcircuit.network_params import net_update
np.set_printoptions(precision=4, suppress=True, linewidth=100)

class Network:
    def __init__(self, sim_dict, net_dict, stim_dict=None, test=True):
        net_update(net_dict)
        self.sim_dict = sim_dict
        self.net_dict = net_dict
        if stim_dict is not None:
            self.stim_dict = stim_dict
        else:
            self.stim_dict = None
        self.data_path = sim_dict['data_path']
        if nest.Rank() == 0:
            if os.path.isdir(self.sim_dict['data_path']):
                print('data directory already exists')
            else:
                os.mkdir(self.sim_dict['data_path'])
                print('data directory created')
            print('Data will be written to %s' % self.data_path)
        self.test = test

    def setup_nest(self):
        nest.ResetKernel()
        master_seed = self.sim_dict['master_seed']
        if nest.Rank() == 0:
            print('Master seed: %i ' % master_seed)
        nest.SetKernelStatus(
            {'local_num_threads': self.sim_dict['local_num_threads']}
            )
        N_tp = nest.GetKernelStatus(['total_num_virtual_procs'])[0]
        if nest.Rank() == 0:
            print('Number of total processes: %i' % N_tp)
        rng_seeds = list(
            range(
                master_seed + 1 + N_tp,
                master_seed + 1 + (2 * N_tp)
                )
            )
        grng_seed = master_seed + N_tp
        if nest.Rank() == 0:
            print(
                'Seeds for random number generators of virtual processes: %r'
                % rng_seeds
                )
            print('Global random number generator seed: %i' % grng_seed)
        self.pyrngs = [np.random.RandomState(s) for s in list(range(
            master_seed, master_seed + N_tp))]
        self.sim_resolution = self.sim_dict['sim_resolution']
        kernel_dict = {
            'resolution': self.sim_resolution,
            'grng_seed': grng_seed,
            'rng_seeds': rng_seeds,
            'overwrite_files': self.sim_dict['overwrite_files'],
            'print_time': self.sim_dict['print_time'],
            }
        nest.SetKernelStatus(kernel_dict)

    def create_populations(self):
        self.nr_neurons = self.net_dict['N_full']
        self.synapses = get_total_number_of_synapses(self.net_dict) # not using
        self.K_ext = self.net_dict['K_ext']
        # self.w_ext = get_weight(self.net_dict['PSP_e'], self.net_dict)
        self.weight_mat = get_weight_mtx(self.net_dict)
        self.weight_mat_std = self.net_dict['psp_stds']
        self.dc_extra = self.net_dict['dc_extra']
        ctsp = self.net_dict['ctsp'][self.net_dict['ctsp']['source']]
        # Create cortical populations.
        self.pops = []
        if os.path.isfile(os.path.join(self.data_path, 'population_GIDs.dat')):
            os.remove(os.path.join(self.data_path, 'population_GIDs.dat'))
        pop_file = open(
            os.path.join(self.data_path, 'population_GIDs.dat'), 'w+'
            )
        for i, pop in enumerate(self.net_dict['populations']):
            population = nest.Create(
                self.net_dict['neuron_model'], int(self.nr_neurons[i])
                )
            E_L, V_reset, V_th = ctsp['E_L'][i], ctsp['V_reset'][i], ctsp['V_th'][i]
            C_m, tau_m = ctsp['C_m'][i], ctsp['tau_m'][i]
            # print('{} ctsp: {}, {}, {}, {}, {}'.format(pop, E_L, V_th, C_m, tau_m, V_reset))
            nest.SetStatus(
                population, {
                    'tau_syn_ex': self.net_dict['neuron_params']['tau_syn_ex'],
                    'tau_syn_in': self.net_dict['neuron_params']['tau_syn_in'],
                    'E_L': E_L,
                    'V_th': V_th,
                    'C_m': C_m,
                    'tau_m': tau_m,
                    'V_reset': V_reset,
                    't_ref': self.net_dict['neuron_params']['t_ref'],
                    'I_e': self.dc_extra[i]
                    # 'I_e': self.DC_amp_e[i]
                    }
                )
            self.pops.append(population)
            pop_file.write('%d  %d \n' % (population[0], population[-1]))
        pop_file.close()
        for thread in np.arange(nest.GetKernelStatus('local_num_threads')):
            # Using GetNodes is a work-around until NEST 3.0 is released. It
            # will issue a deprecation warning.
            local_nodes = nest.GetNodes(
                [0], {
                    'model': self.net_dict['neuron_model'],
                    'thread': thread
                    }, local_only=True
                )[0]
            vp = nest.GetStatus(local_nodes)[0]['vp']
            # vp is the same for all local nodes on the same thread
            nest.SetStatus(
                local_nodes, 'V_m', self.pyrngs[vp].normal(
                    self.net_dict['neuron_params']['V0_mean'],
                    self.net_dict['neuron_params']['V0_sd'],
                    len(local_nodes))
                    )

    def create_devices(self):
        self.spike_detector = []
        self.voltmeter = []
        self.weight_recorder = {}
        for i, pop in enumerate(self.pops):
            if 'spike_detector' in self.net_dict['rec_dev']:
                recdict = {
                    'withgid': True,
                    'withtime': True,
                    'to_memory': False,
                    'to_file': True,
                    'label': os.path.join(self.data_path, 'spike_detector')
                    }
                dummy = nest.Create('spike_detector', params=recdict)
                self.spike_detector.append(dummy)
            if 'voltmeter' in self.net_dict['rec_dev']:
                recdictmem = {
                    'interval': self.sim_dict['rec_V_int'],
                    'withgid': True,
                    'withtime': True,
                    'to_memory': False,
                    'to_file': True,
                    'label': os.path.join(self.data_path, 'voltmeter'),
                    'record_from': ['V_m'],
                    }
                volt = nest.Create('voltmeter', params=recdictmem)
                self.voltmeter.append(volt)
        if 'weight_recorder' in self.net_dict['rec_dev']:
            self.copysynapses = []
            wrdict = {
                'to_memory': False,
                'to_file': True,
                'label': os.path.join(self.data_path, 'weight_recorder'),
                }
            # one for recurrent, one for background input
            wr_recrr = nest.Create('weight_recorder', params=wrdict)
            wr_bg = nest.Create('weight_recorder', params=wrdict)
            self.weight_recorder['recurrent'] = wr_recrr
            self.weight_recorder['background'] = wr_bg

        if 'spike_detector' in self.net_dict['rec_dev']:
            if nest.Rank() == 0:
                print('Spike detectors created')
        if 'voltmeter' in self.net_dict['rec_dev']:
            if nest.Rank() == 0:
                print('Voltmeters created')

    def create_thalamic_input(self):
        if self.stim_dict['thalamic_input']:
            if nest.Rank() == 0:
                print('Thalamic input provided')
            self.thalamic_population = nest.Create(
                'parrot_neuron', self.stim_dict['n_thal']
                )
            self.stop_th = (
                self.stim_dict['th_start'] + self.stim_dict['th_duration']
                )
            # plural poisson generator for orientation (to be improved?)
            if self.net_dict['orient_tuning']:
                n_poisson = self.stim_dict['n_thal']
            else:
                n_poisson = 1
            self.poisson_th = [nest.Create('poisson_generator', n_poisson) for x in range(len(self.stim_dict['th_start']))]
            set_thalamus(self.thalamic_population, self.poisson_th,
                           self.stim_dict['th_start'], self.stop_th,
                           self.stim_dict['th_rate'],
                           self.stim_dict['orientation'],
                           self.net_dict)
            # self.poisson_th = nest.Create('poisson_generator')
            # nest.SetStatus(
            #     self.poisson_th, {
            #         'rate': self.stim_dict['th_rate'],
            #         'start': self.stim_dict['th_start'][0],
            #         'stop': self.stop_th[0]
            #     }
            # )
            # nest.Connect(self.poisson_th, self.thalamic_population)
            # only used in fixed_total_number connection
            self.nr_synapses_th = synapses_th_matrix(
                self.net_dict, self.stim_dict
            )
        else:
            if nest.Rank() == 0:
                print('Thalamic input not provided')

    def create_poisson(self):
        if self.net_dict['poisson_input']:
            if nest.Rank() == 0:
                print('Poisson background input created')
            # rate_ext = self.net_dict['bg_rate']
            rate_ext = self.net_dict['bg_rate'] * self.K_ext
            self.poisson = []
            for i, target_pop in enumerate(self.pops):
                poisson = nest.Create('poisson_generator')
                # nest.SetStatus(poisson, {'rate': rate_ext})
                nest.SetStatus(poisson, {'rate': rate_ext[i]})
                self.poisson.append(poisson)

    def create_dc_generator(self):
        if self.stim_dict['dc_input']:
            if nest.Rank() == 0:
                print('DC generator created')
            dc_amp_stim = self.net_dict['K_ext'] * self.stim_dict['dc_amp']
            self.dc = []
            if nest.Rank() == 0:
                print('DC_amp_stim', dc_amp_stim)
            for i, target_pop in enumerate(self.pops):
                dc = nest.Create(
                    'dc_generator', params={
                        'amplitude': dc_amp_stim[i],
                        'start': self.stim_dict['dc_start'],
                        'stop': (
                            self.stim_dict['dc_start'] +
                            self.stim_dict['dc_dur']
                            )
                        }
                    )
                self.dc.append(dc)

    def create_perturbs(self):
        if self.stim_dict['perturbs']['n_repeat'] > 0:
            self.perturbs = []
            duration = self.stim_dict['perturbs']['duration']
            perturbs_type = self.stim_dict['perturbs']['type']
            generator_type = '{}_generator'.format(perturbs_type)
            ctsp = self.net_dict['ctsp'][self.net_dict['ctsp']['source']]
            for i, level in enumerate(self.stim_dict['perturbs']['levels']):
                print('perturbs starts = {}'.format(self.stim_dict['perturbs']['starts'][str(level)]))
                if level <= 0:
                    continue
                for j, start in enumerate(self.stim_dict['perturbs']['starts'][str(level)]):
                    params = {'start': start, 'stop': start + duration}
                    if perturbs_type == 'dc':
                        params['amplitude'] = level
                    elif perturbs_type == 'ac':
                        params['amplitude'] = self.stim_dict['perturbs']['ac']['amplitude']
                        params['frequency'] = self.stim_dict['perturbs']['ac']['frequency']
                        params['phase'] = self.stim_dict['perturbs']['ac']['phase']
                    elif perturbs_type == 'poisson':
                        params['rate'] = level
                    print(params)
                    gen = nest.Create(generator_type, params=params)
                    self.perturbs.append(gen)
                    for k in self.stim_dict['perturbs']['targets']:
                        if perturbs_type == 'poisson':
                            w = calc_psc(self.net_dict['PSP_e'],
                                ctsp['C_m'][k], ctsp['tau_m'][k],
                                self.net_dict['neuron_params']['tau_syn_ex'])
                            nest.Connect(gen, self.pops[k], syn_spec={'weight': w})
                        else:
                            nest.Connect(gen, self.pops[k])


    def create_connections(self):
        if nest.Rank() == 0:
            print('Recurrent connections are established')
        # renew_conn(self.net_dict)
        verify_collect('weight_mat=\n{}\n'.format(self.weight_mat), 'lognormal')
        verify_collect('weight_mat_std=\n{}\n'.format(self.weight_mat_std), 'lognormal')
        # eq_inh_conn(self.net_dict['N_full'], self.net_dict['conn_probs'])
        for i, target_pop in enumerate(self.pops):
            for j, source_pop in enumerate(self.pops):
                p = self.net_dict['conn_probs'][i, j]
                if p <= 0.:
                    continue
                synapse_nr = int(self.synapses[i][j])
                target_name = self.net_dict['populations'][i]
                source_name = self.net_dict['populations'][j]
                w = self.weight_mat[i][j]
                w_sd = abs(w * self.weight_mat_std[i][j])
                delay = self.net_dict['mean_delay_matrix'][i][j]
                delay_sd = self.net_dict['std_delay_matrix'][i][j]
                syn_dict = assign_syn(source_name,
                                        target_name,
                                        w,
                                        w_sd,
                                        delay,
                                        delay_sd,
                                        self)
                connect_recurrent(source_name,
                                    target_name,
                                    synapse_nr,
                                    syn_dict,
                                    source_pop,
                                    target_pop,
                                    self.net_dict,
                                    bernoulli_prob=p)
                verify_collect('{} to {}: (w, w_sd, syn_dict) = {:.4f}, {:.4f}\n{}\n'.format(source_name, target_name, w, w_sd, syn_dict), 'lognormal')


    def connect_poisson(self):
        """ Connects the Poisson generators to the microcircuit."""
        self.bg_parrots = []
        copysynapse = 'static_synapse_wr_bg'
        ctsp = self.net_dict['ctsp'][self.net_dict['ctsp']['source']]
        if 'weight_recorder' in self.net_dict['rec_dev']:
            if copysynapse not in self.copysynapses:
                print('copysynapse = {}'.format(copysynapse))
                nest.CopyModel('static_synapse', copysynapse, \
                    {'weight_recorder': self.weight_recorder['background'][0]})
                self.copysynapses.append(copysynapse)
        if nest.Rank() == 0:
            print('Poisson background input is connected')
        for i, target_pop in enumerate(self.pops):
            conn_dict_poisson = {'rule': 'all_to_all'}
            w = calc_psc(self.net_dict['PSP_e'],
                ctsp['C_m'][i], ctsp['tau_m'][i],
                self.net_dict['neuron_params']['tau_syn_ex'])
            # syn_dict_poisson = {
            #     'model': 'static_synapse',
            #      'weight': {
            #          'distribution': 'normal_clipped',
            #          'mu': w,
            #          'sigma': (
            #              w * self.stim_dict['PSP_sd']
            #              ),
            #          'low': 0.0
            #          },
            #     'delay': self.net_dict['poisson_delay']
            #     }
            syn_dict_poisson = {
                'model': 'static_synapse',
                 'weight': w,
                'delay': self.net_dict['poisson_delay']
                }
            if 'weight_recorder' in self.net_dict['rec_dev']:
                if self.test==True and i != 0:
                    pass
                else:
                    syn_dict_poisson['model'] = copysynapse
            bg_parrot = nest.Create('parrot_neuron', self.net_dict['N_full'][i])
            nest.Connect(self.poisson[i], bg_parrot, conn_spec=conn_dict_poisson)
            nest.Connect(
                bg_parrot, target_pop,
                # self.poisson[i], target_pop,
                conn_spec={'rule': 'one_to_one'},
                # conn_spec=conn_dict_poisson,
                syn_spec=syn_dict_poisson
                )
            self.bg_parrots.append(bg_parrot)

    def connect_thalamus(self):
        """ Connects the thalamic population to the microcircuit."""
        psp = self.stim_dict['PSP_th']
        ctsp = self.net_dict['ctsp'][self.net_dict['ctsp']['source']]
        if nest.Rank() == 0:
            print('Thalamus connection established')
        for i, target_pop in enumerate(self.pops):
            if self.stim_dict['conn_probs_th'][i] == 0.0:
                continue
            C_m, tau_m = ctsp['C_m'][i], ctsp['tau_m'][i]
            tau_syn = self.net_dict['neuron_params']['tau_syn_ex']
            if isinstance(psp, np.ndarray):
                mu = calc_psc(psp[i], C_m, tau_m, tau_syn)
            else:
                mu = calc_psc(psp, C_m, tau_m, tau_syn)
            # conn_dict_th = {
            #     'rule': 'fixed_total_number',
            #     'N': int(self.nr_synapses_th[i])
            #     }
            # lognormal to be implemented?
            syn_dict_th = {
                'weight': {
                    'distribution': 'normal_clipped',
                    'mu': mu,
                    'sigma': (
                        mu * self.stim_dict['PSP_sd']
                        ),
                    'low': 0.0
                    },
                'delay': {
                    'distribution': 'normal_clipped',
                    'mu': self.stim_dict['delay_th'][i],
                    'sigma': self.stim_dict['delay_th_sd'][i],
                    'low': self.sim_resolution
                    }
                }
            connect_tc(
                self.thalamic_population,
                target_pop,
                self.net_dict['populations'][i],
                self.nr_synapses_th[i],
                syn_dict_th,
                self.net_dict,
                self.stim_dict['conn_probs_th'][i])
            # nest.Connect(
            #     self.thalamic_population, target_pop,
            #     conn_spec=conn_dict_th, syn_spec=syn_dict_th
            #     )

    def connect_dc_generator(self):
        """ Connects the DC generator to the microcircuit."""
        if nest.Rank() == 0:
            print('DC Generator connection established')
        for i, target_pop in enumerate(self.pops):
            if self.stim_dict['dc_input']:
                nest.Connect(self.dc[i], target_pop)

    def connect_devices(self):
        """ Connects the recording devices to the microcircuit."""
        if nest.Rank() == 0:
            if ('spike_detector' in self.net_dict['rec_dev'] and
                    'voltmeter' not in self.net_dict['rec_dev']):
                print('Spike detector connected')
            elif ('spike_detector' not in self.net_dict['rec_dev'] and
                    'voltmeter' in self.net_dict['rec_dev']):
                print('Voltmeter connected')
            elif ('spike_detector' in self.net_dict['rec_dev'] and
                    'voltmeter' in self.net_dict['rec_dev']):
                print('Spike detector and voltmeter connected')
            else:
                print('no recording devices connected')
        for i, target_pop in enumerate(self.pops):
            if 'voltmeter' in self.net_dict['rec_dev']:
                nest.Connect(self.voltmeter[i], target_pop)
            if 'spike_detector' in self.net_dict['rec_dev']:
                nest.Connect(target_pop, self.spike_detector[i])

    def setup(self):
        """ Execute subfunctions of the network."""
        self.setup_nest()
        self.create_populations()
        self.create_devices()
        # self.create_thalamic_input()
        self.create_poisson()
        self.create_dc_generator()
        self.create_perturbs()
        self.create_connections()
        if self.net_dict['poisson_input']:
            self.connect_poisson()
        if self.stim_dict['thalamic_input']:
            self.create_thalamic_input()
            self.connect_thalamus()
        if self.stim_dict['dc_input']:
            self.connect_dc_generator()
        self.connect_devices()

    def simulate(self):
        """ Simulates the microcircuit."""
        verify_print(self.data_path)
        nest.Simulate(self.sim_dict['t_sim'])
