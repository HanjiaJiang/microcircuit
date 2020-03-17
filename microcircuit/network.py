import os
from microcircuit.functions import *
from microcircuit.network_params import net_update
np.set_printoptions(precision=4, suppress=True, linewidth=200)

class Network:
    """ Handles the setup of the network parameters and
    provides functions to connect the network and devices.

    Arguments
    ---------
    sim_dict
        dictionary containing all parameters specific to the simulation
        such as the directory the data is stored in and the seeds
        (see: sim_params.py)
    net_dict
         dictionary containing all parameters specific to the neurons
         and the network (see: network_params.py)

    Keyword Arguments
    -----------------
    stim_dict
        dictionary containing all parameter specific to the stimulus
        (see: stimulus_params.py)

    """
    def __init__(self, sim_dict, net_dict, stim_dict=None, spe_dict=None):
        net_update(net_dict, net_dict['g'])
        self.sim_dict = sim_dict
        self.net_dict = net_dict
        if stim_dict is not None:
            self.stim_dict = stim_dict
        else:
            self.stim_dict = None
        self.spe_dict = spe_dict
        self.stp_dict = spe_dict['stp_dict']
        self.data_path = sim_dict['data_path']
        if nest.Rank() == 0:
            if os.path.isdir(self.sim_dict['data_path']):
                print('data directory already exists')
            else:
                os.mkdir(self.sim_dict['data_path'])
                print('data directory created')
            print('Data will be written to %s' % self.data_path)

    def setup_nest(self):
        """ Hands parameters to the NEST-kernel.

        Resets the NEST-kernel and passes parameters to it.
        The number of seeds for the NEST-kernel is computed, based on the
        total number of MPI processes and threads of each.

        """
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
        """ Creates the neuronal populations.

        The neuronal populations are created and the parameters are assigned
        to them. The initial membrane potential of the neurons is drawn from a
        normal distribution. Scaling of the number of neurons and of the
        synapses is performed. If scaling is performed extra DC input is added
        to the neuronal populations.

        """
        self.nr_neurons = self.net_dict['N_full']
        self.synapses = get_total_number_of_synapses(self.net_dict)
        self.K_ext = self.net_dict['K_ext']
        self.w_ext = get_weight(self.net_dict['PSP_e'], self.net_dict)
        self.weight_mat = get_weights(self.net_dict)
        self.weight_mat_std = get_weight_stds(self.net_dict)
        # self.weight_mat = get_weight(
        #     self.net_dict['PSP_mean_matrix'], self.net_dict
        #     )
        # self.weight_mat_std = self.net_dict['PSP_std_matrix']   # ratio of std to mean!
        # print(self.weight_mat)
        # print(self.weight_mat_std)
        if self.net_dict['poisson_input']:
            self.DC_amp_e = np.zeros(len(self.net_dict['populations']))
        else:
            if nest.Rank() == 0:
                print(
                    '''
                    no poisson input provided
                    calculating dc input to compensate
                    '''
                    )
            self.DC_amp_e = compute_DC(self.net_dict, self.w_ext)


        # Create cortical populations.
        self.pops = []
        # HJ 190718
        if os.path.isfile(os.path.join(self.data_path, 'population_GIDs.dat')):
            os.remove(os.path.join(self.data_path, 'population_GIDs.dat'))
        pop_file = open(
            os.path.join(self.data_path, 'population_GIDs.dat'), 'w+'
            )
        for i, pop in enumerate(self.net_dict['populations']):
            population = nest.Create(
                self.net_dict['neuron_model'], int(self.nr_neurons[i])
                )
            E_L, V_th, C_m, tau_m = ctsp_assign(pop, self.net_dict, self.spe_dict)
            nest.SetStatus(
                population, {
                    'tau_syn_ex': self.net_dict['neuron_params']['tau_syn_ex'],
                    'tau_syn_in': self.net_dict['neuron_params']['tau_syn_in'],
                    'E_L': E_L,
                    'V_th': V_th,
                    'C_m': C_m,
                    'tau_m': tau_m,
                    'V_reset':  self.net_dict['neuron_params']['V_reset'],
                    't_ref': self.net_dict['neuron_params']['t_ref'],
                    'I_e': self.DC_amp_e[i]
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
        """ Creates the recording devices.

        Only devices which are given in net_dict['rec_dev'] are created.

        """
        self.spike_detector = []
        self.voltmeter = []
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

        if 'spike_detector' in self.net_dict['rec_dev']:
            if nest.Rank() == 0:
                print('Spike detectors created')
        if 'voltmeter' in self.net_dict['rec_dev']:
            if nest.Rank() == 0:
                print('Voltmeters created')

    def create_thalamic_input(self):
        """ This function creates the thalamic neuronal population if this
        is specified in stimulus_params.py.

        """
        if self.stim_dict['thalamic_input']:
            if nest.Rank() == 0:
                print('Thalamic input provided')
            self.thalamic_population = nest.Create(
                'parrot_neuron', self.stim_dict['n_thal']
                )
            self.thalamic_weight = get_weight(
                self.stim_dict['PSP_th'], self.net_dict
            )
            self.stop_th = (
                self.stim_dict['th_start'] + self.stim_dict['th_duration']
                )
            try:
                self.poisson_th = [nest.Create('poisson_generator', len(self.thalamic_population)) for x in range(len(self.stim_dict['th_start']))]
                set_thalamus_input(self.thalamic_population, self.poisson_th,
                                   self.stim_dict['th_start'], self.stop_th,
                                   self.stim_dict['th_rate'],
                                   self.stim_dict['orientation'],
                                   self.spe_dict)
            except NameError:
                print('\'set_thalamus_input()\' does not exist')
                self.poisson_th = nest.Create('poisson_generator')
                nest.SetStatus(
                    self.poisson_th, {
                        'rate': self.stim_dict['th_rate'],
                        'start': self.stim_dict['th_start'][0],
                        'stop': self.stop_th[0]
                    }
                )
                nest.Connect(self.poisson_th, self.thalamic_population)
            else:
                pass
            #
            self.nr_synapses_th = synapses_th_matrix(
                self.net_dict, self.stim_dict
            )
        else:
            if nest.Rank() == 0:
                print('Thalamic input not provided')

    def create_poisson(self):
        """ Creates the Poisson generators.

        If Poissonian input is provided, the Poissonian generators are created
        and the parameters needed are passed to the Poissonian generator.

        """
        if self.net_dict['poisson_input']:
            if nest.Rank() == 0:
                print('Poisson background input created')
            rate_ext = self.net_dict['bg_rate'] * self.K_ext
            self.poisson = []
            for i, target_pop in enumerate(self.pops):
                poisson = nest.Create('poisson_generator')
                nest.SetStatus(poisson, {'rate': rate_ext[i]})
                self.poisson.append(poisson)

    def create_dc_generator(self):
        """ Creates a DC input generator.

        If DC input is provided, the DC generators are created and the
        necessary parameters are passed to them.

        """
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

    def create_connections(self):
        """ Creates the recurrent connections.

        The recurrent connections between the neuronal populations are created.

        """
        if nest.Rank() == 0:
            print('Recurrent connections are established')
        mean_delays = self.net_dict['mean_delay_matrix']
        std_delays = self.net_dict['std_delay_matrix']
        for i, target_pop in enumerate(self.pops):
            for j, source_pop in enumerate(self.pops):
                synapse_nr = int(self.synapses[i][j])
                target_name = self.net_dict['populations'][i]
                source_name = self.net_dict['populations'][j]
                if synapse_nr >= 0.:
                    weight = self.weight_mat[i][j]
                    try:
                        weight = inh_weight(source_name, weight, self.spe_dict)
                    except NameError:
                        print('\'ins_weight()\' does not exist')
                    else:
                        pass

                    w_sd = abs(weight * self.weight_mat_std[i][j])
                    var_ln = np.log((weight ** 2 + w_sd ** 2) / weight ** 2)
                    mu_ln = np.log(abs(weight)) - var_ln / 2.
                    conn_dict_rec = {
                        'rule': 'fixed_total_number', 'N': synapse_nr
                        }
                    weight_dict = {
                       'distribution': 'lognormal', 'mu': mu_ln,
                       'sigma': np.sign(weight)*np.sqrt(var_ln)
                       }
                    delay_dict = {
                        'distribution': 'normal_clipped',
                        'mu': mean_delays[i][j], 'sigma': std_delays[i][j],
                        'low': self.sim_resolution
                        }

                    try:
                        syn_dict = assign_stp(source_name, target_name, weight_dict, delay_dict, self.stp_dict)
                    except NameError:
                        print('\'assign_stp()\' does not exist')
                        syn_dict = {
                            'model': 'static_synapse',
                            'weight': weight_dict,
                            'delay': delay_dict
                        }
                    else:
                        pass
                    try:
                        nr = connect_by_cluster(source_name, target_name, synapse_nr, syn_dict,
                                           source_pop, target_pop, self.spe_dict, conn_prob=self.net_dict['conn_probs'][i, j])
                    except NameError:
                        print('\'connect_by_cluster()\' does not exist')
                        nest.Connect(
                            source_pop, target_pop,
                            conn_spec=conn_dict_rec,
                            syn_spec=syn_dict
                            )
                    else:
                        pass
        # print('synapse numbers:')
        # print(syn_nr_bernoulli)

    def connect_poisson(self):
        """ Connects the Poisson generators to the microcircuit."""
        if nest.Rank() == 0:
            print('Poisson background input is connected')
        for i, target_pop in enumerate(self.pops):
            conn_dict_poisson = {'rule': 'all_to_all'}
            syn_dict_poisson = {
                'model': 'static_synapse',
                 'weight': self.w_ext,
                'delay': self.net_dict['poisson_delay']
                }
            nest.Connect(
                self.poisson[i], target_pop,
                conn_spec=conn_dict_poisson,
                syn_spec=syn_dict_poisson
                )

    def connect_thalamus(self):
        cell_types = ['Exc', 'PV', 'SOM', 'VIP', 'Exc', 'PV', 'SOM', 'Exc', 'PV', 'SOM', 'Exc', 'PV', 'SOM']
        psp = self.stim_dict['PSP_th']
        """ Connects the thalamic population to the microcircuit."""
        if nest.Rank() == 0:
            print('Thalamus connection established')
        for i, target_pop in enumerate(self.pops):
            C_m = self.net_dict['neuron_params']['C_m']['default']
            tau_m = self.net_dict['neuron_params']['tau_m']['default']
            tau_syn = self.net_dict['neuron_params']['tau_syn_ex']
            if isinstance(psp, np.ndarray):
                mu = calc_psc(psp[i], C_m, tau_m, tau_syn)
            else:
                mu = calc_psc(psp, C_m, tau_m, tau_syn)
            # print('psp={}'.format(psp))
            # print('mu={}'.format(mu))
            conn_dict_th = {
                'rule': 'fixed_total_number',
                'N': int(self.nr_synapses_th[i])
                }
            syn_dict_th = {
                'weight': {
                    'distribution': 'normal_clipped',
                    'mu': mu,
                    'sigma': (
                        mu * self.net_dict['PSP_sd']
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
            try:
                connect_thalamus_orientation(self.thalamic_population, target_pop, self.net_dict['populations'][i],
                                             self.nr_synapses_th[i], syn_dict_th, self.spe_dict, self.stim_dict['conn_probs_th'][i])
            except NameError:
                print('\'connect_thalamus_with_selectivity()\' does not exist')
                nest.Connect(
                    self.thalamic_population, target_pop,
                    conn_spec=conn_dict_th, syn_spec=syn_dict_th
                    )
            else:
                pass


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
        """ Execute subfunctions of the network.

        This function executes several subfunctions to create neuronal
        populations, devices and inputs, connects the populations with
        each other and with devices and input nodes.

        """
        self.setup_nest()
        self.create_populations()
        self.create_devices()
        self.create_thalamic_input()
        self.create_poisson()
        self.create_dc_generator()
        self.create_connections()
        if self.net_dict['poisson_input']:
            self.connect_poisson()
        if self.stim_dict['thalamic_input']:
            self.connect_thalamus()
        if self.stim_dict['dc_input']:
            self.connect_dc_generator()
        self.connect_devices()

    def simulate(self):
        """ Simulates the microcircuit."""
        nest.Simulate(self.sim_dict['t_sim'])
