import numpy as np


def get_mean_delays(mean_delay_exc, mean_delay_inh, number_of_pop):
    """ Creates matrix containing the delay of all connections.

    Arguments
    ---------
    mean_delay_exc
        Delay of the excitatory connections.
    mean_delay_inh
        Delay of the inhibitory connections.
    number_of_pop
        Number of populations.

    Returns
    -------
    mean_delays
        Matrix specifying the mean delay of all connections.

    """

    dim = number_of_pop
    mean_delays = np.zeros((dim, dim))
    mean_delays[:] = mean_delay_inh
    mean_delays[:, [0, 4, 7, 10]] = mean_delay_exc
    return mean_delays


def get_std_delays(std_delay_exc, std_delay_inh, number_of_pop):
    """ Creates matrix containing the standard deviations of all delays.

    Arguments
    ---------
    std_delay_exc
        Standard deviation of excitatory delays.
    std_delay_inh
        Standard deviation of inhibitory delays.
    number_of_pop
        Number of populations in the microcircuit.

    Returns
    -------
    std_delays
        Matrix specifying the standard deviation of all delays.

    """

    dim = number_of_pop
    std_delays = np.zeros((dim, dim))
    std_delays[:] = std_delay_inh
    std_delays[:, [0, 4, 7, 10]] = std_delay_exc
    return std_delays


def get_mean_PSP_matrix(PSP_e, g, number_of_pop):
    """ Creates a matrix of the mean evoked postsynaptic potential.

    The function creates a matrix of the mean evoked postsynaptic
    potentials between the recurrent connections of the microcircuit.
    The weight of the connection from L4E to L23E is doubled.

    Arguments
    ---------
    PSP_e
        Mean evoked potential.
    g
        Relative strength of the inhibitory to excitatory connection.
    number_of_pop
        Number of populations in the microcircuit.

    Returns
    -------
    weights
        Matrix of the weights for the recurrent connections.

    """
    dim = number_of_pop
    weights = np.zeros((dim, dim))
    exc = PSP_e
    inh = PSP_e * g
    weights[:] = inh
    weights[:, [0,4,7,10]] = exc
    weights[0, 4] = exc * 2
    return weights


def get_std_PSP_matrix(PSP_rel, number_of_pop):
    """ Relative standard deviation matrix of postsynaptic potential created.

    The relative standard deviation matrix of the evoked postsynaptic potential
    for the recurrent connections of the microcircuit is created.

    Arguments
    ---------
    PSP_rel
        Relative standard deviation of the evoked postsynaptic potential.
    number_of_pop
        Number of populations in the microcircuit.

    Returns
    -------
    std_mat
        Matrix of the standard deviation of postsynaptic potentials.

    """
    dim = number_of_pop
    std_mat = np.zeros((dim, dim))
    std_mat[:, :] = PSP_rel
    return std_mat

net_dict = {
    # Neuron model.
    'neuron_model': 'iaf_psc_exp',
    # The default recording device is the spike_detector. If you also
    # want to record the membrane potentials of the neurons, add
    # 'voltmeter' to the list.
    'rec_dev': ['spike_detector'],
    # Names of the simulated populations.
    'populations': ['L23_Exc', 'L23_PV', 'L23_SOM', 'L23_VIP', 'L4_Exc', 'L4_PV', 'L4_SOM', 'L5_Exc', 'L5_PV', 'L5_SOM', 'L6_Exc', 'L6_PV', 'L6_SOM'],
    # Number of neurons in the different populations. The order of the
    # elements corresponds to the names of the variable 'populations'.
    # 190717: round up to fit 8 clusters
    'N_full': np.array(
        [1688, 136, 48, 48, 1656, 88, 48, 1096, 112, 104, 1288, 64, 56]), # mouse column (rounded)
        # [5096, 520, 64, 88, 4088, 288, 64, 3264, 544, 144, 4424, 288, 104]),  # rat column (rounded)
    # rat C2 barrel column totally 18976 cells
    # 'N_full': np.array([5099, 521, 67, 88, 4089, 287, 61, 3267, 540, 140, 4425, 290, 102]),
    # Connection probabilities. The first index corresponds to the targets
    # and the second to the sources.
    'conn_probs': # 190707
        np.array([
               [0.0872, 0.3173, 0.4612, 0.0448, 0.1056, 0.4011, 0.0374, 0.0234, 0.09  , 0.1864, 0.    , 0.    , 0.    ],
               [0.3763, 0.3453, 0.2142, 0.0683, 0.0802, 0.012 , 0.0257, 0.0257, 0.1937, 0.2237, 0.0001, 0.0001, 0.0062],
               [0.2288, 0.1342, 0.1242, 0.2618, 0.0033, 0.0097, 0.0363, 0.0003, 0.0222, 0.018 , 0.    , 0.    , 0.    ],
               [0.0224, 0.0516, 0.0567, 0.0274, 0.0021, 0.0086, 0.0142, 0.0002, 0.0008, 0.0051, 0.    , 0.0001, 0.0048],

               [0.0128, 0.0668, 0.049 , 0.0584, 0.1764, 0.4577, 0.2761, 0.0059, 0.0232, 0.0427, 0.    , 0.0017, 0.0212],
               [0.0317, 0.0121, 0.0198, 0.0428, 0.0937, 0.3487, 0.4068, 0.0072, 0.0231, 0.0369, 0.0009, 0.002 , 0.0157],
               [0.033 , 0.0144, 0.0198, 0.2618, 0.2906, 0.4432, 0.0386, 0.0087, 0.0257, 0.0384, 0.001 , 0.0018, 0.0198],

               [0.0841, 0.0528, 0.072 , 0.0539, 0.0844, 0.0546, 0.0621, 0.0957, 0.1871, 0.1575, 0.0094, 0.0139, 0.0418],
               [0.0705, 0.1211, 0.0444, 0.0165, 0.0315, 0.0225, 0.0183, 0.0846, 0.3574, 0.2594, 0.0029, 0.0102, 0.0212],
               [0.0998, 0.0072, 0.0089, 0.2618, 0.0343, 0.0225, 0.0209, 0.0587, 0.1182, 0.0427, 0.0038, 0.0124, 0.0262],

               [0.    , 0.0018, 0.0028, 0.0068, 0.0297, 0.0125, 0.0084, 0.0381, 0.017 , 0.0128, 0.021 , 0.3249, 0.3014],
               [0.0025, 0.0001, 0.0003, 0.002 , 0.0045, 0.0016, 0.0004, 0.0149, 0.    , 0.0031, 0.1865, 0.3535, 0.2968],
               [0.0021, 0.    , 0.0002, 0.2618, 0.0004, 0.0014, 0.0003, 0.0141, 0.    , 0.0019, 0.1062, 0.3321, 0.0379]]),
    # Number of external connections to the different populations.
    # The order corresponds to the order in 'populations'.
    'K_ext': np.array([2000, 2000, 1500, 500, 2000, 2000, 1500, 2000, 2000, 1500, 2000, 2000, 1500]),
    # Mean amplitude of excitatory postsynaptic potential (in mV).
    'PSP_e': 0.5,
    # Relative standard deviation of the postsynaptic potential.
    'PSP_sd': 1.0,
    # Relative inhibitory synaptic strength (in relative units).
    'g': -4,
    # Rate of the Poissonian spike generator (in Hz).
    'bg_rate': 4.0,
    # Turn Poisson input on or off (True or False).
    'poisson_input': True,
    # Delay of the Poisson generator (in ms).
    'poisson_delay': 1.5,
    # Mean delay of excitatory connections (in ms).
    'mean_delay_exc': 1.2, #1.5,
    # Mean delay of inhibitory connections (in ms).
    'mean_delay_inh': 0.7, #0.75,
    # Relative standard deviation of the delay of excitatory and
    # inhibitory connections (in relative units).
    'rel_std_delay': 0.5,
    # Parameters of the neurons.
    # 190701 cell-type specific parameters
    'neuron_params': {
        # Membrane potential average for the neurons (in mV).
        'V0_mean': -62.0, #-58.0,
        # Standard deviation of the average membrane potential (in mV).
        'V0_sd': 5.0, #10.0,
        # Reset membrane potential of the neurons (in mV).
        'E_L': {'default': -67.0, 'Exc': -63.3, 'PV': -66.8, 'SOM': -61.6, 'VIP': -65.7}, #-67.0,
        # Threshold potential of the neurons (in mV).
        'V_th': {'default': -40.0, 'Exc': -41.0, 'PV': -40.5, 'SOM': -40.3, 'VIP': -41.2}, # Gentet, Petersen, 2012 (in vivo)
        # 'V_th': {'default': -40.0, 'Exc': -45.6, 'PV': -42.9, 'SOM': -45.0, 'VIP': -43.7}, # Neske, Patrick, Connors, 2015 (in vitro)
        # Membrane potential after a spike (in mV).
        'V_reset': -67.0, #-65.0,
        # Membrane capacitance (in pF).
        'C_m': {'default': 200.0, 'Exc': 322.0, 'PV': 86.2, 'SOM': 134.0, 'VIP': 86.5}, #200.0, #250.0,
        # Membrane time constant (in ms).
        'tau_m': {'default': 10.0, 'Exc': 13.0, 'PV': 3.6, 'SOM': 11.8, 'VIP': 10.9}, #7.0, #10.0,
        # Time constant of postsynaptic excitatory currents (in ms).
        'tau_syn_ex': 2.1, # 1.74, # 0.5,
        # Time constant of postsynaptic inhibitory currents (in ms).
        'tau_syn_in': 3.2, # 4.6, # 0.5,
        # Time constant of external postsynaptic excitatory current (in ms).
        'tau_syn_E': 0.5,   # not using
        # Refractory period of the neurons after a spike (in ms).
        't_ref': 2.0},
    'animal': 'rat',
    'renew_conn': False,
    'w_dict': {
        'psp_mtx':
            np.array([[0.70, 0.78, 0.47, 0.23],
                      [0.34, 0.95, 0.38, 0.23],
                      [0.70, 0.63, 0.68, 0.23],
                      [0.70, 2.27, 0.40, 0.53]]),
            # np.full((4, 4), 0.5), # previous
        'psp_std_mtx':
            np.array([[0.8958, 1.2372, 0.7228, 1.0000],
                      [0.4540, 1.3421, 1.0000, 1.0000],
                      [1.0520, 0.9618, 1.2379, 1.0000],
                      [1.0520, 1.3124, 0.8739, 1.3884]])
            # np.full((4, 4), 1.0) # previous
        }
    }


def net_update(n_dict, g):
    updated_dict = {
        # PSP mean matrix.
        'PSP_mean_matrix': get_mean_PSP_matrix(
            n_dict['PSP_e'], n_dict['g'], len(n_dict['populations'])
        ),
        # PSP std matrix.
        'PSP_std_matrix': get_std_PSP_matrix(
            n_dict['PSP_sd'], len(n_dict['populations'])
        ),
        # mean delay matrix.
        'mean_delay_matrix': get_mean_delays(
            n_dict['mean_delay_exc'], n_dict['mean_delay_inh'],
            len(n_dict['populations'])
        ),
        # std delay matrix.
        'std_delay_matrix': get_std_delays(
            n_dict['mean_delay_exc'] * n_dict['rel_std_delay'],
            n_dict['mean_delay_inh'] * n_dict['rel_std_delay'],
            len(n_dict['populations'])
        ),
    }
    n_dict.update(updated_dict)
