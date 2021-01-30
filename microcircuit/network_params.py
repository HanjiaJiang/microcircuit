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

def get_psp_mtx(mtx_e, mtx_i, flg_e, flg_i, g=None, psp_mean=0.7543, psp_rel_sd=1.2704):
    mtx = np.zeros((16, 16))
    ipsp_pv2e = mtx_i[0, 1]
    if flg_e is False:
        if isinstance(g, int) or isinstance(g, float):
            mtx_e = np.full((4, 4), psp_mean)    # mean
        else:
            mtx_e = np.full((4, 4), psp_rel_sd)    # s.d. (relative)
    # assign by layer
    for a, row in enumerate(mtx_e):
        for b, epsp in enumerate(row):
            mtx[a*4:a*4+4, b*4:b*4+4] = epsp
            if isinstance(g, int) or isinstance(g, float):
                if flg_i:
                    mtx[a*4:a*4+4, b*4+1:b*4+4] = mtx_i[:, 1:]*np.abs(g/(ipsp_pv2e/epsp))
                else:
                    mtx[a*4:a*4+4, b*4+1:b*4+4] = epsp*g
    mtx = np.delete(mtx, [7, 11, 15], 0)
    mtx = np.delete(mtx, [7, 11, 15], 1)
    # print(repr(mtx))
    return mtx

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
        [1691, 88, 72, 90, 1656, 85, 48, 1095, 109, 105, 1288, 56, 67]), # mouse column (fairly distribution among 3 classes)
    # 'N_full': np.array([5099, 521, 67, 88, 4089, 287, 61, 3267, 540, 140, 4425, 290, 102]),
    # Connection probabilities. The first index corresponds to the targets
    # and the second to the sources.
    'conn_probs': # 200426
        np.array([ [0.0839, 0.3053, 0.4438, 0.0522, 0.1039, 0.0192, 0.0429, 0.0232, 0.0891, 0.1845, 0.    , 0.    , 0.    ],
                   [0.3621, 0.3323, 0.2061, 0.0806, 0.0042, 0.0155, 0.0298, 0.0254, 0.1918, 0.2215, 0.0001, 0.0001, 0.0054],
                   [0.2201, 0.4057, 0.0254, 0.2519, 0.0038, 0.0111, 0.0417, 0.0004, 0.022 , 0.0199, 0.    , 0.    , 0.    ],
                   [0.0262, 0.0574, 0.0662, 0.0318, 0.0024, 0.0097, 0.0162, 0.0002, 0.0009, 0.0056, 0.    , 0.0001, 0.005 ],

                   [0.0126, 0.0333, 0.0562, 0.0663, 0.1668, 0.4327, 0.261 , 0.0058, 0.0264, 0.0491, 0.    , 0.0021, 0.0232],
                   [0.0378, 0.0152, 0.0216, 0.0503, 0.0886, 0.3297, 0.3846, 0.009 , 0.0262, 0.0446, 0.0013, 0.0026, 0.0175],
                   [0.0379, 0.0172, 0.0227, 0.0458, 0.0859, 0.419 , 0.0264, 0.01  , 0.0303, 0.0441, 0.0017, 0.0021, 0.0217],

                   [0.0832, 0.0523, 0.0713, 0.0589, 0.0826, 0.0658, 0.0714, 0.091 , 0.178 , 0.1498, 0.0093, 0.0167, 0.0477],
                   [0.0698, 0.1199, 0.0439, 0.0186, 0.0362, 0.0288, 0.0216, 0.0804, 0.34  , 0.2468, 0.0047, 0.0122, 0.0243],
                   [0.0988, 0.0071, 0.0098, 0.0184, 0.0394, 0.0285, 0.024 , 0.0558, 0.1124, 0.0355, 0.0061, 0.014 , 0.0299],

                   [0.    , 0.0018, 0.0031, 0.0075, 0.0291, 0.0145, 0.0094, 0.0374, 0.0184, 0.0157, 0.0199, 0.3083, 0.286 ],
                   [0.0028, 0.0001, 0.0002, 0.002 , 0.0052, 0.0022, 0.0005, 0.0171, 0.    , 0.0032, 0.177 , 0.3355, 0.2817],
                   [0.0022, 0.    , 0.0002, 0.0011, 0.0047, 0.0019, 0.0003, 0.0161, 0.    , 0.0021, 0.1208, 0.3151, 0.0292]]),

    # Number of external connections to the different populations.
    # The order corresponds to the order in 'populations'.
    'K_ext': np.array([1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]),
    # Mean amplitude of excitatory postsynaptic potential (in mV).
    'PSP_e': 0.49,
    # Relative standard deviation of the postsynaptic potential.
    'PSP_sd': 0.27,
    # Relative inhibitory synaptic strength (in relative units).
    'g': -8,
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
    'ctsp': {
        'source': 'neske-lyr',
        'default': {
            'tau_m': 10.0, 'C_m': 200.0, 'E_L': -67.0, 'V_reset': -67.0, 'V_th': -40.0,
        },
        'neske-lyr':{
            'tau_m': [10.5, 3.1, 11.8, 10.9, 10.5, 3.1, 11.8, 12.1, 4.01, 11.7, 12.1, 4.01, 11.7],
            'C_m': [229.8, 93.9, 123.3, 86.5, 229.8, 93.9, 123.3, 269.2, 81.0, 146.8, 269.2, 81.0, 146.8],
            'E_L': [-67.4, -66.4, -59.9, -65.7, -67.4, -66.4, -59.9, -63.2, -67.1, -63.2, -63.2, -67.1, -63.2],
            'V_reset': [-67.4, -66.4, -59.9, -65.7, -67.4, -66.4, -59.9, -63.2, -67.1, -63.2, -63.2, -67.1, -63.2],
            'V_th': [-41.5, -41.6, -41.8, -43.7, -41.5, -41.6, -41.8, -45.2, -42.3, -48.1, -45.2, -42.3, -48.1],
        },
        'neske-avg':{
            'tau_m': [11.1, 3.6, 11.8, 10.9, 11.1, 3.6, 11.8, 11.1, 3.6, 11.8, 11.1, 3.6, 11.8],
            'C_m': [240.2, 85.5, 133.4, 86.5, 240.2, 85.5, 133.4, 240.2, 85.5, 133.4, 240.2, 85.5, 133.4],
            'E_L': [-65.9, -66.8, -61.5, -65.7, -65.9, -66.8, -61.5, -65.9, -66.8, -61.5, -65.9, -66.8, -61.5],
            'V_reset': [-65.9, -66.8, -61.5, -65.7, -65.9, -66.8, -61.5, -65.9, -66.8, -61.5, -65.9, -66.8, -61.5],
            'V_th': [-42.8, -42.0, -44.8, -43.7, -42.8, -42.0, -44.8, -42.8, -42.0, -44.8, -42.8, -42.0, -44.8],
        },
        'neske-old':{
            'tau_m': [13.0, 3.6, 11.8, 10.9, 13.0, 3.6, 11.8, 13.0, 3.6, 11.8, 13.0, 3.6, 11.8],
            'C_m': [322.0, 86.2, 134.0, 86.5, 322.0, 86.2, 134.0, 322.0, 86.2, 134.0, 322.0, 86.2, 134.0],
            'E_L': [-63.3, -66.8, -61.6, -65.7, -63.3, -66.8, -61.6, -63.3, -66.8, -61.6, -63.3, -66.8, -61.6],
            'V_reset': [-63.3, -66.8, -61.6, -65.7, -63.3, -66.8, -61.6, -63.3, -66.8, -61.6, -63.3, -66.8, -61.6],
            'V_th': [-45.6, -42.0, -45.0, -43.7, -45.6, -42.0, -45.0, -45.6, -42.0, -45.0, -45.6, -42.0, -45.0],
        }
    },
    'neuron_params': {
        # Membrane potential average for the neurons (in mV).
        'V0_mean': -62.0, #-58.0,
        # Standard deviation of the average membrane potential (in mV).
        'V0_sd': 5.0, #10.0,
        # Reset membrane potential of the neurons (in mV).
        'E_L': {'default': -67.0, 'Exc': -63.3, 'PV': -66.8, 'SOM': -61.6, 'VIP': -65.7},
        # Threshold potential of the neurons (in mV).
        'V_th': {'default': -40.0, 'Exc': -45.6, 'PV': -42.0, 'SOM': -45.0, 'VIP': -43.7},
        # Membrane potential after a spike (in mV).
        'V_reset': {'default': -67.0, 'Exc': -63.3, 'PV': -66.8, 'SOM': -61.6, 'VIP': -65.7},
        # Membrane capacitance (in pF).
        'C_m': {'default': 200.0, 'Exc': 322.0, 'PV': 86.2, 'SOM': 134.0, 'VIP': 86.5},
        # Membrane time constant (in ms).
        'tau_m': {'default': 10.0, 'Exc': 13.0, 'PV': 3.6, 'SOM': 11.8, 'VIP': 10.9},
        # Time constant of postsynaptic excitatory currents (in ms).
        'tau_syn_ex': 2.0, #2.1, # Allen Institue,
        # Time constant of postsynaptic inhibitory currents (in ms).
        'tau_syn_in': 4.0, #3.2, # Allen Institue,
        # Time constant of external postsynaptic excitatory current (in ms).
        'tau_syn_E': 0.5,   # not using
        # Refractory period of the neurons after a spike (in ms).
        't_ref': 2.0},
    # EPSPs: Lefort et al., 2009, Neuron
    'lyr-spe-epsp': {
        'use': False,
        'means':
            np.array([[0.7001, 0.7800, 0.4673, 0.3783],
                      [0.3433, 0.9500, 0.3767, 0.3783],
                      [0.7006, 0.6270, 0.6763, 0.2267],
                      [0.5813, 2.2700, 0.4025, 0.5300]]),
            # np.full((4, 4), 0.5), # previous
        'stds':
            np.array([[0.8958, 1.2372, 0.7581, 1.3884],
                      [0.5243, 1.3421, 0.9566, 1.3884],
                      [1.0520, 0.9618, 1.2379, 1.3884],
                      [0.8240, 1.3124, 0.8739, 1.3884]])
            # np.full((4, 4), 1.0) # previous
        },
    # IPSPs from different interneurons: Ma et al., 2012, J. Neurosci.
    'cell-spe-epsp': {
        'use': False,
        'means':
            np.array([[0.0000, -2.7000, -1.0400, -1.7460],
                      [0.0000, -1.7900, -1.3800, -1.7460],
                      [0.0000, -1.8200, -1.2100, -1.7460],
                      [0.0000, -2.1033, -1.2100, -1.7460]]),
        'stds':
            np.array([[0.0000, 0.2519, 0.1731, 0.1744],
                      [0.0000, 0.1061, 0.1594, 0.1744],
                      [0.0000, 0.1813, 0.1662, 0.1744],
                      [0.0000, 0.1798, 0.1662, 0.1744]])
        },
    # compensate w by U
    'U-compensate': True,
    # dc input
    'dc_extra': np.zeros(13).astype(float),
    # max firing rate
    'fmax': False,
    # STP
    'stp_dict': {},
    # selectivity
    'orient_tuning': False,
    'sel_inh_src': ['PV', 'SOM'],
    'sel_inh_trg': ['PV', 'SOM'],
    'k_th': 0.8,
    'k_e2e': 0.8,
    'k_e2i': 0.8,
    'k_i2e': 0.2,
    # recurrent weight distribution
    'recurrent_weight_distribution': 'lognormal',
    # in vivo EPSPs
    'in-vivo-epsp': {
        'thalamocortical': (0.49, 0.27),    # Bruno, Sakmann, 2006, Cortex is ...
        'intracortical': (0.43, 0.16),      # Jouhanneau et al., 2015, In vivo ...
    },
    }


def net_update(n_dict):
    updated_dict = {
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
        # epsp and ipsp
        'psp_means': get_psp_mtx(
            n_dict['lyr-spe-epsp']['means'],
            n_dict['cell-spe-epsp']['means'],
            n_dict['lyr-spe-epsp']['use'],
            n_dict['cell-spe-epsp']['use'],
            g=n_dict['g']),
            # psp_mean=.43,
            # psp_rel_sd=.16),
        'psp_stds': get_psp_mtx(
            n_dict['lyr-spe-epsp']['stds'],
            n_dict['cell-spe-epsp']['stds'],
            n_dict['lyr-spe-epsp']['use'],
            n_dict['cell-spe-epsp']['use']),
            # psp_mean=.43,
            # psp_rel_sd=.16),
    }
    n_dict.update(updated_dict)

# net_update(net_dict)
