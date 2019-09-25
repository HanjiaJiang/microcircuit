# -*- coding: utf-8 -*-
#
# network_params.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

'''
pynest microcircuit parameters
------------------------------

Network parameters for the microcircuit.

Hendrik Rothe, Hannah Bos, Sacha van Albada; May 2016
'''

import numpy as np

# HJ
# from scan_params import *

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
    #HJ
    mean_delays[:] = mean_delay_inh
    mean_delays[:, [0, 4, 7, 10]] = mean_delay_exc
    #mean_delays[:, 0:dim:2] = mean_delay_exc
    #mean_delays[:, 1:dim:2] = mean_delay_inh
    #HJ#
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
    #HJ
    std_delays[:] = std_delay_inh
    std_delays[:, [0, 4, 7, 10]] = std_delay_exc
    #std_delays[:, 0:dim:2] = std_delay_exc
    #std_delays[:, 1:dim:2] = std_delay_inh
    #HJ#
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
    #HJ
    weights[:] = inh
    weights[:, [0,4,7,10]] = exc
    weights[0, 4] = exc * 2
    #weights[:, 0:dim:2] = exc
    #weights[:, 1:dim:2] = inh
    #weights[0, 2] = exc * 2
    #HJ#
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
    'populations': ['L23_E', 'L23_PV', 'L23_SOM', 'L23_VIP', 'L4_E', 'L4_PV', 'L4_SOM', 'L5_E', 'L5_PV', 'L5_SOM', 'L6_E', 'L6_PV', 'L6_SOM'],
    # Number of neurons in the different populations. The order of the
    # elements corresponds to the names of the variable 'populations'.
    # 190717: round up to fit 8 clusters
    'N_full': np.array(
        [5096, 520, 64, 88, 4088, 288, 64, 3264, 544, 144, 4424, 288, 104]),
    # 190616 model of c2 barrel column totally 18976 cells
    # 'N_full': np.array([5099, 521, 67, 88, 4089, 287, 61, 3267, 540, 140, 4425, 290, 102]),
    # Mean rates of the different populations in the non-scaled version
    # of the microcircuit. Necessary for the scaling of the network.
    # The order corresponds to the order in 'populations'.
    'full_mean_rates':
        np.array([0.971, 2.868, 2.868, 2.868, 4.746, 5.396, 5.396, 8.142, 9.078, 9.078, 0.991, 7.523, 7.523]),
    # Connection probabilities. The first index corresponds to the targets
    # and the second to the sources.
    'conn_probs':
    # 190707
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
        # np.array(
        #     [[0.1009, 0.1689, 0.0437, 0.0818, 0.0323, 0.,     0.0076, 0.],
        #      [0.1346, 0.1371, 0.0316, 0.0515, 0.0755, 0.,     0.0042, 0.],
        #      [0.0077, 0.0059, 0.0497, 0.135,  0.0067, 0.0003, 0.0453, 0.],
        #      [0.0691, 0.0029, 0.0794, 0.1597, 0.0033, 0.,     0.1057, 0.],
        #      [0.1004, 0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0.],
        #      [0.0548, 0.0269, 0.0257, 0.0022, 0.06,   0.3158, 0.0086, 0.],
        #      [0.0156, 0.0066, 0.0211, 0.0166, 0.0572, 0.0197, 0.0396, 0.2252],
        #      [0.0364, 0.001,  0.0034, 0.0005, 0.0277, 0.008,  0.0658, 0.1443]]
        #     ),
    # Number of external connections to the different populations.
    # The order corresponds to the order in 'populations'.
    'K_ext': np.array([1600, 1500, 1500, 1500, 2100, 1900, 1900, 2000, 1900, 1900, 2900, 2100, 2100]),  # test
    # 'K_ext': np.array([2000, PV_ext_scan, SOM_ext_scan, VIP_ext_scan, 2000, PV_ext_scan, SOM_ext_scan, 2000, PV_ext_scan, SOM_ext_scan, 2000, PV_ext_scan, SOM_ext_scan]),    # layer-independent
    # Factor to scale the indegrees.
    'K_scaling': 1.0,
    # Factor to scale the number of neurons.
    'N_scaling': 1.0,
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
        'E_L': {'default': -67.0, 'PC': -63.3, 'PV': -66.8, 'SOM': -61.6, 'VIP': -65.7}, #-67.0,
        # Threshold potential of the neurons (in mV).
        'V_th': {'default': -40.0, 'PC': -41.0, 'PV': -40.5, 'SOM': -40.3, 'VIP': -41.2}, #-40.0,  # -50.0
        # Membrane potential after a spike (in mV).
        'V_reset': -67.0, #-65.0,
        # Membrane capacitance (in pF).
        'C_m': {'default': 200.0, 'PC': 322.0, 'PV': 86.2, 'SOM': 134.0, 'VIP': 86.5}, #200.0, #250.0,
        # Membrane time constant (in ms).
        'tau_m': {'default': 10.0, 'PC': 13.0, 'PV': 3.6, 'SOM': 11.8, 'VIP': 10.9}, #7.0, #10.0,
        # Time constant of postsynaptic excitatory currents (in ms).
        'tau_syn_ex': 1.74, #1.0, # 0.5,
        # Time constant of postsynaptic inhibitory currents (in ms).
        'tau_syn_in': 4.6, #2.0, # 0.5,
        # Time constant of external postsynaptic excitatory current (in ms).
        'tau_syn_E': 1.74,  # 0.5,
        # Refractory period of the neurons after a spike (in ms).
        't_ref': 2.0}
    }

updated_dict = {
    # PSP mean matrix.
    'PSP_mean_matrix': get_mean_PSP_matrix(
        net_dict['PSP_e'], net_dict['g'], len(net_dict['populations'])
        ),
    # PSP std matrix.
    'PSP_std_matrix': get_std_PSP_matrix(
        net_dict['PSP_sd'], len(net_dict['populations'])
        ),
    # mean delay matrix.
    'mean_delay_matrix': get_mean_delays(
        net_dict['mean_delay_exc'], net_dict['mean_delay_inh'],
        len(net_dict['populations'])
        ),
    # std delay matrix.
    'std_delay_matrix': get_std_delays(
        net_dict['mean_delay_exc'] * net_dict['rel_std_delay'],
        net_dict['mean_delay_inh'] * net_dict['rel_std_delay'],
        len(net_dict['populations'])
        ),
    }


net_dict.update(updated_dict)
