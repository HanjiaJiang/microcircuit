from microcircuit.network_params import net_dict
import copy
cell_types = ['Exc', 'PV', 'SOM', 'VIP']

# Allen data
allen_stp_dict = {}
allen_stps = \
    [[-0.273996, -0.455301, -0.137356, -0.008266],
     [-0.185856, -0.365362, -0.130458, -0.066462],
     [0.165965, -0.423736, -0.147765, -0.096198],
     [-0.007874, -0.327010, 0.132345, -0.128896]]
tau_arr = [
    [21, 30, 9, 2],
    [16, 27, 10, 7],
    [15, 27, 10, 7],
    [4, 20, 12, 9],
]
stp_dict_template = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 0.01,
}
for pre_type in cell_types:
    allen_stp_dict[pre_type] = {}

for i, post_type in enumerate(cell_types):
    for j, pre_type in enumerate(cell_types):
        tmp_dict = copy.deepcopy(stp_dict_template)
        if pre_type != 'Exc':
            tmp_dict['tau_psc'] = net_dict['neuron_params']['tau_syn_in']
        if allen_stps[i][j] >= 0:
            tmp_dict['tau_fac'] = tau_arr[i][j]
        else:
            tmp_dict['tau_rec'] = tau_arr[i][j]
        # print(tmp_dict)
        allen_stp_dict[pre_type][post_type] = copy.deepcopy(tmp_dict)
# print(allen_stp_dict)

# Doiron data
doiron_dep_exc = {
    'model': 'tsodyks_synapse',
    'U': 0.75,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 800.0,
}

doiron_dep_pv = {
    'model': 'tsodyks_synapse',
    'U': 0.9,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_in'],
    'tau_rec': 800.0,
}

doiron_fac_exc2som = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 200.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 0.01,
}

doiron_stp = {
    'Exc': {
        'Exc': doiron_dep_exc,
        'PV': doiron_dep_exc,
        'SOM': doiron_fac_exc2som
    },
    'PV': {
        'Exc': doiron_dep_pv,
        'PV': doiron_dep_pv,
        'SOM': doiron_dep_pv,
        'VIP': doiron_dep_pv
    }
}

# test stp
U_test = 0.5
dep_syn = {
    'model': 'tsodyks_synapse',
    'U': U_test,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 10.0,
}

fac_syn = {
    'model': 'tsodyks_synapse',
    'U': U_test,
    'tau_fac': 10.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 0.01,
}

test_stp = {
    'Exc': {
        'Exc': dep_syn,
        'PV': dep_syn,
        'SOM': fac_syn,
        'VIP': dep_syn
    },
    'PV': {
        'Exc': dep_syn,
        'PV': dep_syn,
        'SOM': dep_syn,
        'VIP': dep_syn
    },
    'SOM': {
        'Exc': dep_syn,
        'PV': dep_syn,
        'SOM': dep_syn,
        'VIP': fac_syn
    },
    'VIP': {
        'Exc': dep_syn,
        'PV': dep_syn,
        'SOM': dep_syn,
        'VIP': dep_syn
    },
}