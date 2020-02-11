from microcircuit.network_params import net_dict
import copy
cell_types = ['Exc', 'PV', 'SOM', 'VIP']

'''
Static synapses
'''
static_template = {
    'model': 'static_synapse'
}
no_stp = {}
for i, pre_type in enumerate(cell_types):
    no_stp[pre_type] = {}
    for j, post_type in enumerate(cell_types):
        no_stp[pre_type][post_type] = copy.deepcopy(static_template)

'''
Allen data
'''
allen_stp = {}
# relative changes from 1st to 5th PSP
allen_1to5 = \
    [[-0.273996, -0.455301, -0.137356, -0.008266],
     [-0.185856, -0.365362, -0.130458, -0.066462],
     [0.165965, -0.423736, -0.147765, -0.096198],
     [-0.007874, -0.327010, 0.132345, -0.128896]]
tau_arr = [
    [16, 18, 0.1, 0.1],
    [16, 24, 7, 0.1],
    [25, 16, 0.1, 0.1],
    [0.1, 9, 57, 0.1],
]
stp_dict_template = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 0.01,
}
for pre_type in cell_types:
    allen_stp[pre_type] = {}
for i, post_type in enumerate(cell_types):
    for j, pre_type in enumerate(cell_types):
        tmp_dict = copy.deepcopy(stp_dict_template)
        if pre_type != 'Exc':
            tmp_dict['tau_psc'] = net_dict['neuron_params']['tau_syn_in']
        #     if post_type == 'VIP':
        #         tmp_dict = copy.deepcopy(static_template)
        # if pre_type == 'VIP':
        #     tmp_dict = copy.deepcopy(static_template)
        if 'tau_fac' in tmp_dict:
            if allen_1to5[i][j] >= 0:
                tmp_dict['tau_fac'] = tau_arr[i][j]
            else:
                tmp_dict['tau_rec'] = tau_arr[i][j]
        # print(tmp_dict)
        allen_stp[pre_type][post_type] = copy.deepcopy(tmp_dict)
# print(allen_stp)

'''
Doiron data
'''
exc_strong = {
    'model': 'tsodyks_synapse',
    'U': 0.75,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 800.0,
}

pv_strong = {
    'model': 'tsodyks_synapse',
    'U': 0.9,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_in'],
    'tau_rec': 800.0,
}

exc_weak = {
    'model': 'tsodyks_synapse',
    'U': 0.75,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 100.0,
}

pv_weak = {
    'model': 'tsodyks_synapse',
    'U': 0.9,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_in'],
    'tau_rec': 100.0,
}

som_fac = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 200.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 0.01,
}

doiron_stp = {
    'Exc': {
        'Exc': exc_strong,
        'PV': exc_strong,
        'SOM': som_fac
    },
    'PV': {
        'Exc': pv_strong,
        'PV': pv_strong,
        'SOM': pv_strong,
        'VIP': pv_strong
    }
}

doiron_stp_weak = {
    'Exc': {
        'Exc': exc_weak,
        'PV': exc_weak,
        'SOM': som_fac
    },
    'PV': {
        'Exc': pv_weak,
        'PV': pv_weak,
        'SOM': pv_weak,
        'VIP': pv_weak
    }
}

'''
For testing in ins_models.py
'''
dep_syn = {
    'model': 'tsodyks_synapse',
    'U': 0.9,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 10.0,
}

fac_syn = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 10.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 0.01,
}

test_stp = {
    'Exc': {
        'Exc': dep_syn,
        'PV': dep_syn,
        'SOM': fac_syn,
        'VIP': static_template
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
        'Exc': static_template,
        'PV': static_template,
        'SOM': static_template,
        'VIP': dep_syn
    },
}