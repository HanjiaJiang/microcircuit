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
self-defined
'''
dep_e_low = {
    'model': 'tsodyks_synapse',
    'U': 0.75,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 50.0,
}
dep_e_high = {
    'model': 'tsodyks_synapse',
    'U': 0.75,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 200.0,
}
dep_i_low = {
    'model': 'tsodyks_synapse',
    'U': 0.9,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_in'],
    'tau_rec': 50.0,
}
dep_i_high = {
    'model': 'tsodyks_synapse',
    'U': 0.9,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_in'],
    'tau_rec': 200.0,
}

fac_e_low = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 100.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 0.01,
}
fac_e_high = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 200.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 0.01,
}
fac_i_low = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 100.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_in'],
    'tau_rec': 0.01,
}
fac_i_high = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 200.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_in'],
    'tau_rec': 0.01,
}
custom_stp = {
    'Exc': {
        'Exc': dep_e_high,
        'PV': static_template,
        'SOM': fac_e_high,
        'VIP': static_template
    },
    'PV': {
        'Exc': static_template,
        'PV': static_template,
        'SOM': static_template,
        'VIP': static_template
    },
    'SOM': {
        'Exc': static_template,
        'PV': static_template,
        'SOM': static_template,
        'VIP': static_template
    },
    'VIP': {
        'Exc': static_template,
        'PV': static_template,
        'SOM': static_template,
        'VIP': static_template
    },
}

'''
BBP data
'''
e1 = {
    'model': 'tsodyks_synapse',
    'U': 0.09,
    'tau_rec': 138.0,
    'tau_fac': 670.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex']
}
e2 = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_rec': 671.0,
    'tau_fac': 17.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex']
}
e3 = {
    'model': 'tsodyks_synapse',
    'U': 0.29,
    'tau_rec': 329.0,
    'tau_fac': 326.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex']
}
i1 = {
    'model': 'tsodyks_synapse',
    'U': 0.016,
    'tau_rec': 45.0,
    'tau_fac': 376.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex']
}
i2 = {
    'model': 'tsodyks_synapse',
    'U': 0.25,
    'tau_rec': 706.0,
    'tau_fac': 21.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex']
}
i3 = {
    'model': 'tsodyks_synapse',
    'U': 0.32,
    'tau_rec': 144.0,
    'tau_fac': 62.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex']
}
bbp_stp = {
    'Exc': {
        'Exc': e2,
        'PV': e2,
        'SOM': e1,
        'VIP': e1
    },
    'PV': {
        'Exc': i3,
        'PV': i3,
        'SOM': i3,
        'VIP': i3
    },
    'SOM': {
        'Exc': i2,
        'PV': i2,
        'SOM': i2,
        'VIP': i2
    },
    'VIP': {
        'Exc': i2,
        'PV': i2,
        'SOM': i2,
        'VIP': i2
    },
}

def create_neuron(subtype, n_dict):
    nid = nest.Create(n_dict['neuron_model'])
    nest.SetStatus(nid, {
        'tau_syn_ex': n_dict['neuron_params']['tau_syn_ex'],
        'tau_syn_in': n_dict['neuron_params']['tau_syn_in'],
        'E_L': n_dict['neuron_params']['E_L'][subtype],
        'V_th': n_dict['neuron_params']['V_th'][subtype],
        'C_m': n_dict['neuron_params']['C_m'][subtype],
        'tau_m': n_dict['neuron_params']['tau_m'][subtype],
        'V_reset':  n_dict['neuron_params']['V_reset'],
        't_ref': n_dict['neuron_params']['t_ref']
        })
    return nid

if __name__ == '__main__':
    # create cells and connections
    pre_cells = []
    post_cells = []
    for cell_type in cell_types:
        pre_cells.append(create_neuron(cell_type, net_dict))
        post_cells.append(create_neuron(cell_type, net_dict))
    for i, pre_type in enumerate(cell_types):
        for j, post_type in enumerate(cell_types):
            nest.Connect(pre_cells[i], post_cells[i], syn_spec=custom_stp[pre_type][post_type])
    # create
