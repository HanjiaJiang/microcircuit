from microcircuit.network_params import net_dict

# stp dicts
doiron_dep_exc = {
    'model': 'tsodyks_synapse',
    'U': 0.75,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 15.0,
}

doiron_dep_pv = {
    'model': 'tsodyks_synapse',
    'U': 0.9,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_in'],
    'tau_rec': 15.0,
}

doiron_fac_exc2som = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 30.0,
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
normal_syn = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 10.0,
}

dep_syn = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 0.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 20.0,
}

fac_syn = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 20.0,
    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
    'tau_rec': 0.01,
}

test_stp = {
    'Exc': {
        'Exc': dep_syn,
        'PV': dep_syn,
        'SOM': fac_syn,
        'VIP': normal_syn
    },
    'PV': {
        'Exc': dep_syn,
        'PV': dep_syn,
        'SOM': dep_syn,
        'VIP': dep_syn
    },
    'SOM': {
        'Exc': normal_syn,
        'PV': normal_syn,
        'SOM': normal_syn,
        'VIP': fac_syn
    },
    'VIP': {
        'Exc': normal_syn,
        'PV': normal_syn,
        'SOM': normal_syn,
        'VIP': normal_syn
    },
}