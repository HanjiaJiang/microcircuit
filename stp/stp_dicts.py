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