import copy
import pickle

# load fitted STPs
fitted_stp = {}
try:
    with open('microcircuit/stp/stp_fitted.pickle', 'rb') as h:
        fitted_stp = pickle.load(h)
except FileNotFoundError:
    print('stp_dicts.py: stp_fitted.pickle not found')

'''
Doiron data
'''
doiron_e = {
    'model': 'tsodyks_synapse',
    'U': 0.75,
    'tau_fac': 0.0,
    'tau_rec': 800.0,
}

doiron_pv = {
    'model': 'tsodyks_synapse',
    'U': 0.9,
    'tau_fac': 0.0,
    'tau_rec': 800.0,
}

doiron_e_weak = {
    'model': 'tsodyks_synapse',
    'U': 0.75,
    'tau_fac': 0.0,
    'tau_rec': 100.0,
}

doiron_pv_weak = {
    'model': 'tsodyks_synapse',
    'U': 0.9,
    'tau_fac': 0.0,
    'tau_rec': 100.0,
}

doiron_e2som = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 200.0,
    'tau_rec': 0.01,
}

doiron_stp = {
    'Exc': {
        'Exc': doiron_e,
        'PV': doiron_e,
        'SOM': doiron_e2som
    },
    'PV': {
        'Exc': doiron_pv,
        'PV': doiron_pv,
        'SOM': doiron_pv,
        'VIP': doiron_pv
    }
}

doiron_stp_weak = {
    'Exc': {
        'Exc': doiron_e_weak,
        'PV': doiron_e_weak,
        'SOM': doiron_e2som
    },
    'PV': {
        'Exc': doiron_pv_weak,
        'PV': doiron_pv_weak,
        'SOM': doiron_pv_weak,
        'VIP': doiron_pv_weak
    }
}

'''
BBP data
'''
e1 = {
    'model': 'tsodyks_synapse',
    'U': 0.09,
    'tau_rec': 138.0,
    'tau_fac': 670.0
}
e2 = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_rec': 671.0,
    'tau_fac': 17.0
}
e3 = {
    'model': 'tsodyks_synapse',
    'U': 0.29,
    'tau_rec': 329.0,
    'tau_fac': 326.0
}
i1 = {
    'model': 'tsodyks_synapse',
    'U': 0.016,
    'tau_rec': 45.0,
    'tau_fac': 376.0
}
i2 = {
    'model': 'tsodyks_synapse',
    'U': 0.25,
    'tau_rec': 706.0,
    'tau_fac': 21.0
}
i3 = {
    'model': 'tsodyks_synapse',
    'U': 0.32,
    'tau_rec': 144.0,
    'tau_fac': 62.0
}
bbp_stp = {
    'Exc': {
        'Exc': e2,
        'PV': e2,
        'SOM': e1,
        'VIP': e1
    },
    'PV': {
        'Exc': i3
    },
    'SOM': {
        'Exc': i2
    },
    'VIP': {
        'Exc': i2
    },
}
