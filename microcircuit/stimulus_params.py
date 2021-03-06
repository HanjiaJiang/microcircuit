import numpy as np
from microcircuit.network_params import net_dict
from microcircuit.sim_params import sim_dict

stim_dict = {
    # Turn thalamic input on or off (True or False).
    'thalamic_input': False,
    # Turn DC input on or off (True or False).
    'dc_input': False,
    # Number of thalamic neurons.
    'n_thal': 200,  # Constantinople, Bruno, 2013; Oberlaender et al., 2011
    # Mean amplitude of the thalamic postsynaptic potential (in mV).
    'PSP_th': 0.49, #0.15,  # Bruno, 2006, Cortex Is Driven...
    # Standard deviation of the postsynaptic potential (in relative units).
    'PSP_sd': 0.27, #0.1,   # Bruno, 2006, Cortex Is Driven...
    # Start of the thalamic input (in ms).
    'th_start': np.array([np.nan]),
    # 'th_start': np.arange(2000.0, sim_dict['t_sim'], stim_duration*2),
    # Duration of the thalamic input (in ms).
    'th_duration': 10.0,
    # Rate of the thalamic input (in Hz).
    'th_rate': 120.0,
    # Start of the DC generator (in ms).
    'dc_start': 0.0,
    # Duration of the DC generator (in ms).
    'dc_dur': 1000.0,
    # Connection probabilities of the thalamus to the different populations.
    # Order as in 'populations' in 'network_params.py'
    'conn_probs_th':
        np.array([0.0, 0.0, 0.0, 0.0, 0.0983, 0.0619, 0.0, 0.0, 0.0, 0.0, 0.0512, 0.0196, 0.0]),
    # Mean delay of the thalamic input (in ms).
    'delay_th':
        np.asarray([1.5 for i in list(range(len(net_dict['populations'])))]),
    # Standard deviation of the thalamic delay (in ms).
    'delay_th_sd':
        np.asarray([0.75 for i in list(range(len(net_dict['populations'])))]),
    # Amplitude of the DC generator (in pA).
    'dc_amp': np.ones(len(net_dict['populations'])) * 0.3,
    # stimulus orientation 190614
    'orientation': 0.0,  # from -pi/2 to pi/2
    # perturbation: for inhibition vs. disinhibition, paradox effects, etc.
    'perturbs': {
            'type': 'dc',
            'levels': [0.],     # strength or offset levels
            'duration': 1000.,  # stimulus duration
            'n_repeat': 0,      # n of repeat for each level
            'targets': [1],     # target populations
            'interval': 1000.,  # inter-stimulus interval
            'starts': {'0.0': [2000.]},
            'ac':{
                'amplitude': 0.1,
                'frequency': 100.,
                'phase': 0.
            }
        },
    }
