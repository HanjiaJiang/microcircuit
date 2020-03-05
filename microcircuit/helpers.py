import numpy as np
from microcircuit.conn import conn_barrel_integrate
from microcircuit.raw_data import rat_dict, mouse_dict, bbp, exp, dia, allen, dia_allen
np.set_printoptions(precision=4, linewidth=100, suppress=True)

def compute_DC(net_dict, w_ext):
    """ Computes DC input if no Poisson input is provided to the microcircuit.

    Parameters
    ----------
    net_dict
        Parameters of the microcircuit.
    w_ext
        Weight of external connections.

    Returns
    -------
    DC
        DC input, which compensates lacking Poisson input.
    """
    DC = (
        net_dict['bg_rate'] * net_dict['K_ext'] *
        w_ext * net_dict['neuron_params']['tau_syn_E'] * 0.001
        )
    return DC


def get_total_number_of_synapses(net_dict):
    """ Returns the total number of synapses between all populations.

    The first index (rows) of the output matrix is the target population
    and the second (columns) the source population. If a scaling of the
    synapses is intended this is done in the main simulation script and the
    variable 'K_scaling' is ignored in this function.

    Parameters
    ----------
    net_dict
        Dictionary containing parameters of the microcircuit.
    N_full
        Number of neurons in all populations.
    number_N
        Total number of populations.
    conn_probs
        Connection probabilities of the eight populations.
    scaling
        Factor that scales the number of neurons.

    Returns
    -------
    K
        Total number of synapses with
        dimensions [len(populations), len(populations)].

    """
    N_full = net_dict['N_full']
    number_N = len(N_full)
    conn_probs = net_dict['conn_probs']
    scaling = net_dict['N_scaling']
    prod = np.outer(N_full, N_full)

    # HJ
    if net_dict['renew_conn'] is True:
        if net_dict['animal'] == 'mouse':
            animal_dict = mouse_dict
        else:
            animal_dict = rat_dict
        conn_probs = conn_barrel_integrate(animal_dict, bbp, exp, allen, dia_allen)

    n_syn_temp = np.log(1. - conn_probs)/np.log((prod - 1.) / prod)
    N_full_matrix = np.column_stack(
        (N_full for i in list(range(number_N)))
        )
    # If the network is scaled the indegrees are calculated in the same
    # fashion as in the original version of the circuit, which is
    # written in sli.
    K = (((n_syn_temp * (
        N_full_matrix * scaling).astype(int)) / N_full_matrix).astype(int))
    print('synapse numbers with multiple connections:')
    # print(K)
    return K


def synapses_th_matrix(net_dict, stim_dict):
    """ Computes number of synapses between thalamus and microcircuit.

    This function ignores the variable, which scales the number of synapses.
    If this is intended the scaling is performed in the main simulation script.

    Parameters
    ----------
    net_dict
        Dictionary containing parameters of the microcircuit.
    stim_dict
        Dictionary containing parameters of stimulation settings.
    N_full
        Number of neurons in the eight populations.
    number_N
        Total number of populations.
    conn_probs
        Connection probabilities of the thalamus to the eight populations.
    scaling
        Factor that scales the number of neurons.
    T_full
        Number of thalamic neurons.

    Returns
    -------
    K
        Total number of synapses.

    """
    N_full = net_dict['N_full']
    number_N = len(N_full)
    scaling = net_dict['N_scaling']
    conn_probs = stim_dict['conn_probs_th']
    T_full = stim_dict['n_thal']
    prod = (T_full * N_full).astype(float)
    n_syn_temp = np.log(1. - conn_probs)/np.log((prod - 1.)/prod)
    K = (((n_syn_temp * (N_full * scaling).astype(int))/N_full).astype(int))
    return K


def adj_w_ext_to_K(K_full, K_scaling, w, w_from_PSP, DC, net_dict, stim_dict):
    """ Adjustment of weights to scaling is performed.

    The recurrent and external weights are adjusted to the scaling
    of the indegrees. Extra DC input is added to compensate the scaling
    and preserve the mean and variance of the input.

    Parameters
    ----------
    K_full
        Total number of connections between the eight populations.
    K_scaling
        Scaling factor for the connections.
    w
        Weight matrix of the connections of the eight populations.
    w_from_PSP
        Weight of the external connections.
    DC
        DC input to the eight populations.
    net_dict
        Dictionary containing parameters of the microcircuit.
    stim_dict
        Dictionary containing stimulation parameters.
    tau_syn_E
        Time constant of the external postsynaptic excitatory current.
    full_mean_rates
        Mean rates of the eight populations in the full scale version.
    K_ext
        Number of external connections to the eight populations.
    bg_rate
        Rate of the Poissonian spike generator.

    Returns
    -------
    w_new
        Adjusted weight matrix.
    w_ext_new
        Adjusted external weight.
    I_ext
        Extra DC input.

    """
    tau_syn_E = net_dict['neuron_params']['tau_syn_E']
    full_mean_rates = net_dict['full_mean_rates']
    w_mean = w_from_PSP
    K_ext = net_dict['K_ext']
    bg_rate = net_dict['bg_rate']
    w_new = w / np.sqrt(K_scaling)
    I_ext = np.zeros(len(net_dict['populations']))
    x1_all = w * K_full * full_mean_rates
    x1_sum = np.sum(x1_all, axis=1)
    if net_dict['poisson_input']:
        x1_ext = w_mean * K_ext * bg_rate
        w_ext_new = w_mean / np.sqrt(K_scaling)
        I_ext = 0.001 * tau_syn_E * (
            (1. - np.sqrt(K_scaling)) * x1_sum + (
                1. - np.sqrt(K_scaling)) * x1_ext) + DC
    else:
        w_ext_new = w_from_PSP / np.sqrt(K_scaling)
        I_ext = 0.001 * tau_syn_E * (
            (1. - np.sqrt(K_scaling)) * x1_sum) + DC
    return w_new, w_ext_new, I_ext

