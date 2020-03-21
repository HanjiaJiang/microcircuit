import nest
import numpy as np
import copy
from stp.stp_dicts import cell_types, allen_stp, doiron_stp, doiron_stp_weak
from microcircuit.conn import conn_barrel_integrate
from microcircuit.raw_data import rat_dict, mouse_dict, bbp, exp, dia, allen, dia_allen
np.set_printoptions(precision=2, linewidth=500, suppress=True)

'''
Max firing rate:
    use refractory period t_ref to constrain firing rate.

Short-term plasticity:
    Tsodyks synapse, see
    Ashok Litwin-Kumar, Robert Rosenbaum, Brent Doiron, 2016
    Tsodyks et al. 1997

Relative strengths:
    adjustment of relative PV and SOM inhibitory strength

Selectivity:
    orientation tuning by self-defined cosine functions

Cell-type specific parameters:
    Garrett T. Neske, Saundra L. Patrick, and Barry W. Connors, 2015
'''

special_dict = {
    'fmax': False,
    # cell-type specific parameters
    'ctsp': True,
    # STP
    'stp_dict': doiron_stp_weak,
    # relative inhibitory strengths
    'adjust_inh': True,
    'som_power': 1.0,
    'pv_power': 1.0,
    # selectivity
    'orient_tuning': False,
    'sel_inh_src': ['PV', 'SOM'],
    'sel_inh_trg': ['PV', 'SOM'],
    'k_th': 0.8,
    'k_e2e': 0.8,
    'k_e2i': 0.8,
    'k_i2e': 0.2
}


def assign_stp(source_name, target_name, weight_dict, delay_dict, stp_dict):
    syn_dict = {
        'model': 'static_synapse',
        'weight': weight_dict,
        'delay': delay_dict
    }
    for pre_type in stp_dict.keys():
        for post_type in stp_dict[pre_type].keys():
            if pre_type in source_name and post_type in target_name:
                syn_dict = copy.deepcopy(stp_dict[pre_type][post_type])
                syn_dict['weight'] = weight_dict
                syn_dict['delay'] = delay_dict
    # print(source_name, target_name)
    # print(syn_dict)
    return syn_dict


def set_fmax(names, population, spe_dict):
    # Fmax:
    # Two Dynamically Distinct Inhibitory Networks in
    # Layer 4 of the Neocortex
    if spe_dict['fmax'] is True:
        if 'Exc' in names:
            nest.SetStatus(population, {'t_ref': 35.7})
            # if 'L5' in pop:
            #     nest.SetStatus(population, {'I_e': 375.0})
            # elif 'L2' in pop:
            #     nest.SetStatus(population, {'I_e': 300.0})
        elif 'PV' in names:
            nest.SetStatus(population, {'t_ref': 6.2})
        elif 'SOM' in names:
            nest.SetStatus(population, {'t_ref': 9.5})
        # The Largest Group of Superficial Neocortical GABAergic
        # Interneurons Expresses Ionotropic Serotonin Receptors
        elif 'VIP' in names:
            nest.SetStatus(population, {'t_ref': 15.3})


def inh_weight(source_name, weight, spe_dict):
    # HJ: Short-Term Plasticity of Unitary Inhibitory-to-Inhibitory
    # Synapses Depends on the Presynaptic Interneuron Subtype
    if spe_dict['adjust_inh'] is True:
        if 'SOM' in source_name:
            weight *= spe_dict['som_power']
        if 'PV' in source_name:
            weight *= spe_dict['pv_power']
        return weight


# clumsy, to be improved...
def set_thalamus_input(th_pop,
                       poisson_pops,
                       start_times,
                       stop_times,
                       rate_0,
                       stim_theta,
                       spe_dict):
    # tuning
    if spe_dict['orient_tuning'] is True:
        # each poisson population is for 1 stimulus (in e.g. 20 repetitions)
        for i, pop in enumerate(poisson_pops):
            # within each population,
            # rates are distributed according to stimulus angle
            for j, node in enumerate(pop):
                theta = -np.pi / 2.0 + np.pi * ((j + 0.5) / float(len(pop)))
                rate = rate_0 * (1.0 + spe_dict['k_th']
                                 * np.cos(2.0 * (theta - stim_theta)))
                nest.SetStatus([node], {
                    'rate': rate,
                    'start': start_times[i],
                    'stop': stop_times[i]
                })
            # and then connected to the thalamus (th_pop) one by one
            nest.Connect(pop, th_pop, conn_spec={'rule': 'one_to_one'})
    # no tuning
    else:
        for i, pop in enumerate(poisson_pops):
            nest.SetStatus(pop, {
                'rate': rate_0,
                'start': start_times[i],
                'stop': stop_times[i]
            })
            nest.Connect(pop, th_pop, conn_spec={'rule': 'one_to_one'})


def connect_thalamus_orientation(th_pop,
                                 target_pop,
                                 target_name,
                                 nr_synapses,
                                 syn_dict_th,
                                 spe_dict,
                                 bernoulli_prob=None):
    # do it only if connection is not 0
    if nr_synapses > 0:
        if spe_dict['orient_tuning'] and 'Exc' in target_name:
            nr_cluster = 8
            len_th = len(th_pop)
            len_target = len(target_pop)
            p_0 = nr_synapses / float(len_th * len_target)

            # thalamus clustering
            len_th_cluster = int(len_th / nr_cluster)
            if len_target % nr_cluster != 0:
                len_th_cluster += 1
            th_cluster_list = []
            for x in range(nr_cluster):
                head_idx = x * len_th_cluster
                tail_idx = (x + 1) * len_th_cluster
                if tail_idx <= len(th_pop):
                    th_cluster_list.append(th_pop[head_idx:tail_idx])
                else:
                    th_cluster_list.append(th_pop[head_idx:])

            # target clustering
            len_target_cluster = int(len_target / nr_cluster)
            if len_target % nr_cluster != 0:
                len_target_cluster += 1
            target_cluster_list = []
            for y in range(nr_cluster):
                head_idx = y * len_target_cluster
                tail_idx = (y + 1) * len_target_cluster
                # print('head,tail={0},{1}'.format(head_idx,tail_idx))
                if tail_idx <= len(target_pop):
                    target_cluster_list.append(target_pop[head_idx:tail_idx])
                else:
                    target_cluster_list.append(target_pop[head_idx:])

            for i, th_cluster in enumerate(th_cluster_list):
                for j, target_cluster in enumerate(target_cluster_list):
                    theta_th = -np.pi / 2.0 + np.pi * ((i + 0.5) / float(nr_cluster))
                    theta_target = -np.pi / 2.0 + np.pi * ((j + 0.5) / float(nr_cluster))
                    p = p_0 * (1.0 + spe_dict['k_e2e'] * np.cos(2.0 * (theta_th - theta_target)))
                    conn_nr = int(round(len(th_cluster) * len(target_cluster) * p))
                    # print('theta_th={0}, theta_target={1}, conn_nr={2}'
                    #       .format(theta_th, theta_target, conn_nr))
                    nest.Connect(th_cluster, target_cluster,
                                 conn_spec={
                                     'rule': 'fixed_total_number',
                                     'N': conn_nr,
                                 },
                                 syn_spec=syn_dict_th
                                 )
        else:
            # connect by probability (Bernoulli)
            if bernoulli_prob is not None:
                nest.Connect(th_pop, target_pop,
                conn_spec={'rule': 'pairwise_bernoulli', 'p': bernoulli_prob},
                syn_spec=syn_dict_th)
                # print('p={}'.format(bernoulli_prob))
            # connect by synapse number
            else:
                nest.Connect(
                    th_pop, target_pop,
                    conn_spec={
                        'rule': 'fixed_total_number',
                        'N': nr_synapses,
                    },
                    syn_spec=syn_dict_th
                )
                # print('n={}'.format(nr_synapses))


def connect_by_cluster(source_name,
                       target_name,
                       synapse_nr,
                       syn_dict,
                       source_pop,
                       target_pop,
                       spe_dict,
                       conn_prob=None):
    # For orientation clustering
    nr_cluster = 8  # number of clusters
    k = 0.0     # modulation constant
    if 'Exc' in source_name and 'Exc' in target_name:
        k = spe_dict['k_e2e']
    for inh_target in spe_dict['sel_inh_trg']:
        if 'Exc' in source_name and inh_target in target_name:
            k = spe_dict['k_e2i']
    for inh_source in spe_dict['sel_inh_src']:
        if inh_source in source_name and 'Exc' in target_name:
            k = spe_dict['k_i2e']
    if spe_dict['orient_tuning'] is True and k != 0.0:
        # do it only if connection is not 0
        if synapse_nr > 0:
            len_source = len(source_pop)
            len_target = len(target_pop)
            p_0 = synapse_nr / float(len_source * len_target)

            # source clustering
            len_source_cluster = int(len_source / nr_cluster)
            if len_source % nr_cluster != 0:
                len_source_cluster += 1
                print(source_name + ' source not exact, len_source={}'.
                      format(len_source))
            source_cluster_list = []
            for x in range(nr_cluster):
                head_idx = x * len_source_cluster
                tail_idx = (x + 1) * len_source_cluster
                if tail_idx <= len(source_pop):
                    source_cluster_list.append(source_pop[head_idx:tail_idx])
                else:
                    source_cluster_list.append(source_pop[head_idx:])

            # target clustering
            len_target_cluster = int(len_target / nr_cluster)
            if len_target % nr_cluster != 0:
                len_target_cluster += 1
                print(target_name + ' target not exact, len_target={}'.format(len_target))
            target_cluster_list = []
            for y in range(nr_cluster):
                head_idx = y * len_target_cluster
                tail_idx = (y + 1) * len_target_cluster
                if tail_idx <= len(target_pop):
                    target_cluster_list.append(target_pop[head_idx:tail_idx])
                else:
                    target_cluster_list.append(target_pop[head_idx:])

            conn_nr_sum = 0
            for i, source_cluster in enumerate(source_cluster_list):
                for j, target_cluster in enumerate(target_cluster_list):
                    theta_source = -np.pi / 2.0 + np.pi * ((i + 0.5) / float(nr_cluster))
                    theta_target = -np.pi / 2.0 + np.pi * ((j + 0.5) / float(nr_cluster))
                    p = p_0 * (1.0 + k * np.cos(2.0 * (theta_source - theta_target)))
                    conn_nr = int(round(len(source_cluster) * len(target_cluster) * p))
                    # print('theta_source={:.2f}, theta_target={:.2f}, conn_nr={:.2f}'
                    #       .format(theta_source, theta_target, conn_nr))
                    nest.Connect(source_cluster, target_cluster,
                                 conn_spec={
                                     'rule': 'fixed_total_number',
                                     'N': conn_nr,
                                 },
                                 syn_spec=syn_dict
                                 )
                    conn_nr_sum += conn_nr
            print(source_name + ' to ' + target_name + ' verify synapse_nr: {0}, {1}'.format(
                synapse_nr, conn_nr_sum))
    else:
        if isinstance(conn_prob, float):
            conn_dict_rec = {
                'rule': 'pairwise_bernoulli', 'p': conn_prob
            }
            synapse_nr = len(source_pop)*len(target_pop)*conn_prob
        else:
            conn_dict_rec = {
                'rule': 'fixed_total_number', 'N': synapse_nr
            }
        nest.Connect(
            source_pop, target_pop,
            conn_spec=conn_dict_rec,
            syn_spec=syn_dict
        )
    return synapse_nr


def ctsp_assign(pop, net_dict, spe_dict):
    E_L = net_dict['neuron_params']['E_L']['default']
    V_th = net_dict['neuron_params']['V_th']['default']
    C_m = net_dict['neuron_params']['C_m']['default']
    tau_m = net_dict['neuron_params']['tau_m']['default']
    if spe_dict['ctsp'] is True:
        for celltype in ['Exc', 'PV', 'SOM', 'VIP']:
            if celltype in pop:
                E_L = net_dict['neuron_params']['E_L'][celltype]
                V_th = net_dict['neuron_params']['V_th'][celltype]
                C_m = net_dict['neuron_params']['C_m'][celltype]
                tau_m = net_dict['neuron_params']['tau_m'][celltype]
                break
    return E_L, V_th, C_m, tau_m


def get_weight(psp_val, net_dict):
    C_m = net_dict['neuron_params']['C_m']['default']
    tau_m = net_dict['neuron_params']['tau_m']['default']
    tau_syn_ex = net_dict['neuron_params']['tau_syn_ex']

    PSC_e_over_PSP_e = (((C_m) ** (-1) * tau_m * tau_syn_ex / (
        tau_syn_ex - tau_m) * ((tau_m / tau_syn_ex) ** (
            - tau_m / (tau_m - tau_syn_ex)) - (tau_m / tau_syn_ex) ** (
                - tau_syn_ex / (tau_m - tau_syn_ex)))) ** (-1))
    PSC_e = (PSC_e_over_PSP_e * psp_val)
    return PSC_e


# calculate celltype-specific psc
def calc_psc(psp_val, C_m, tau_m, tau_syn):
    PSC_e_over_PSP_e = (((C_m) ** (-1) * tau_m * tau_syn / (
        tau_syn - tau_m) * ((tau_m / tau_syn) ** (
            - tau_m / (tau_m - tau_syn)) - (tau_m / tau_syn) ** (
                - tau_syn / (tau_m - tau_syn)))) ** (-1))
    PSC_e = (PSC_e_over_PSP_e * psp_val)
    return PSC_e


# get the psc matrix
def get_weights(net_dict, dim=13, lyr_gps=None):
    if lyr_gps is None:
        lyr_gps = [[0, 1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]
    pscs = np.zeros((dim, dim))
    for i in range(dim):    # post
        for j in range(dim):    # pre
            a = 0
            b = 0
            for idx, gps in enumerate(lyr_gps):
                if i in gps:
                    a = idx
                if j in gps:
                    b = idx
            psp = net_dict['w_dict']['psp_mtx'][a, b]
            for celltype in ['Exc', 'PV', 'SOM', 'VIP']:
                if celltype in net_dict['populations'][i]:
                    psc = calc_psc(psp,
                    net_dict['neuron_params']['C_m']['default'],
                    net_dict['neuron_params']['tau_m']['default'],
                    net_dict['neuron_params']['tau_syn_ex'])
                    if j in [0, 4, 7, 10]:
                        pscs[i, j] = psc
                    else:
                        pscs[i, j] = psc*net_dict['g']
    return pscs


# get the psc std matrix
def get_weight_stds(net_dict, dim=13, lyr_gps=None):
    if lyr_gps is None:
        lyr_gps = [[0, 1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]
    std_ratios = np.zeros((dim, dim))
    # print(net_dict['psp_std_mtx'])
    for i, pre_lyr in enumerate(lyr_gps):
        for j, post_lyr in enumerate(lyr_gps):
            std_ratio = net_dict['w_dict']['psp_std_mtx'][j, i]
            std_ratios[post_lyr[0]:post_lyr[-1]+1, pre_lyr[0]:pre_lyr[-1]+1] = std_ratio
    return std_ratios


def eq_inh_conn(n_full, conn, lyr_gps=None):
    if lyr_gps is None:
        lyr_gps = [[0, 1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]
    conn_out = np.zeros(conn.shape)
    for lyr_pre in lyr_gps:
        for lyr_post in lyr_gps:
            # connectivity to be calculated
            e2e = e2i = i2e = i2i = 0.0
            e2i_max = e2i_sum = 0.0
            i2e_max = i2e_sum = 0.0
            i2i_max = i2i_sum = 0.0

            # index of pre and post excitatory
            exc_pre = lyr_pre[0]
            exc_post = lyr_post[0]

            # re-calculate connectivity
            # i to e
            for inh_pre in lyr_pre[1:]:
                i2e_max += n_full[exc_post]*n_full[inh_pre]
                i2e_sum += n_full[exc_post]*n_full[inh_pre]*conn[exc_post, inh_pre]

            # e to i
            for inh_post in lyr_post[1:]:
                e2i_max += n_full[exc_pre]*n_full[inh_post]
                e2i_sum += n_full[exc_pre]*n_full[inh_post]*conn[inh_post, exc_pre]

            # i to i
            for inh_pre in lyr_pre[1:]:
                for inh_post in lyr_post[1:]:
                    i2i_max += n_full[inh_pre]*n_full[inh_post]
                    i2i_sum += n_full[inh_pre]*n_full[inh_post]*conn[inh_post, inh_pre]

            e2e = conn[exc_post, exc_pre]
            if e2i_max != 0.0:
                e2i = e2i_sum/e2i_max
            if i2e_max != 0.0:
                i2e = i2e_sum/i2e_max
            if i2i_max != 0.0:
                i2i = i2i_sum/i2i_max

            # re-assign connectivity
            conn_out[exc_post, exc_pre] = e2e
            for inh_pre in lyr_pre[1:]:
                conn_out[exc_post, inh_pre] = i2e
            for inh_post in lyr_post[1:]:
                conn_out[inh_post, exc_pre] = e2i
            for inh_pre in lyr_pre[1:]:
                for inh_post in lyr_post[1:]:
                    conn_out[inh_post, inh_pre] = i2i

    return conn_out


# original helpers.py functions
def compute_DC(net_dict, w_ext):
    DC = (
        net_dict['bg_rate'] * net_dict['K_ext'] *
        w_ext * net_dict['neuron_params']['tau_syn_E'] * 0.001
        )
    return DC


def get_total_number_of_synapses(net_dict):
    N_full = net_dict['N_full']
    number_N = len(N_full)
    conn_probs = net_dict['conn_probs']
    prod = np.outer(N_full, N_full)

    #
    if net_dict['renew_conn'] is True:
        if net_dict['animal'] == 'mouse':
            animal_dict = mouse_dict
        else:
            animal_dict = rat_dict
        conn_probs = conn_barrel_integrate(animal_dict, bbp, exp, allen, dia_allen)

    n_syn_temp = np.log(1. - conn_probs)/np.log((prod - 1.) / prod)
    N_full_matrix = np.column_stack((N_full for i in list(range(number_N))))
    K = (((n_syn_temp * (N_full_matrix).astype(int)) / N_full_matrix).astype(int))
    return K


def synapses_th_matrix(net_dict, stim_dict):
    N_full = net_dict['N_full']
    conn_probs = stim_dict['conn_probs_th']
    T_full = stim_dict['n_thal']
    prod = (T_full * N_full).astype(float)
    n_syn_temp = np.log(1. - conn_probs)/np.log((prod - 1.)/prod)
    K = (((n_syn_temp * (N_full).astype(int))/N_full).astype(int))
    return K
