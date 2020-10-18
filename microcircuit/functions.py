import os
import nest
import numpy as np
import copy
from microcircuit.conn import conn_barrel_integrate
import microcircuit.raw_data as raw
np.set_printoptions(precision=2, linewidth=500, suppress=True)

'''
Verification
'''
verify_dict = {}
def verify_collect(in_str, tag):
    if tag in verify_dict.keys():
        verify_dict[tag] += in_str
    else:
        verify_dict[tag] = in_str

def verify_print(path=None):
    for key, value in verify_dict.items():
        try:
            fpath = os.path.join(path, 'verify-{}.txt'.format(key))
        except TypeError:
            fpath = 'verify-{}.txt'.format(key)
        with open(fpath, 'w') as f:
            f.write(value)
            f.close()

'''
Main functions
'''
# assign synapse dictionary
def assign_syn(source_name, target_name, w, w_sd, delay, delay_sd, network):
    syn_dict = {'model': 'static_synapse'}
    net_dict = network.net_dict
    stp_dict = network.net_dict['stp_dict']
    resol = network.sim_resolution
    for pre_type in stp_dict.keys():
        for post_type in stp_dict[pre_type].keys():
            if pre_type in source_name and post_type in target_name:
                # overwrite
                syn_dict = copy.deepcopy(stp_dict[pre_type][post_type])
                if 'Exc' in pre_type:
                    syn_dict['tau_psc'] = net_dict['neuron_params']['tau_syn_ex']
                else:
                    syn_dict['tau_psc'] = net_dict['neuron_params']['tau_syn_in']
    # make up for the release probability U
    if net_dict['U-compensate'] and 'U' in syn_dict.keys():
        w /= syn_dict['U']
        w_sd /= syn_dict['U']
    # var_n, mu_n: variance and mean of the underlying normal distribution
    var_n = np.log((w ** 2 + w_sd ** 2) / w ** 2)
    mu_n = np.log(abs(w)) - var_n / 2.
    weight_dict = {
           'distribution': 'lognormal', 'mu': mu_n,
           'sigma': np.sign(w)*np.sqrt(var_n)
       }
    delay_dict = {
            'distribution': 'normal_clipped',
            'mu': delay, 'sigma': delay_sd,
            'low': resol
        }
    syn_dict['weight'] = weight_dict
    syn_dict['delay'] = delay_dict
    # weight recorder
    x = net_dict['populations'].index(source_name)
    y = net_dict['populations'].index(target_name)
    if 'weight_recorder' in net_dict['rec_dev']:
        if network.test==True and (x != 0 or y != 0):
            pass
        else:
            copysynapse = syn_dict['model'] + '_' + source_name
            if copysynapse not in network.copysynapses:
                print(copysynapse)
                nest.CopyModel(syn_dict['model'], copysynapse, {'weight_recorder': network.weight_recorder[x][0]})
                # print('wr = {}'.format(network.weight_recorder[x][0]))
                network.copysynapses.append(copysynapse)
            syn_dict['model'] = copysynapse
    return syn_dict

# set max firing rate by refractory (not using)
def set_fmax(names, population, net_dict):
    # Fmax:
    # Two Dynamically Distinct Inhibitory Networks in
    # Layer 4 of the Neocortex
    if net_dict['fmax'] is True:
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

# set thalamic cells
# clumsy, to be improved...
def set_thalamus(th_pop,
                   poisson_pops,
                   start_times,
                   stop_times,
                   rate_0,
                   stim_theta,
                   net_dict):
    # tuning
    if net_dict['orient_tuning'] is True:
        # limits are +- pi/2
        base_theta = min(stim_theta, np.pi/2)
        base_theta = max(base_theta, -np.pi/2)
        # set conjugating theta to switch in between
        if stim_theta > 0:
            conj_theta = base_theta - np.pi/2
        else:
            conj_theta = base_theta + np.pi/2
        do_conj = False
        # each poisson population is for 1 stimulus (in e.g. 20 repetitions)
        for i, pop in enumerate(poisson_pops):
            # print('set_thalamus_poisson():\npop len={}'.format(len(pop)))
            # switching
            if do_conj:
                current_theta = conj_theta
                do_conj = False
            else:
                current_theta = base_theta
                do_conj = True
            # print('stim no. {}, theta={}'.format(i, current_theta))
            # within each population,
            # rates are distributed according to stimulus angle
            for j, node in enumerate(pop):
                theta = -np.pi / 2.0 + np.pi * ((j + 0.5) / float(len(pop)))
                rate = rate_0 * (1.0 + net_dict['k_th']
                                 * np.cos(2.0 * (theta - current_theta)))
                # print('node {}, rate = {}'.format(j, rate))
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
            # print('set_thalamus_poisson():\npop len={}'.format(len(pop)))
            nest.SetStatus(pop, {
                'rate': rate_0,
                'start': start_times[i],
                'stop': stop_times[i]
            })
            nest.Connect(pop, th_pop)
            # nest.Connect(pop, th_pop, conn_spec={'rule': 'one_to_one'})

# connect thalamocortical
def connect_tc(th_pop,
                 target_pop,
                 target_name,
                 total_conn_nr,
                 syn_dict_th,
                 net_dict,
                 bernoulli_prob=None,
                 nr_cluster=8):
    if net_dict['orient_tuning'] and 'Exc' in target_name:
        # size of population
        len_th = len(th_pop)
        len_target = len(target_pop)

        # determine p_0 by type of connection
        if isinstance(bernoulli_prob, float):
            p_0 = bernoulli_prob
        else:
            p_0 = total_conn_nr / float(len_th * len_target)

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
            if tail_idx <= len(target_pop):
                target_cluster_list.append(target_pop[head_idx:tail_idx])
            else:
                target_cluster_list.append(target_pop[head_idx:])

        # connection
        syn_sum = 0
        for i, th_cluster in enumerate(th_cluster_list):
            for j, target_cluster in enumerate(target_cluster_list):
                theta_th = -np.pi / 2.0 + np.pi * ((i + 0.5) / float(nr_cluster))
                theta_target = -np.pi / 2.0 + np.pi * ((j + 0.5) / float(nr_cluster))
                p = p_0 * (1.0 + net_dict['k_e2e'] * np.cos(2.0 * (theta_th - theta_target)))
                conn_nr = int(round(len(th_cluster) * len(target_cluster) * p))
                if isinstance(bernoulli_prob, float):
                    conn_dict = {'rule': 'pairwise_bernoulli', 'p': p}
                else:
                    conn_dict = {'rule': 'fixed_total_number', 'N': conn_nr}
                syn_sum += conn_nr
                # print('theta_th={}, theta_target={}, p_0={}, p={}'
                #       .format(theta_th, theta_target, p_0, p))
                nest.Connect(th_cluster, target_cluster,
                             conn_spec=conn_dict,
                             syn_spec=syn_dict_th
                             )
        # print('thalamus to {} expected vs. clustered conn. nr.: {} vs., {}'.
        # format(target_name, int(p_0*len_th*len_target), syn_sum))
    else:
        # connect by probability (Bernoulli)
        if isinstance(bernoulli_prob, float):
            nest.Connect(th_pop, target_pop,
            conn_spec={'rule': 'pairwise_bernoulli', 'p': bernoulli_prob},
            syn_spec=syn_dict_th)
            verify_collect('{}, p={:.4f}\n'.format(target_name, bernoulli_prob), 'thalamus')
        # connect by synapse number
        else:
            nest.Connect(
                th_pop, target_pop,
                conn_spec={
                    'rule': 'fixed_total_number',
                    'N': total_conn_nr,
                },
                syn_spec=syn_dict_th
            )

# connect recurrent
def connect_recurrent(source_name,
                       target_name,
                       total_conn_nr,
                       syn_dict,
                       source_pop,
                       target_pop,
                       net_dict,
                       bernoulli_prob=None,
                       nr_cluster=8):
    # get sizes and probability
    len_source = len(source_pop)
    len_target = len(target_pop)
    if isinstance(bernoulli_prob, float):
        p_0 = bernoulli_prob
    else:
        p_0 = total_conn_nr / float(len_source * len_target)

    # get k (modulation constant)
    k = 0.0
    if 'Exc' in source_name and 'Exc' in target_name:
        k = net_dict['k_e2e']
    for inh_target in net_dict['sel_inh_trg']:
        if 'Exc' in source_name and inh_target in target_name:
            k = net_dict['k_e2i']
    for inh_source in net_dict['sel_inh_src']:
        if inh_source in source_name and 'Exc' in target_name:
            k = net_dict['k_i2e']

    # connection
    if net_dict['orient_tuning'] is True and k != 0.0:
        if p_0 != 0.0:
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

            # connection
            conn_nr_sum = 0
            for i, source_cluster in enumerate(source_cluster_list):
                for j, target_cluster in enumerate(target_cluster_list):
                    # calculate p and connection number
                    theta_source = -np.pi / 2.0 + np.pi * ((i + 0.5) / float(nr_cluster))
                    theta_target = -np.pi / 2.0 + np.pi * ((j + 0.5) / float(nr_cluster))
                    p = p_0 * (1.0 + k * np.cos(2.0 * (theta_source - theta_target)))
                    conn_nr = int(round(len(source_cluster) * len(target_cluster) * p))
                    # if 'Exc' in source_name and 'PV' in target_name:
                    #     print('theta_source={:.2f}, theta_target={:.2f}, p={:.4f}'
                    #           .format(theta_source, theta_target, p))

                    # define conn_dict
                    if isinstance(bernoulli_prob, float):
                        conn_dict = {'rule': 'pairwise_bernoulli', 'p': p}
                    else:
                        conn_dict = {'rule': 'fixed_total_number', 'N': conn_nr}
                    nest.Connect(source_cluster, target_cluster,
                                 conn_spec=conn_dict,
                                 syn_spec=syn_dict)
                    conn_nr_sum += conn_nr
            print(source_name + ' to ' + target_name + ' expected vs. clustered conn. nr.: {} vs. {}'.format(
                int(p_0*len_source*len_target), conn_nr_sum))
    else:
        if isinstance(bernoulli_prob, float):
            # verify how?
            conn_dict_rec = {
                'rule': 'pairwise_bernoulli', 'p': bernoulli_prob
            }
            total_conn_nr = int(len_source*len_target*bernoulli_prob)
        else:
            conn_dict_rec = {
                'rule': 'fixed_total_number', 'N': total_conn_nr
            }
        nest.Connect(
            source_pop, target_pop,
            conn_spec=conn_dict_rec,
            syn_spec=syn_dict
        )
    return total_conn_nr

# assign ctsp
def assign_ctsp(pop, net_dict):
    E_L = net_dict['neuron_params']['E_L']['default']
    V_th = net_dict['neuron_params']['V_th']['default']
    C_m = net_dict['neuron_params']['C_m']['default']
    tau_m = net_dict['neuron_params']['tau_m']['default']
    V_reset = net_dict['neuron_params']['V_reset']['default']
    if net_dict['ctsp'] is True:
        for celltype in ['Exc', 'PV', 'SOM', 'VIP']:
            if celltype in pop:
                E_L = net_dict['neuron_params']['E_L'][celltype]
                V_th = net_dict['neuron_params']['V_th'][celltype]
                C_m = net_dict['neuron_params']['C_m'][celltype]
                tau_m = net_dict['neuron_params']['tau_m'][celltype]
                V_reset = net_dict['neuron_params']['V_reset'][celltype]
                break
    # print('pop={}, E_L={}, V_th={}, C_m={}, tau_m={}'.format(pop, E_L, V_th, C_m, tau_m))
    return E_L, V_th, C_m, tau_m, V_reset


# calculate celltype-specific psc
def calc_psc(psp_val, C_m, tau_m, tau_syn):
    PSC_e_over_PSP_e = (((C_m) ** (-1) * tau_m * tau_syn / (
        tau_syn - tau_m) * ((tau_m / tau_syn) ** (
            - tau_m / (tau_m - tau_syn)) - (tau_m / tau_syn) ** (
                - tau_syn / (tau_m - tau_syn)))) ** (-1))
    PSC_e = (PSC_e_over_PSP_e * psp_val)
    verify_collect('calc_psc(): psp_val={:.2f}, C_m={:.2f}, tau_m={:.2f}, tau_syn={:.2f}, PSC_e={:.2f}\n'.format(psp_val, C_m, tau_m, tau_syn, PSC_e), 'psc')
    return PSC_e


# get psc from single psp
def get_weight(psp_val, net_dict):
    C_m = net_dict['neuron_params']['C_m']['default']
    tau_m = net_dict['neuron_params']['tau_m']['default']
    tau_syn_ex = net_dict['neuron_params']['tau_syn_ex']
    return calc_psc(psp_val, C_m, tau_m, tau_syn_ex)


# get the psc matrix
def get_weight_mtx(net_dict):
    np.set_printoptions(precision=2, linewidth=500, suppress=True)
    lyrs = [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]
    types = ['Exc', 'PV', 'SOM', 'VIP', 'Exc', 'PV', 'SOM', 'Exc', 'PV', 'SOM', 'Exc', 'PV', 'SOM']
    pscs = np.zeros((len(types), len(types)))
    taus = np.zeros(pscs.shape)
    verify_collect('psps=\n{}\n'.format(net_dict['psp_means']), 'psc')
    for i, trg_type in enumerate(types):    # post
        for j, src_type in enumerate(types):    # pre
            psp = net_dict['psp_means'][i, j]
            tau_syn = net_dict['neuron_params']['tau_syn_ex']
            # for inhibitory connections
            if src_type != 'Exc':
                tau_syn = net_dict['neuron_params']['tau_syn_in']
            taus[i, j] = tau_syn
            if net_dict['ctsp'] and net_dict['ctsp_dependent_psc']:
                type = trg_type
            else:
                type = 'default'
            pscs[i, j] = calc_psc(psp,
                net_dict['neuron_params']['C_m'][type],
                net_dict['neuron_params']['tau_m'][type],
                tau_syn)
    verify_collect('pscs=\n{}\n'.format(pscs), 'psc')
    verify_collect('taus=\n{}\n'.format(taus), 'psc')
    return pscs


# equalize inhibitory connectivity; for a network without specific connectivity
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

# connectivity integration
def renew_conn(conn_probs, exp_csv=None):
    if exp_csv is None:
        exp=raw.exp
    else:
        exp=np.loadtxt(exp_csv, delimiter=',')
    print('conn. map before =\n{}'.format(conn_probs))
    conn_probs = conn_barrel_integrate(raw.mouse_dict, raw.bbp, exp, raw.allen, (exp>0)*1)  # flags: use exp. data if > 0
    np.savetxt(exp_csv.replace('raw', 'conn'), conn_probs, fmt='%.4f', delimiter=',')
    print('conn. map renewed =\n{}'.format(conn_probs))
    return conn_probs


'''
from helpers.py
'''
def compute_DC(net_dict, w_ext):
    DC = (
        net_dict['bg_rate'] * net_dict['K_ext'] *
        w_ext * net_dict['neuron_params']['tau_syn_E'] * 0.001
        )
    return DC

# no using (using Bernoulli now)
def get_total_number_of_synapses(net_dict):
    N_full = net_dict['N_full']
    number_N = len(N_full)
    conn_probs = net_dict['conn_probs']
    prod = np.outer(N_full, N_full)
    n_syn_temp = np.log(1. - conn_probs)/np.log((prod - 1.) / prod)
    N_full_matrix = np.column_stack((N_full for i in list(range(number_N))))
    K = (((n_syn_temp * (N_full_matrix).astype(int)) / N_full_matrix).astype(int))
    return K

# no using (using Bernoulli now)
def synapses_th_matrix(net_dict, stim_dict):
    N_full = net_dict['N_full']
    conn_probs = stim_dict['conn_probs_th']
    T_full = stim_dict['n_thal']
    prod = (T_full * N_full).astype(float)
    n_syn_temp = np.log(1. - conn_probs)/np.log((prod - 1.)/prod)
    K = (((n_syn_temp * (N_full).astype(int))/N_full).astype(int))
    return K
