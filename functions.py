import nest
import numpy as np

from scan_params import *

Fmax = False

# short-term plasticity 190605
STP = True
PV_depress = True

# relative strengths 190629
adjust_inh = True
SOM_g_relative = 1.0
PV_g_relative = 1.0

# selectivity 190614
tuning_orientation = False
selective_inh_source = ['PV', 'SOM']
selective_inh_target = ['PV', 'SOM']
k_th = 0.8
k_exc_to_exc = 0.8
k_exc_to_inh = 0.4
k_inh_to_exc = 0.0

# cell-type-specific parameters 190704
ctsp = True

if STP:
    def assign_syn_dict(source_name, target_name, weight_dict, delay_dict, net_dict):
        syn_dict = {
            'model': 'static_synapse',
            'weight': weight_dict,
            'delay': delay_dict
        }
        if 'PC' in source_name:
            if 'PC' in target_name:
                syn_dict = {
                    'model': 'tsodyks_synapse',
                    'U': 0.75,
                    'tau_fac': 0.01,
                    'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
                    'tau_rec': 800.0,
                    'weight': weight_dict,
                    'delay': delay_dict
                }
            elif 'PV' in target_name:
                if PV_depress is True:
                    syn_dict = {
                        'model': 'tsodyks_synapse',
                        'U': 0.75,
                        'tau_fac': 0.01,
                        'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
                        'tau_rec': 800.0,
                        'weight': weight_dict,
                        'delay': delay_dict
                    }
            elif 'SOM' in target_name:
                if som_facilitate is True:
                    tau_fac = 200.0
                    syn_dict = {
                        'model': 'tsodyks_synapse',
                        'U': 0.5,
                        'tau_fac': tau_fac,
                        'tau_psc': net_dict['neuron_params']['tau_syn_ex'],
                        'tau_rec': 0.01,
                        'weight': weight_dict,
                        'delay': delay_dict
                    }
        elif 'PV' in source_name:
            syn_dict = {
                'model': 'tsodyks_synapse',
                'U': 0.9,
                'tau_fac': 0.01,
                'tau_psc': net_dict['neuron_params']['tau_syn_in'],
                'tau_rec': 800.0,
                'weight': weight_dict,
                'delay': delay_dict
            }
        return syn_dict


# not using
def set_fmax(names, population):
    # Fmax:
    # Two Dynamically Distinct Inhibitory Networks in
    # Layer 4 of the Neocortex
    if Fmax:
        if 'PC' in names:
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

if adjust_inh:
    def inh_weight(source_name, weight):
        # HJ: Short-Term Plasticity of Unitary Inhibitory-to-Inhibitory
        # Synapses Depends on the Presynaptic Interneuron Subtype
        # 190629
        if 'SOM' in source_name:
            weight *= SOM_g_relative
        if 'PV' in source_name:
            weight *= PV_g_relative
        return weight


if tuning_orientation is True:
    def set_thalamus_input(th_pop, poisson_pops, start_times, stop_times, rate_0, stim_theta):
        for j, pop in enumerate(poisson_pops):  # for each stimulation; e.g. repeat 20 times
            for i, node in enumerate(pop):
                theta = -np.pi/2.0+np.pi*((i+0.5)/float(len(pop))) # from -pi/2 to pi/2
                rate = rate_0*(1.0+k_th*np.cos(2.0*(theta - stim_theta)))
                # if j == 0:
                #     print(theta, rate)
                nest.SetStatus([node], {
                    'rate': rate,
                    'start': start_times[j],
                    'stop': stop_times[j]
                })
            nest.Connect(pop, th_pop, conn_spec={'rule': 'one_to_one'})


    def connect_thalamus_orientation(th_pop, target_pop, target_name,
                                     nr_synapses, syn_dict_th):
        # do it only if connection is not 0
        if nr_synapses > 0:
            if 'PC' in target_name:
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
                    head_idx = x*len_th_cluster
                    tail_idx = (x+1)*len_th_cluster
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

                # print('pop len = {0}, cluster len = {1}, p_0={2}'
                #       .format(len_target, len_target_cluster, p_0))
                # for cluster in target_cluster_list:
                #     print('cluster len = {0}'.format(len(cluster)))

                for i, th_cluster in enumerate(th_cluster_list):
                    for j, target_cluster in enumerate(target_cluster_list):
                        theta_th = -np.pi/2.0+np.pi*((i+0.5)/float(nr_cluster))
                        theta_target = -np.pi/2.0+np.pi*((j+0.5)/float(nr_cluster))
                        p = p_0*(1.0+k_exc_to_exc*np.cos(2.0*(theta_th - theta_target)))
                        conn_nr = int(round(len(th_cluster)*len(target_cluster)*p))
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
                nest.Connect(
                    th_pop, target_pop,
                    conn_spec={
                        'rule': 'fixed_total_number',
                        'N': nr_synapses,
                    },
                    syn_spec=syn_dict_th
                )


    def connect_by_cluster(source_name, target_name, synapse_nr, syn_dict, source_pop, target_pop):
        nr_cluster = 8
        k = 0.0
        if 'PC' in source_name and 'PC' in target_name:
            k = k_exc_to_exc
        for inh_target in selective_inh_target:
            if 'PC' in source_name and inh_target in target_name:
                k = k_exc_to_inh
        for inh_source in selective_inh_source:
            if inh_source in source_name and 'PC' in target_name:
                k = k_inh_to_exc
        if k != 0.0:
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
                        theta_source = -np.pi / 2.0 + np.pi * ((i+0.5) / float(nr_cluster))
                        theta_target = -np.pi / 2.0 + np.pi * ((j+0.5) / float(nr_cluster))
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
                print(source_name + ' to ' + target_name + ' verify synapse_nr: {0}, {1}'.format(synapse_nr, conn_nr_sum))
        else:
            conn_dict_rec = {
                'rule': 'fixed_total_number', 'N': synapse_nr
            }
            nest.Connect(
                source_pop, target_pop,
                conn_spec=conn_dict_rec,
                syn_spec=syn_dict
            )

if ctsp is True:
    def ctsp_assign(pop, net_dict, E_L, V_th, C_m, tau_m):
        for celltype in ['PC', 'PV', 'SOM', 'VIP']:
            if celltype in pop:
                E_L = net_dict['neuron_params']['E_L'][celltype]
                V_th = net_dict['neuron_params']['V_th'][celltype]
                C_m = net_dict['neuron_params']['C_m'][celltype]
                tau_m = net_dict['neuron_params']['tau_m'][celltype]
                break
        return E_L, V_th, C_m, tau_m


    # def get_weight_ctsp(PSP_val, net_dict, target_name):
    #     """ Computes weight to elicit a change in the membrane potential.
    #
    #     This function computes the weight which elicits a change in the membrane
    #     potential of size PSP_val. To implement this, the weight is calculated to
    #     elicit a current that is high enough to implement the desired change in the
    #     membrane potential.
    #
    #     Parameters
    #     ----------
    #     PSP_val
    #         Evoked postsynaptic potential.
    #     net_dict
    #         Dictionary containing parameters of the microcircuit.
    #
    #     Returns
    #     -------
    #     PSC_e
    #         Weight value(s).
    #
    #     """
    #     cell_type = 'default'
    #     for name in ['PC', 'PV', 'SOM', 'VIP']:
    #         if name in target_name:
    #             cell_type = name
    #             break
    #
    #     C_m = net_dict['neuron_params']['C_m'][cell_type]
    #     tau_m = net_dict['neuron_params']['tau_m'][cell_type]
    #     tau_syn_ex = net_dict['neuron_params']['tau_syn_ex']
    #
    #     PSC_e_over_PSP_e = (((C_m) ** (-1) * tau_m * tau_syn_ex / (
    #             tau_syn_ex - tau_m) * ((tau_m / tau_syn_ex) ** (
    #             - tau_m / (tau_m - tau_syn_ex)) - (tau_m / tau_syn_ex) ** (
    #                                            - tau_syn_ex / (
    #                                                tau_m - tau_syn_ex)))) ** (-1))
    #     PSC_e = (PSC_e_over_PSP_e * PSP_val)
    #     # print(PSC_e)
    #     return PSC_e

