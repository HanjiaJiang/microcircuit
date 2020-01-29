import os
import sys
import numpy as np
import pickle
from microcircuit.network_params import net_dict
from microcircuit.sim_params import sim_dict
from microcircuit.stimulus_params import stim_dict
from microcircuit.functions import special_dict
from stp.stp_dicts import allen_stp, doiron_stp, doiron_stp_weak


def set_constant():
    sim_dict['t_sim'] = 102000.0
    net_dict['g'] = -4
    net_dict['bg_rate'] = 4.0
    net_dict['animal'] = 'mouse'
    net_dict['renew_conn'] = False
    net_dict['neuron_params']['tau_syn_ex'] = 2.1
    net_dict['neuron_params']['tau_syn_inh'] = 3.2
    stim_dict['thalamic_input'] = False
    stim_dict['th_start'] = np.arange(1500.0, sim_dict['t_sim'], 500.0)
    special_dict['orient_tuning'] = False
    special_dict['stp_dict'] = doiron_stp_weak
    net_dict['K_ext'] = np.array([2000, 2000, 1500, 550,
                                  2000, 2000, 1500,
                                  2000, 2000, 1500,
                                  2000, 2000, 1500])
    net_dict['conn_probs'] = \
        np.array([[0.0872, 0.3173, 0.4612, 0.0443, 0.1056, 0.4011, 0.0374, 0.0234, 0.09  , 0.1864, 0.    , 0.    , 0.    ],
           [0.3763, 0.3453, 0.2142, 0.0683, 0.0802, 0.0135, 0.026 , 0.0257, 0.1937, 0.2237, 0.0001, 0.0001, 0.0051],
           [0.2288, 0.4216, 0.0263, 0.2618, 0.0033, 0.0097, 0.0363, 0.0003, 0.0222, 0.018 , 0.    , 0.    , 0.    ],
           [0.0222, 0.0487, 0.0561, 0.027 , 0.0021, 0.0085, 0.0141, 0.0002, 0.0008, 0.0051, 0.    , 0.0001, 0.0047],
           [0.0128, 0.0668, 0.049 , 0.0578, 0.1764, 0.4577, 0.2761, 0.0059, 0.0229, 0.0427, 0.    , 0.0019, 0.0212],
           [0.0329, 0.0132, 0.0188, 0.0438, 0.0937, 0.3487, 0.4068, 0.0078, 0.0228, 0.0389, 0.0011, 0.0024, 0.016 ],
           [0.033 , 0.015 , 0.0198, 0.2618, 0.2906, 0.4432, 0.028 , 0.0087, 0.0263, 0.0384, 0.0016, 0.0019, 0.0198],
           [0.0841, 0.0528, 0.072 , 0.0534, 0.0844, 0.0573, 0.0621, 0.0957, 0.1871, 0.1575, 0.0094, 0.0146, 0.0418],
           [0.0705, 0.1211, 0.0444, 0.0169, 0.0315, 0.025 , 0.0188, 0.0846, 0.3574, 0.2594, 0.0041, 0.0107, 0.0213],
           [0.0998, 0.0072, 0.0089, 0.2618, 0.0343, 0.0248, 0.0209, 0.0587, 0.1182, 0.0373, 0.0054, 0.0122, 0.0262],
           [0.    , 0.0017, 0.0029, 0.007 , 0.0297, 0.0133, 0.0086, 0.0381, 0.0162, 0.0138, 0.021 , 0.3249, 0.3014],
           [0.0026, 0.0001, 0.0002, 0.0019, 0.0047, 0.002 , 0.0004, 0.015 , 0.    , 0.0028, 0.1865, 0.3535, 0.2968],
           [0.0021, 0.    , 0.0002, 0.2618, 0.0043, 0.0018, 0.0003, 0.0141, 0.    , 0.0019, 0.1955, 0.3321, 0.0307]])


def params_single(path):
    set_constant()
    # sim_dict['master_seed'] = 55
    # net_dict['renew_conn'] = True
    # stim_dict['orientation'] = 0.0
    para_dict = {
        'net_dict': net_dict,
        'sim_dict': sim_dict,
        'stim_dict': stim_dict,
        'special_dict': special_dict
    }
    with open(path, 'wb') as h:
        pickle.dump(para_dict, h)


def get_conn_probs(list_n=10):
    cwd = os.getcwd()
    conn_probs_list = []
    conn_folder = os.path.join(cwd, 'conn_probs')
    for file in os.listdir(conn_folder):
        if file.endswith(".npy") and 'conn_probs' in file:
            tmp = np.load(os.path.join(conn_folder, file))
            if tmp.shape == net_dict['conn_probs'].shape:
                conn_probs_list.append(tmp)
    if len(conn_probs_list) < list_n:
        return_list = conn_probs_list
    else:
        return_list = conn_probs_list[:list_n]
    return return_list


def read_levels(in_str):
    out_str = os.path.basename(in_str).split('.')[0]
    out_list = np.array(out_str.split('_')).astype(int)
    return out_str, out_list


def g_bg(out_list):
    set_constant()
    max_str, max_list = read_levels(out_list[-1])
    # net_dict['K_ext'] = np.array([2000, 2000, 1600, 600,
    #                               2000, 2000, 1600,
    #                               2000, 2000, 1600,
    #                               2000, 2000, 1600])
    stp_list = [doiron_stp_weak, allen_stp]
    for i, output in enumerate(out_list):
        levels_str, levels_list = read_levels(output)
        special_dict['stp_dict'] = stp_list[levels_list[0]]
        net_dict['g'] = -4.0 + -4.0 * (float(levels_list[1]) / (max_list[1]))
        net_dict['bg_rate'] = 4.0 + 4.0*(float(levels_list[2]) / (max_list[2]))
        sim_dict['data_path'] = os.path.join(os.path.dirname(output), levels_str)
        para_dict = {
            'net_dict': net_dict,
            'sim_dict': sim_dict,
            'stim_dict': stim_dict,
            'special_dict': special_dict
        }
        with open(output, 'wb') as handle:
            pickle.dump(para_dict, handle)
        handle.close()

def conn_stp_som_vip(out_list):
    # set constant parameters
    set_constant()

    # get the size of each dimension from the largest (#_#_#_#)
    max_str, max_list = read_levels(out_list[-1])
    print(max_list)

    # load conn_probs; default return the whole list
    conn_probs_list = get_conn_probs()
    for conn_probs in conn_probs_list:
        print(conn_probs)

    # stp list
    stp_list = [allen_stp, doiron_stp, doiron_stp_weak]

    for i, output in enumerate(out_list):
        # get levels
        levels_str = os.path.basename(output).split('.')[0]
        levels_list = np.array(levels_str.split('_')).astype(int)

        # assign to dictionary
        net_dict['conn_probs'] = conn_probs_list[levels_list[0]]
        # print(levels_list[1])
        special_dict['stp_dict'] = stp_list[levels_list[1]]
        som = (float(levels_list[2] + 1) / (
                    max_list[2] + 1)) * 2000.0  # assign som and vip strengths
        vip = (float(levels_list[3] + 1) / (max_list[3] + 1)) * 2000.0  # according to levels
        sim_dict['data_path'] = os.path.join(os.path.dirname(output), levels_str)
        net_dict['K_ext'] = np.array([2000, 2000, som, vip,
                                      2000, 2000, som,
                                      2000, 2000, som,
                                      2000, 2000, som])
        para_dict = {
            'net_dict': net_dict,
            'sim_dict': sim_dict,
            'stim_dict': stim_dict,
            'special_dict': special_dict
        }

        with open(output, 'wb') as handle:
            pickle.dump(para_dict, handle)
        handle.close()


if __name__ == "__main__":
    np.set_printoptions(precision=4, suppress=True)

    # get output names from system input
    output_list = sys.argv[1:]

    g_bg(output_list)
    # conn_stp_som_vip(out_list)