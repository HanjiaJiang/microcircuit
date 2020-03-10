import os
import sys
import numpy as np
import pickle
from microcircuit.network_params import net_dict
from microcircuit.sim_params import sim_dict
from microcircuit.stimulus_params import stim_dict
from microcircuit.functions import special_dict
import microcircuit.functions as func
from stp.stp_dicts import no_stp, allen_stp, doiron_stp, doiron_stp_weak
np.set_printoptions(suppress=True, precision=4)

# set constant parameters
def set_constant(th_starts=None, th_rate=None):
    net_dict['g'] = -4
    net_dict['bg_rate'] = 8.0
    net_dict['animal'] = 'mouse'
    # net_dict['conn_probs'] = funcs.eq_inh_conn(net_dict['N_full'], net_dict['conn_probs'])
    special_dict['orient_tuning'] = False
    special_dict['stp_dict'] = doiron_stp_weak
    net_dict['K_ext'] = np.array([2000, 2000, 1500, 600,
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

    # thalamic input
    if th_starts is not None and isinstance(th_starts, np.ndarray) \
            and th_rate is not None:
        stim_dict['thalamic_input'] = True
        stim_dict['th_rate'] = th_rate
        stim_dict['th_start'] = th_starts


# print summary of parameters
def print_summary(all_dict):
    if not os.path.isdir(all_dict['sim_dict']['data_path']):
        os.mkdir(all_dict['sim_dict']['data_path'])
    with open(os.path.join(all_dict['sim_dict']['data_path'], 'params_summary.txt'), 'w') as f:
        f.write('master seed = {}'.format(all_dict['sim_dict']['master_seed']))
        f.write('ctsp = {}\n\n'.format(all_dict['special_dict']['ctsp']))
        f.write('stp = \n{}\n\n'.format(all_dict['special_dict']['stp_dict']))
        f.write('thalamic input = {}\n\n'.format(all_dict['stim_dict']['thalamic_input']))
        f.write('orient_tuning = {}\n\n'.format(all_dict['special_dict']['orient_tuning']))
        f.write('g = {}\n\n'.format(all_dict['net_dict']['g']))
        f.write('bg_rate = {}\n\n'.format(all_dict['net_dict']['bg_rate']))
        f.write('K_ext = \n{}\n\n'.format(all_dict['net_dict']['K_ext']))
        f.write('w_dict = \n{}\n\n'.format(all_dict['net_dict']['w_dict']))
        f.write('conn_probs = \n{}\n\n'.format(all_dict['net_dict']['conn_probs']))


# save assigned parameters to pickle
def save_pickle(pickle_str, all_dict):
    lvls_str = os.path.basename(pickle_str).split('.')[0]
    all_dict['sim_dict']['data_path'] = os.path.join(os.path.dirname(pickle_str), lvls_str)
    with open(pickle_str, 'wb') as handle:
        pickle.dump(all_dict, handle)
    handle.close()


# read parameters from given input string
def read_levels(in_str):
    out_str = os.path.basename(in_str).split('.')[0]    # filename with .
    out_list = np.array(out_str.split('_')).astype(int)
    return out_str, out_list


# set parameters for single-run
def params_single(path):
    set_constant(th_starts=np.arange(3000.0, 8000.0, 2000.0), th_rate=240.0)

    # properties
    # net_dict['conn_probs'] = funcs.eq_inh_conn(net_dict['N_full'], net_dict['conn_probs'])
    # special_dict['stp_dict'] = no_stp
    # special_dict['ctsp'] = False

    para_dict = {
        'net_dict': net_dict,
        'sim_dict': sim_dict,
        'stim_dict': stim_dict,
        'special_dict': special_dict
    }
    with open(path, 'wb') as h:
        pickle.dump(para_dict, h)

    print_summary(para_dict)


# get different connection probability map
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


# set psps
def test_list_psp():
    w_dict_normal = {
        'psp_mtx':
            np.full((4, 4), 0.5),
        'psp_std_mtx':
            np.full((4, 4), 1.0)}
    w_dict_specific = {
        'psp_mtx':
        np.array([[0.70, 0.78, 0.47, 0.23],
                  [0.34, 0.95, 0.38, 0.23],
                  [0.70, 0.63, 0.68, 0.23],
                  [0.70, 2.27, 0.40, 0.53]]),
        'psp_std_mtx':
        np.array([[0.8958, 1.2372, 0.7228, 1.0000],
                  [0.4540, 1.3421, 1.0000, 1.0000],
                  [1.0520, 0.9618, 1.2379, 1.0000],
                  [1.0520, 1.3124, 0.8739, 1.3884]])}
    return [w_dict_normal, w_dict_specific]


# double-parameter
def set_g_bg(all_dict, lvls):
    all_dict['net_dict']['g'] = -float(lvls[0])
    all_dict['net_dict']['bg_rate'] = float(lvls[1])
    return all_dict

def set_ins(all_dict, lvls):
    all_dict['net_dict']['K_ext'] = np.array([2000, 2000, lvls[0], lvls[1],
                                              2000, 2000, lvls[0],
                                              2000, 2000, lvls[0],
                                              2000, 2000, lvls[0]])


# single-parameter
def set_stp(all_dict, lvl):
    stp_list = [no_stp, doiron_stp_weak, allen_stp]
    all_dict['special_dict']['stp_dict'] = stp_list[lvl]

def set_ctsp(all_dict, lvl):
    if lvl == 0:
        all_dict['special_dict']['ctsp'] = False
    else:
        all_dict['special_dict']['ctsp'] = True


# main loop
def set_main(out_list, f1, f2, f_double):
    origin_seed = sim_dict['master_seed']
    all_dict = {
        'net_dict': net_dict,
        'sim_dict': sim_dict,
        'stim_dict': stim_dict,
        'special_dict': special_dict
    }
    for i, out in enumerate(out_list):
        all_dict['sim_dict']['master_seed'] = origin_seed + i
        lvls_list = read_levels(out)[1]
        f1(all_dict, lvls_list[0])
        f2(all_dict, lvls_list[1])
        f_double(all_dict, lvls_list[2:])
        save_pickle(out, all_dict)
        print_summary(all_dict)


if __name__ == "__main__":
    # get output names from system input
    output_list = sys.argv[1:]

    # set constant parameters
    set_constant()

    # set the network with parameters
    set_main(output_list, set_ctsp, set_stp, set_ins)
