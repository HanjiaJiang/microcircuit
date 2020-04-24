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
import copy
np.set_printoptions(suppress=True, precision=4)


# set layer-specific thalamic input
def set_thalamic(para_dict, th_starts=None, th_rate=None, orient=False, duration=10):
    th_dict = {}
    if type(th_starts) is list and len(th_starts) > 0 and type(th_rate) is float:
        para_dict['stim_dict']['thalamic_input'] = True
        para_dict['stim_dict']['th_rate'] = th_rate
        para_dict['stim_dict']['th_start'] = np.array(th_starts).astype(float)
        para_dict['stim_dict']['th_duration'] = duration
        para_dict['stim_dict']['n_thal'] = 200 # Constantinople, Bruno, 2013; Oberlaender et al., 2011
        # Bruno, Simons, 2002; Oberlaender et al., 2011; Sermet et al., 2019; Constantinople, Bruno, 2013
        para_dict['stim_dict']['conn_probs_th'] = np.array([0.058, 0.058, 0.0, 0.0, 0.4, 0.4, 0.0, 0.259, 0.259, 0.0, 0.09, 0.09, 0.0])
        para_dict['special_dict']['orient_tuning'] = orient
    print_summary(para_dict)

# set constant parameters
def set_constant():
    # testing area
    net_dict['neuron_params']['tau_syn_ex'] = 2.37  # Feldmeyer, Sakmann, 2002, J of Physiology
    net_dict['neuron_params']['tau_syn_in'] = 8.3   # Gupta, Markram, 2000, Science
    # net_dict['renew_conn'] = True

    net_dict['g'] = -6
    net_dict['bg_rate'] = 4.0
    special_dict['stp_dict'] = doiron_stp_weak
    net_dict['K_ext'] = np.array([2000, 1800, 800, 300,
                                  2000, 1800, 800,
                                  2000, 2000, 500,
                                  2000, 2000, 500])
    # 20-04-21
    # 1st: L6 conn by averaging; 2nd: L6 conn as L5
    net_dict['conn_probs'] = \
    np.array([ [0.1014, 0.3688, 0.536 , 0.0522, 0.1174, 0.0192, 0.0429, 0.0247, 0.0948, 0.1963, 0.    , 0.    , 0.    ],
               [0.4374, 0.4014, 0.249 , 0.0806, 0.0042, 0.0155, 0.0298, 0.027 , 0.2039, 0.2355, 0.0001, 0.0001, 0.0054],
               [0.2659, 0.4901, 0.0306, 0.3043, 0.0038, 0.0111, 0.0417, 0.0004, 0.0234, 0.0199, 0.    , 0.    , 0.    ],
               [0.0262, 0.0574, 0.0662, 0.0318, 0.0024, 0.0097, 0.0162, 0.0002, 0.0009, 0.0056, 0.    , 0.0001, 0.005 ],
               [0.0142, 0.0333, 0.0562, 0.0663, 0.2058, 0.5339, 0.322 , 0.0065, 0.0264, 0.0491, 0.    , 0.0021, 0.0232],
               [0.0378, 0.0152, 0.0216, 0.0503, 0.1093, 0.4068, 0.4746, 0.009 , 0.0262, 0.0446, 0.0013, 0.0026, 0.0175],
               [0.0379, 0.0172, 0.0227, 0.0458, 0.1059, 0.5169, 0.0326, 0.01  , 0.0303, 0.0441, 0.0017, 0.0021, 0.0217],
               [0.0885, 0.0556, 0.0758, 0.0589, 0.0932, 0.0658, 0.0714, 0.1098, 0.2146, 0.1807, 0.0103, 0.0167, 0.0477],
               [0.0742, 0.1275, 0.0467, 0.0186, 0.0362, 0.0288, 0.0216, 0.097 , 0.41  , 0.2976, 0.0047, 0.0122, 0.0243],
               [0.1051, 0.0076, 0.0098, 0.0184, 0.0394, 0.0285, 0.024 , 0.0673, 0.1355, 0.0428, 0.0061, 0.014 , 0.0299],
               [0.    , 0.0018, 0.0031, 0.0075, 0.0306, 0.0145, 0.0094, 0.0418, 0.0184, 0.0157, 0.0241, 0.3737, 0.3468],
               [0.0028, 0.0001, 0.0002, 0.002 , 0.0052, 0.0022, 0.0005, 0.0171, 0.    , 0.0032, 0.2145, 0.4067, 0.3414],
               [0.0022, 0.    , 0.0002, 0.0011, 0.0047, 0.0019, 0.0003, 0.0161, 0.    , 0.0021, 0.1464, 0.382 , 0.0353]])
    # net_dict['conn_probs'] = np.load('conn_probs/conn_probs_200423.npy')
    # net_dict['conn_probs'] = np.load('conn_probs/conn_probs_l6-as-l5.npy')


# print summary of parameters
def print_summary(all_dict):
    if not os.path.isdir(all_dict['sim_dict']['data_path']):
        os.mkdir(all_dict['sim_dict']['data_path'])
    with open(os.path.join(all_dict['sim_dict']['data_path'], 'params_summary.txt'), 'w') as f:
        f.write('master seed = {}\n'.format(all_dict['sim_dict']['master_seed']))
        f.write('ctsp = {}\n\n'.format(all_dict['special_dict']['ctsp']))
        f.write('stp = \n{}\n\n'.format(all_dict['special_dict']['stp_dict']))
        f.write('thalamic input = {}\n\n'.format(all_dict['stim_dict']['thalamic_input']))
        if all_dict['stim_dict']['thalamic_input'] is True:
            f.write('    th_rate = {}, th_start = {}\n\n'.format(all_dict['stim_dict']['th_rate'], all_dict['stim_dict']['th_start']))
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


# read parameters from given input string
def read_levels(in_str):
    out_str = os.path.basename(in_str).split('.')[0]    # filename with .
    out_list = np.array(out_str.split('_')).astype(int)
    return out_str, out_list


# list generation
# connection probability map
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


# psp list
def list_psp():
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

def set_ins(all_dict, lvls):
    all_dict['net_dict']['K_ext'] = np.array([2000, 2000, lvls[0], lvls[1],
                                              2000, 2000, lvls[0],
                                              2000, 2000, lvls[0],
                                              2000, 2000, lvls[0]])

def set_K_ext(all_dict, lvls):
    all_dict['net_dict']['K_ext'] = np.array([lvls[0], lvls[1], lvls[2], lvls[3],
                                              lvls[0], lvls[1], lvls[2],
                                              lvls[0], lvls[1], lvls[2],
                                              lvls[0], lvls[1], lvls[2]])


# single-parameter
def set_stp_config(all_dict, lvl):
    stp_list = [copy.deepcopy(doiron_stp_weak)]
    for key_a in doiron_stp_weak.keys():
        for key_b in doiron_stp_weak[key_a].keys():
            tmp_stp = copy.deepcopy(doiron_stp_weak)
            tmp_stp[key_a][key_b] = {'model': 'static_synapse'}
            stp_list.append(tmp_stp)
    if lvl < len(stp_list):
        all_dict['special_dict']['stp_dict'] = stp_list[lvl]

def set_stp(all_dict, lvl):
    stp_list = [no_stp, doiron_stp_weak, allen_stp]
    all_dict['special_dict']['stp_dict'] = stp_list[lvl]

def set_ctsp(all_dict, lvl):
    if lvl == 0:
        all_dict['special_dict']['ctsp'] = False
    else:
        all_dict['special_dict']['ctsp'] = True

def set_seed(all_dict, lvl):
    all_dict['sim_dict']['master_seed'] = 55 + lvl

def set_ctsp_stp(all_dict, lvl):
    if lvl < 4:
        ctsp = [False, False, True, True]
        stp = [no_stp, doiron_stp_weak, no_stp, doiron_stp_weak]
        all_dict['special_dict']['ctsp'] = ctsp[lvl]
        all_dict['special_dict']['stp_dict'] = stp[lvl]

def set_conn_ctsp_stp(all_dict, lvl):
    if lvl < 8:
        if lvl < 4:
            all_dict['net_dict']['conn_probs'] = func.eq_inh_conn(all_dict['net_dict']['N_full'], all_dict['net_dict']['conn_probs'])
        ctsp = [False, False, True, True, False, False, True, True]
        stp = [no_stp, doiron_stp_weak, no_stp, doiron_stp_weak, no_stp, doiron_stp_weak, no_stp, doiron_stp_weak]
        all_dict['special_dict']['ctsp'] = ctsp[lvl]
        all_dict['special_dict']['stp_dict'] = stp[lvl]

def set_vip(all_dict, lvl):
    all_dict['net_dict']['K_ext'][3] = lvl


# main loop
def set_main(out_list, f1, f2, f3=None, f4=None):
    all_dict = {
        'net_dict': net_dict,
        'sim_dict': sim_dict,
        'stim_dict': stim_dict,
        'special_dict': special_dict
    }
    for i, out in enumerate(out_list):
        lvls_list = read_levels(out)[1]
        f1(all_dict, lvls_list[0])
        f2(all_dict, lvls_list[1])
        if f3 is not None:
            f3(all_dict, lvls_list[2:])
        save_pickle(out, all_dict)


# set parameters for single-run
def params_single(path):
    set_constant()
    # net_dict['conn_probs'] = func.eq_inh_conn(net_dict['N_full'], net_dict['conn_probs'])
    # special_dict['stp_dict'] = no_stp
    # special_dict['ctsp'] = False
    all_dict = {
        'net_dict': net_dict,
        'sim_dict': sim_dict,
        'stim_dict': stim_dict,
        'special_dict': special_dict
    }
    with open(path, 'wb') as h:
        pickle.dump(all_dict, h)


if __name__ == "__main__":
    # get output names from system input
    output_list = sys.argv[1:]

    # set constant parameters
    set_constant()

    # set the network with parameters
    set_main(output_list, set_conn_ctsp_stp, set_vip, set_g_bg)
    # set_main(output_list, set_K_ext, f2=None)
    # set_main(output_list, set_ins, set_g_bg)
    # set_main(output_list, set_stp_config, set_seed, set_g_bg)
