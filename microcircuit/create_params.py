import os
import sys
import copy
import json
import pickle
import numpy as np
from microcircuit.network_params import net_dict
from microcircuit.sim_params import sim_dict
from microcircuit.stimulus_params import stim_dict
from microcircuit.functions import special_dict
import microcircuit.functions as func
import microcircuit.stp.stp_dicts as stps
np.set_printoptions(suppress=True, precision=4)

'''
objects
'''
stp = stps.fitted_stp
class ScanParams:
    def __init__(self, n_dict, si_dict, st_dict, sp_dict):
        self.net_dict = copy.deepcopy(n_dict)
        self.sim_dict = copy.deepcopy(si_dict)
        self.stim_dict = copy.deepcopy(st_dict)
        self.special_dict = copy.deepcopy(sp_dict)

'''
preliminary settings
'''
# set layer-specific thalamic input
def set_thalamic(para_dict, th_starts=None, th_rate=None, orient=False, duration=10):
    th_dict = {}
    if type(th_starts) is list and len(th_starts) > 0 and type(th_rate) is float:
        para_dict['stim_dict']['thalamic_input'] = True
        para_dict['stim_dict']['th_rate'] = th_rate
        para_dict['stim_dict']['th_start'] = np.array(th_starts).astype(float)
        para_dict['stim_dict']['th_duration'] = duration
        # Bruno, Simons, 2002; Oberlaender et al., 2011; Sermet et al., 2019; Constantinople, Bruno, 2013
        para_dict['stim_dict']['conn_probs_th'] = np.array([0.058, 0.058, 0.0, 0.0, 0.4, 0.4, 0.0, 0.259, 0.259, 0.0, 0.09, 0.09, 0.0])
    para_dict['special_dict']['orient_tuning'] = orient

# set constant parameters
def set_constant():
    # testing area
    # special_dict['ctsp'] = False
    # net_dict['ipsp']['use'] = True
    # net_dict['renew_conn'] = True
    # sim_dict['master_seed'] = 60

    net_dict['g'] = -8
    net_dict['bg_rate'] = 4.0
    special_dict['stp_dict'] = stp
    exc, pv, som, vip = 1000, 2000, 1000, 1000
    net_dict['K_ext'] = np.array([exc, pv, som, vip, exc, pv, som, exc, pv, som, exc, pv, som])

    # 20-04-26
    # L6 conn by averaging
    # 20-05-26 VIP-SOM conn adjusted
    net_dict['conn_probs'] = \
        np.array([ [0.0839, 0.3053, 0.4438, 0.0522, 0.1039, 0.0192, 0.0429, 0.0232, 0.0891, 0.1845, 0.    , 0.    , 0.    ],
                   [0.3621, 0.3323, 0.2061, 0.0806, 0.0042, 0.0155, 0.0298, 0.0254, 0.1918, 0.2215, 0.0001, 0.0001, 0.0054],
                   [0.2201, 0.4057, 0.0254, 0.2519, 0.0038, 0.0111, 0.0417, 0.0004, 0.022 , 0.0199, 0.    , 0.    , 0.    ],
                   [0.0262, 0.0574, 0.0662, 0.0318, 0.0024, 0.0097, 0.0162, 0.0002, 0.0009, 0.0056, 0.    , 0.0001, 0.005 ],

                   [0.0126, 0.0333, 0.0562, 0.0663, 0.1668, 0.4327, 0.261 , 0.0058, 0.0264, 0.0491, 0.    , 0.0021, 0.0232],
                   [0.0378, 0.0152, 0.0216, 0.0503, 0.0886, 0.3297, 0.3846, 0.009 , 0.0262, 0.0446, 0.0013, 0.0026, 0.0175],
                   [0.0379, 0.0172, 0.0227, 0.2519, 0.0859, 0.419 , 0.0264, 0.01  , 0.0303, 0.0441, 0.0017, 0.0021, 0.0217],

                   [0.0832, 0.0523, 0.0713, 0.0589, 0.0826, 0.0658, 0.0714, 0.091 , 0.178 , 0.1498, 0.0093, 0.0167, 0.0477],
                   [0.0698, 0.1199, 0.0439, 0.0186, 0.0362, 0.0288, 0.0216, 0.0804, 0.34  , 0.2468, 0.0047, 0.0122, 0.0243],
                   [0.0988, 0.0071, 0.0098, 0.2519, 0.0394, 0.0285, 0.024 , 0.0558, 0.1124, 0.0355, 0.0061, 0.014 , 0.0299],

                   [0.    , 0.0018, 0.0031, 0.0075, 0.0291, 0.0145, 0.0094, 0.0374, 0.0184, 0.0157, 0.0199, 0.3083, 0.286 ],
                   [0.0028, 0.0001, 0.0002, 0.002 , 0.0052, 0.0022, 0.0005, 0.0171, 0.    , 0.0032, 0.177 , 0.3355, 0.2817],
                   [0.0022, 0.    , 0.0002, 0.2519, 0.0047, 0.0019, 0.0003, 0.0161, 0.    , 0.0021, 0.1208, 0.3151, 0.0292]])

    # net_dict['conn_probs'] = func.eq_inh_conn(net_dict['N_full'], net_dict['conn_probs'])

'''
tools
'''
def print_all(all_dict):
    os.system('mkdir -p ' + all_dict['sim_dict']['data_path'])
    for k1, v1 in all_dict.items():
        with open(os.path.join(all_dict['sim_dict']['data_path'], 'params_{}.txt'.format(k1)), 'w') as f:
            for k2, v2 in v1.items():
                if isinstance(v2, dict):
                    f.write(k2 + ':\n')
                    for k3, v3 in v2.items():
                        f.write(k3 + '=\n')
                        f.write('{}\n'.format(v3))
                        # if isinstance(v3, np.ndarray):
                        #     f.write('hash={}\n'.format(hash(bytes(v3))))
                    f.write('\n')
                else:
                    f.write(k2 + ':\n')
                    f.write('{}\n'.format(v2))
                    # if isinstance(v2, np.ndarray):
                    #     f.write('hash={}\n'.format(hash(bytes(v2))))
                    f.write('\n')
            f.close()

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

'''
scan setting
'''
# multi-parameter
def set_g_bg(all_dict, lvls):
    all_dict['net_dict']['g'] = -float(lvls[0])
    all_dict['net_dict']['bg_rate'] = float(lvls[1])

def set_indg_somvip(all_dict, lvls):
    all_dict['net_dict']['K_ext'] = np.array([2000, 2000, lvls[0], lvls[1],
                                              2000, 2000, lvls[0],
                                              2000, 2000, lvls[0],
                                              2000, 2000, lvls[0]])

def set_indg_all(all_dict, lvls):
    all_dict['net_dict']['K_ext'] = np.array([lvls[0], lvls[1], lvls[2], lvls[3],
                                              lvls[0], lvls[1], lvls[2],
                                              lvls[0], lvls[1], lvls[2],
                                              lvls[0], lvls[1], lvls[2]])


# single-parameter
def set_stp_config(all_dict, lvl):
    stp_list = [copy.deepcopy(stp)]
    for key_a in stp.keys():
        for key_b in stp[key_a].keys():
            tmp_stp = copy.deepcopy(stp)
            tmp_stp[key_a][key_b] = {'model': 'static_synapse'}
            stp_list.append(tmp_stp)
    if lvl < len(stp_list):
        all_dict['special_dict']['stp_dict'] = stp_list[lvl]

def set_ctsp(all_dict, lvl):
    if lvl == 0:
        all_dict['special_dict']['ctsp'] = False
    else:
        all_dict['special_dict']['ctsp'] = True

def set_seed(all_dict, lvl):
    all_dict['sim_dict']['master_seed'] = 55 + lvl

def set_properties(all_dict, lvl):
    if lvl < 8:
        # set equal conn for INs
        if lvl < 4:
            all_dict['net_dict']['conn_probs'] = func.eq_inh_conn(net_dict['N_full'], net_dict['conn_probs'])
        else:
            all_dict['net_dict']['conn_probs'] = copy.deepcopy(net_dict['conn_probs'])
        ctsp = [False, False, True, True, False, False, True, True]
        stps = [{}, stp, {}, stp, {}, stp, {}, stp]
        all_dict['special_dict']['ctsp'] = ctsp[lvl]
        all_dict['special_dict']['stp_dict'] = stps[lvl]

def set_indg_vip(all_dict, lvl):
    all_dict['net_dict']['K_ext'][3] = lvl


'''
main creator functions
'''
def set_main(out_list, scanparams, f1, f2, f3=None, f4=None):
    all_dict = {
        'net_dict': scanparams.net_dict,
        'sim_dict': scanparams.sim_dict,
        'stim_dict': scanparams.stim_dict,
        'special_dict': scanparams.special_dict
    }
    for i, out in enumerate(out_list):
        lvls_list = read_levels(out)[1]
        print('lvls_list={}'.format(lvls_list))
        f1(all_dict, lvls_list[:])
        # f2(all_dict, lvls_list[1])
        # if f3 is not None:
        #     f3(all_dict, lvls_list[2:])
        save_pickle(out, all_dict)


# set parameters for single-run
def set_single(pickle_path):
    set_constant()
    scanparams = ScanParams(net_dict, sim_dict, stim_dict, special_dict)
    all_dict = {
        'net_dict': scanparams.net_dict,
        'sim_dict': scanparams.sim_dict,
        'stim_dict': scanparams.stim_dict,
        'special_dict': scanparams.special_dict
    }
    with open(pickle_path, 'wb') as h:
        pickle.dump(all_dict, h)


if __name__ == "__main__":
    # get output names from system input
    output_list = sys.argv[1:]

    # set constant parameters
    set_constant()

    # create ScanParams object
    scanparams = ScanParams(net_dict, sim_dict, stim_dict, special_dict)

    # set the network with parameters
    # set_main(output_list, scanparams, set_properties, set_indg_vip, set_g_bg)
    set_main(output_list, scanparams, set_indg_all, f2=None)
    # set_main(output_list, scanparams, set_indg_somvip, set_g_bg)
    # set_main(output_list, scanparams, set_stp_config, set_seed, set_g_bg)
