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
from microcircuit.stp.stp_dicts import stps
np.set_printoptions(suppress=True, precision=4)

'''
objects
'''
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
def set_constant(constants=None):
    if constants is None:
        constants = ['6-6', '2', '1', '0']  # conn, stp, layer-specific epsp and ipsp

    net_dict['g'] = -8
    net_dict['bg_rate'] = 4.0

    # psps
    if int(constants[2]) == 1:
        if int(constants[3]) == 1:
            net_dict['ipsp']['use'] = True
    else:
        net_dict['epsp']['means'] = np.full((4, 4), 0.5)
        net_dict['epsp']['stds'] = np.full((4, 4), 1.0)

    # stp
    if 0 < int(constants[1]) < len(stps):
        stp_name = list(stps)[int(constants[1])]
        special_dict['stp_dict'] = stps[stp_name]
        print('stp used = {}'.format(stp_name))
    else:
        print('stp not found; using static synapse')

    # conn
    exc, pv, som, vip = 1000, 2000, 1000, 1000
    net_dict['K_ext'] = np.array([exc, pv, som, vip, exc, pv, som, exc, pv, som, exc, pv, som])

    # connectivity map
    # net_dict['conn_probs'] = func.renew_conn(net_dict['conn_probs'], 'microcircuit/conn_probs/raw_2020-5.csv')
    net_dict['conn_probs'] = np.loadtxt('microcircuit/conn_probs/conn_2020-{}.csv'.format(constants[0]), delimiter=',')
    # net_dict['conn_probs'] = np.loadtxt('microcircuit/conn_probs/conn_bbp.csv', delimiter=',')
    adjust_vip_conn(net_dict['conn_probs'])
    # net_dict['conn_probs'] = func.eq_inh_conn(net_dict['N_full'], net_dict['conn_probs'])

'''
tools
'''
def adjust_vip_conn(conn_probs):
    if isinstance(conn_probs, np.ndarray):
        # vip-to-som change to layer-no-specific
        conn_probs[[6, 9, 12], 3] = conn_probs[2, 3]

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
    lvls_str = os.path.basename(pickle_str).split('.')[0] + '/'
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
    stp_names = list(all_dict['special_dict']['stp_dict'])
    if lvl < len(stp_names):
        for k, v in all_dict['special_dict']['stp_dict'][stp_names[lvl]].items():
            all_dict['special_dict']['stp_dict'][stp_names[lvl]][k] = {'model': 'static_synapse'}

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
            all_dict['net_dict']['conn_probs'] = func.eq_inh_conn(all_dict['net_dict']['N_full'], all_dict['net_dict']['conn_probs'])
        # else:
        #     all_dict['net_dict']['conn_probs'] = copy.deepcopy(net_dict['conn_probs'])
        ctsp = [False, False, True, True, False, False, True, True]
        all_dict['special_dict']['ctsp'] = ctsp[lvl]
        if lvl % 2 == 0:
            all_dict['special_dict']['stp_dict'] = {}
    else:
        print('set_properties(): lvl too large')

def set_indg_vip(all_dict, lvl):
    all_dict['net_dict']['K_ext'][3] = lvl


'''
main creator functions
'''
def set_main(out_list, scanparams, f1, f2=None, f3=None, f4=None):
    all_dict = {
        'net_dict': scanparams.net_dict,
        'sim_dict': scanparams.sim_dict,
        'stim_dict': scanparams.stim_dict,
        'special_dict': scanparams.special_dict
    }
    for i, out in enumerate(out_list):
        lvls_list = read_levels(out)[1]
        # print('lvls_list={}'.format(lvls_list))
        f1(all_dict, lvls_list[:])
        if f2 is not None:
            f2(all_dict, lvls_list[1])
        if f3 is not None:
            f3(all_dict, lvls_list[2:])
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
    output_list = sys.argv[5:]
    constant_list = sys.argv[1:5]

    # set constant parameters
    set_constant(constant_list)

    # create ScanParams object
    scanparams = ScanParams(net_dict, sim_dict, stim_dict, special_dict)

    # set the network with parameters
    # set_main(output_list, scanparams, set_properties, set_indg_vip, set_g_bg)
    set_main(output_list, scanparams, set_indg_all)
    # set_main(output_list, scanparams, set_indg_somvip, set_g_bg)
    # set_main(output_list, scanparams, set_stp_config, set_seed, set_g_bg)
