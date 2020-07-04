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
    def __init__(self):
        self.net_dict = copy.deepcopy(net_dict)
        self.sim_dict = copy.deepcopy(sim_dict)
        self.stim_dict = copy.deepcopy(stim_dict)
        self.special_dict = copy.deepcopy(special_dict)
        self.stps = copy.deepcopy(stps)

    def save_pickle(self, pickle_path):
        all_dict = {
            'net_dict': self.net_dict,
            'sim_dict': self.sim_dict,
            'stim_dict': self.stim_dict,
            'special_dict': self.special_dict
        }
        with open(pickle_path, 'wb') as h:
            pickle.dump(all_dict, h)

    def set_path(self, pickle_path, lvls_str):
        self.sim_dict['data_path'] = os.path.join(os.path.dirname(pickle_path), lvls_str)

    def do_single(self, pickle_path):
        self.set_constant()
        self.save_pickle(pickle_path)

    def set_constant(self):
        self.net_dict['g'] = -8
        self.net_dict['bg_rate'] = 4.0
        self.net_dict['epsp']['use'] = True
        self.net_dict['ipsp']['use'] = False
        self.special_dict['stp_dict'] = self.stps['stp_fitted_02.pickle']
        # self.special_dict['stp_dict'] = self.stps['stp_fitted_02.pickle']
        self.set_indgs([1000, 2000, 800, 1000])
        self.load_conn('6-6')

    def set_g(self, g):
        self.net_dict['g'] = -int(np.abs(g))

    def set_bg(self, bg):
        self.net_dict['bg_rate'] = float(bg)

    def set_epsp(self, lyr_specific):
        if int(lyr_specific) == 0:
            self.net_dict['epsp']['use'] = False
        else:
            self.net_dict['epsp']['use'] = True

    def set_ipsp(self, lyr_specific):
        if int(lyr_specific) == 0:
            self.net_dict['ipsp']['use'] = False
        else:
            self.net_dict['ipsp']['use'] = True

    def set_ctsp(self, ctsp):
        if int(ctsp) == 0:
            self.special_dict['ctsp'] = False
        else:
            self.special_dict['ctsp'] = True

    def set_stp(self, stp):
        stp = int(stp)
        if 0 <= stp < len(self.stps):
            stp_name = list(self.stps)[stp]
            self.special_dict['stp_dict'] = self.stps[stp_name]
            print('stp used = {}'.format(stp_name))

    def set_indgs(self, indgs):
        if len(indgs) == 4:
            exc, pv, som, vip = indgs[0], indgs[1], indgs[2], indgs[3]
            self.net_dict['K_ext'] = np.array([exc, pv, som, vip, exc, pv, som, exc, pv, som, exc, pv, som])

    def renew_conn(self, raw):
        self.net_dict['conn_probs'] = func.renew_conn(net_dict['conn_probs'], 'microcircuit/conn_probs/raw_{}.csv'.format(raw))

    def load_conn(self, conn):
        self.net_dict['conn_probs'] = np.loadtxt('microcircuit/conn_probs/conn_{}.csv'.format(conn), delimiter=',')

    def adjust_vip_conn(self, adjust):
        if int(adjust) != 0:
            # vip-to-som change to layer-no-specific
            self.net_dict['conn_probs'][[6, 9, 12], 3] = self.net_dict['conn_probs'][2, 3]

    def set_ucomp(self, input):
        if int(input) == 0:
            self.net_dict['U-compensate'] = False
        else:
            self.net_dict['U-compensate'] = True

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

# read parameters from given input string
def read_levels(in_str):
    out_str = os.path.basename(in_str).split('.')[0]    # filename with .
    out_list = np.array(out_str.split('_')).astype(int)
    return out_str, out_list

# run scanning
if __name__ == "__main__":
    # get output names from system input
    outs = sys.argv[5:]
    constants = sys.argv[1:5]

    # create ScanParams object
    scanparams = ScanParams()
    scanparams.set_constant()

    # constant parameters
    scanparams.load_conn(constants[0])
    scanparams.set_epsp(constants[1])
    scanparams.set_ucomp(constants[2])
    scanparams.adjust_vip_conn(constants[3])

    # parameters to be scanned
    for out in outs:
        lvls_str, lvls = read_levels(out)
        scanparams.set_path(out, lvls_str)
        scanparams.set_g(lvls[0])
        scanparams.set_bg(lvls[1])
        scanparams.set_ctsp(lvls[2])
        scanparams.set_stp(lvls[3])
        # scanparams.set_indgs(lvls)
        scanparams.save_pickle(out)
