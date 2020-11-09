import os
import sys
import copy
import json
import pickle
import numpy as np
from microcircuit.network_params import net_dict, net_update
from microcircuit.sim_params import sim_dict
from microcircuit.stimulus_params import stim_dict
import microcircuit.functions as func
from microcircuit.stp.stp_dicts import stps
np.set_printoptions(suppress=True, precision=4)

'''
objects
'''
class ScanParams:
    def __init__(self, indgs=[750,1500,500,1000]):
        self.net_dict = copy.deepcopy(net_dict)
        self.sim_dict = copy.deepcopy(sim_dict)
        self.stim_dict = copy.deepcopy(stim_dict)
        self.stps = copy.deepcopy(stps)
        self.set_constant(indgs)

    def load_pickle(self, p_path):
        with open(p_path, 'rb') as p:
            paradict = pickle.load(p)
        self.net_dict = paradict['net_dict']
        self.sim_dict = paradict['sim_dict']
        self.stim_dict = paradict['stim_dict']

    def set_constant(self, indgs):
        self.net_dict['g'] = -8.0
        self.net_dict['bg_rate'] = 4.0
        # self.net_dict['stp_dict'] = {}
        self.net_dict['stp_dict'] = copy.deepcopy(self.stps['stp_fitted_02.pickle'])
        self.set_indgs(indgs)
        self.load_conn(conn='7-15')
        self.set_vip2som(True)

    def set_weight(self, pre, post, factor):
        for i, prepop in enumerate(self.net_dict['populations']):
            for j, postpop in enumerate(self.net_dict['populations']):
                if pre in prepop and post in postpop:
                    self.net_dict['psp_means'][j, i] *= float(factor)

    def set_conn(self, pre, post, factor):
        for i, prepop in enumerate(self.net_dict['populations']):
            for j, postpop in enumerate(self.net_dict['populations']):
                if pre in prepop and post in postpop:
                    self.net_dict['conn_probs'][j, i] *= float(factor)

    def del_item(self, dict_stp, keys=None, keysets=None):
        if keys is not None:
            for key in keys:
                del dict_stp[key]
        if keysets is not None:
            for sets in keysets:
                del dict_stp[sets[0]][sets[1]]

    def save_pickle(self, pickle_path=None):
        if pickle_path is None:
            pickle_path = 'para_dict.pickle'
        all_dict = {
            'net_dict': self.net_dict,
            'sim_dict': self.sim_dict,
            'stim_dict': self.stim_dict,
        }
        with open(pickle_path, 'wb') as h:
            pickle.dump(all_dict, h)

    def set_path(self, pickle_path):
        self.sim_dict['data_path'] = pickle_path.replace('.pickle', '/')

    def set_g(self, g):
        self.net_dict['g'] = -(np.abs(g))

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
            self.net_dict['ctsp'] = False
        else:
            self.net_dict['ctsp'] = True

    def set_stp(self, stp):
        stp = int(stp)
        if 0 <= stp < len(self.stps):
            stp_name = list(self.stps)[stp]
            self.net_dict['stp_dict'] = self.stps[stp_name]
            print('stp used = {}'.format(stp_name))

    def set_exc(self, exc):
        self.net_dict['K_ext'][[0, 4, 7, 10]] = int(exc)

    def set_pv(self, pv):
        self.net_dict['K_ext'][[1, 5, 8, 11]] = int(pv)

    def set_som(self, som):
        self.net_dict['K_ext'][[2, 6, 9, 12]] = int(som)

    def set_vip(self, vip):
        self.net_dict['K_ext'][3] = int(vip)

    def set_indgs(self, indgs):
        exc, pv, som, vip = indgs[0], indgs[1], indgs[2], indgs[3]
        self.net_dict['K_ext'] = np.array([exc, pv, som, vip, exc, pv, som, exc, pv, som, exc, pv, som])

    def renew_conn(self, raw):
        self.net_dict['conn_probs'] = func.renew_conn(net_dict['conn_probs'], 'microcircuit/conn_probs/raw_{}.csv'.format(raw))

    def load_conn(self, conn='7-15'):
        self.net_dict['conn_probs'] = np.loadtxt('microcircuit/conn_probs/conn_{}.csv'.format(conn), delimiter=',')

    def set_vip2som(self, adjust):
        if int(adjust) != 0:
            # vip-to-som all the same across layers
            print('adjust vip-to-som conn.')
            self.net_dict['conn_probs'][[6, 9, 12], 3] = self.net_dict['conn_probs'][2, 3]
        else:
            self.load_conn()

    def set_ucomp(self, input):
        if int(input) == 0:
            self.net_dict['U-compensate'] = False
        else:
            self.net_dict['U-compensate'] = True

    # read parameters from given input string
    def read_levels(self, in_str):
        out_str = os.path.basename(in_str).replace('.pickle', '')    # filename with .
        out_list = np.array(out_str.split('_')).astype(float)
        return out_str, out_list

'''
preliminary settings
'''
def set_perturb(para_dict, perturb_type, n_rep, pops, levels, start, duration, intrv, amp=0.1, freq=10.):
    if n_rep == 0:
        return 0.
    para_dict['stim_dict']['perturbs']['type'] = perturb_type
    para_dict['stim_dict']['perturbs']['levels'] = levels
    para_dict['stim_dict']['perturbs']['n_repeat'] = n_rep
    para_dict['stim_dict']['perturbs']['duration'] = duration
    para_dict['stim_dict']['perturbs']['targets'] = pops
    # handle starts
    starts_dict, n_levels, seg = {}, len(levels), duration+intrv
    lvl_str = ''
    for level in levels:
        lvl_str += str(int(level)) + '-'
    starts_fn = 'perturb_{:.0f}_{}_{}_{:.0f}_{:.0f}.json'.format(start, lvl_str, n_rep, duration, intrv)
    starts_fn = os.path.join(para_dict['sim_dict']['data_path'], starts_fn)
    if os.path.isfile(starts_fn):
        with open(starts_fn, 'r') as jf:
            starts_dict = json.load(jf)
    else:
        for level in levels:
            starts_dict[str(level)] = []
        for i in range(n_rep):
            tmps = np.array([start + i*n_levels*seg + j*seg for j in range(n_levels)])
            np.random.shuffle(tmps)
            for k, tmp in enumerate(tmps):
                starts_dict[str(levels[k])].append(tmp)
        with open(starts_fn, 'w') as jf:
            json.dump(starts_dict, jf)
    para_dict['stim_dict']['perturbs']['starts'] = starts_dict
    para_dict['stim_dict']['perturbs']['interval'] = intrv
    if perturb_type == 'ac':
        para_dict['stim_dict']['perturbs']['ac']['amplitude'] = amp
        para_dict['stim_dict']['perturbs']['ac']['frequency'] = freq
    return n_rep*n_levels*(duration+intrv)

# set layer-specific thalamic input
def set_thalamic(para_dict, th_starts=None, th_rate=None, orient=False, duration=10, conn_probs=None):
    th_dict = {}
    if type(th_starts) is list and len(th_starts) > 0 and type(th_rate) is float:
        if isinstance(conn_probs, np.ndarray) is not True:
            # Bruno, Simons, 2002; Oberlaender et al., 2011; Sermet et al., 2019; Constantinople, Bruno, 2013
            conn_probs = np.array([0.062, 0.062, 0.0, 0.0, 0.4, 0.4, 0.0, 0.259, 0.259, 0.0, 0.09, 0.09, 0.0])
        para_dict['stim_dict']['thalamic_input'] = True
        para_dict['stim_dict']['th_rate'] = th_rate
        para_dict['stim_dict']['th_start'] = np.array(th_starts).astype(float)
        para_dict['stim_dict']['th_duration'] = duration
        para_dict['stim_dict']['conn_probs_th'] = conn_probs
    para_dict['net_dict']['orient_tuning'] = orient

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
                    f.write('\n')
                else:
                    f.write(k2 + ':\n')
                    f.write('{}\n'.format(v2))
                    f.write('\n')
            f.close()

# run scanning
if __name__ == "__main__":
    # get output names from system input
    outs = sys.argv[1:]

    # create ScanParams object
    scanparams = ScanParams()

    # parameters to be scanned
    for out in outs:
        print(out)
        scanparams.set_path(out)
        scanparams.save_pickle(out)
