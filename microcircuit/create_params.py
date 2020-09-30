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
    def __init__(self, indgs=[750,1500,500,1250]):
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
        self.net_dict['stp_dict'] = {}
        self.net_dict['stp_dict'] = copy.deepcopy(self.stps['stp_fitted_02.pickle'])
        self.set_indgs(indgs)
        self.load_conn('7-15')
        self.vip2som(True)
        net_update(self.net_dict)

    # def do_single(self, pickle_path, indgs=None):
    #     # self.set_constant(indgs)
    #     # self.set_weight('Exc', 'Exc', 1.25)
    #     # self.set_weight('SOM', 'PV', 1.25)
    #     self.save_pickle(pickle_path)

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
        self.net_dict['g'] = -int(np.abs(int(g)))

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

    def load_conn(self, conn):
        self.net_dict['conn_probs'] = np.loadtxt('microcircuit/conn_probs/conn_{}.csv'.format(conn), delimiter=',')

    def vip2som(self, adjust):
        if int(adjust) != 0:
            # vip-to-som all the same across layers
            print('adjust vip conn.')
            self.net_dict['conn_probs'][[6, 9, 12], 3] = self.net_dict['conn_probs'][2, 3]

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
def set_paradox(para_dict, paradox_type, n, pops, offsets, start, duration, intrv, amp=0.1, freq=10.):
    para_dict['stim_dict']['paradox']['type'] = paradox_type
    para_dict['stim_dict']['paradox']['offsets'] = offsets
    para_dict['stim_dict']['paradox']['n'] = n
    para_dict['stim_dict']['paradox']['duration'] = duration
    para_dict['stim_dict']['paradox']['targets'] = pops
    # handle starts
    starts_dict, n_offset, seg = {}, len(offsets), duration+intrv
    for offset in offsets:
        starts_dict[str(offset)] = []
    for i in range(n):
        tmps = np.array([start + i*n_offset*seg + j*seg for j in range(n_offset)])
        np.random.shuffle(tmps)
        for k, tmp in enumerate(tmps):
            starts_dict[str(offsets[k])].append(tmp)
    # for i, offset in enumerate(offsets):
    #     starts_dict[str(offset)] = [start + i*n*(duration+intrv) + j*(duration+intrv) for j in range(n)]
    para_dict['stim_dict']['paradox']['starts'] = starts_dict
    para_dict['stim_dict']['paradox']['intrv'] = intrv
    if paradox_type == 'ac':
        para_dict['stim_dict']['paradox']['amplitude'] = amp
        para_dict['stim_dict']['paradox']['frequency'] = freq
    return n*len(offsets)*(duration+intrv)

# set layer-specific thalamic input
def set_thalamic(para_dict, th_starts=None, th_rate=None, orient=False, duration=10):
    th_dict = {}
    if type(th_starts) is list and len(th_starts) > 0 and type(th_rate) is float:
        para_dict['stim_dict']['thalamic_input'] = True
        para_dict['stim_dict']['th_rate'] = th_rate
        para_dict['stim_dict']['th_start'] = np.array(th_starts).astype(float)
        para_dict['stim_dict']['th_duration'] = duration
        # Bruno, Simons, 2002; Oberlaender et al., 2011; Sermet et al., 2019; Constantinople, Bruno, 2013
        para_dict['stim_dict']['conn_probs_th'] = np.array([0.062, 0.062, 0.0, 0.0, 0.4, 0.4, 0.0, 0.259, 0.259, 0.0, 0.09, 0.09, 0.0])
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

# run scanning
if __name__ == "__main__":
    # get output names from system input
    outs = sys.argv[1:]

    # create ScanParams object
    scanparams = ScanParams()
    # scanparams.set_constant()

    # constant parameters
    # scanparams.vip2som(True)

    # parameters to be scanned
    for out in outs:
        print(out)
        scanparams.set_path(out)
        scanparams.save_pickle(out)
