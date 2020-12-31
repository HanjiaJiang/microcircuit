import os
import sys
import copy
import json
import pickle
import numpy as np
from microcircuit.network_params import net_dict
from microcircuit.sim_params import sim_dict
from microcircuit.stimulus_params import stim_dict
import microcircuit.functions as func
from microcircuit.stp.stp_dicts import stps
np.set_printoptions(suppress=True, precision=4)

'''
objects
'''
class ScanParams:
    def __init__(self, indgs=[750,1500,500,1000], conn='0715', g=-8., bg=4.5, stp=2, vip2som=True):
        # parameter dicts
        self.conn = conn
        self.net_dict = copy.deepcopy(net_dict)
        self.sim_dict = copy.deepcopy(sim_dict)
        self.stim_dict = copy.deepcopy(stim_dict)
        self.stps = copy.deepcopy(stps)
        # set basic parameters
        self.net_dict['g'] = g
        self.net_dict['bg_rate'] = bg
        self.net_dict['stp_dict'] = copy.deepcopy(self.stps[list(self.stps)[stp]])
        self.set_indgs(indgs)
        self.load_conn(fn=conn)
        self.set_vip2som(vip2som)

    def set_lognormal(self, ln):
        self.net_dict['recurrent_weight_distribution'] = 'lognormal' if ln is True else 'normal'

    # save the parameters to pickle (Snakefile)
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

    # load the parameters from pickle (Snakefile)
    def load_pickle(self, p_path):
        with open(p_path, 'rb') as p:
            paradict = pickle.load(p)
        self.net_dict = paradict['net_dict']
        self.sim_dict = paradict['sim_dict']
        self.stim_dict = paradict['stim_dict']

    # scale weight by a factor
    def factor_weight(self, pre, post, factor):
        for i, prepop in enumerate(self.net_dict['populations']):
            for j, postpop in enumerate(self.net_dict['populations']):
                if pre in prepop and post in postpop:
                    self.net_dict['psp_means'][j, i] *= float(factor)

    # scale conn. prob. by a factor
    def factor_conn(self, pre, post, factor):
        for i, prepop in enumerate(self.net_dict['populations']):
            for j, postpop in enumerate(self.net_dict['populations']):
                if pre in prepop and post in postpop:
                    self.net_dict['conn_probs'][j, i] *= float(factor)

    # delete specific item in STP dict
    def del_item(self, dict_stp, keys=None, key_tups=None):
        if keys is not None:
            for key in keys:
                del dict_stp[key]
        if key_tups is not None:
            for tup in key_tups:
                del dict_stp[tup[0]][tup[1]]

    # set data/pickle path
    def set_path(self, pickle_path):
        self.sim_dict['data_path'] = pickle_path.replace('.pickle', '/')

    def set_g(self, g):
        self.net_dict['g'] = -(np.abs(g))

    def set_bg(self, bg):
        self.net_dict['bg_rate'] = float(bg)

    def set_epsp(self, use):
        self.net_dict['epsp']['use'] = False if use == False else True

    def set_ipsp(self, use):
        self.net_dict['ipsp']['use'] = False if use == False else True

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

    def renew_conn(self, fn='0715', extrapolate=True):
        self.net_dict['conn_probs'] = func.renew_conn(
            'microcircuit/conn_probs/raw_{}.csv'.format(fn),
            extrapolate=extrapolate)

    def load_conn(self, fn='0715'):
        self.net_dict['conn_probs'] = np.loadtxt('microcircuit/conn_probs/conn_{}.csv'.format(fn), delimiter=',')

    def set_vip2som(self, adjust):
        if int(adjust) != 0:
            # vip-to-som all the same across layers
            print('adjust vip-to-som conn.')
            self.net_dict['conn_probs'][[6, 9, 12], 3] = self.net_dict['conn_probs'][2, 3]
        else:
            self.load_conn(self.conn)

    # compensate for U in STP
    def set_ucomp(self, ucomp):
        self.net_dict['U-compensate'] = False if ucomp == 0 else True

    # read parameters from given input string
    def read_levels(self, in_str):
        out_str = os.path.basename(in_str).replace('.pickle', '')    # filename with .
        out_list = np.array(out_str.split('_')).astype(float)
        return out_str, out_list

    # perturbs
    def set_perturb(self, para_dict, dict_in, amp=0.1, freq=10.):
        if dict_in['n_repeat'] == 0:
            return 0.
        para_dict['stim_dict']['perturbs'] = copy.deepcopy(dict_in)

        # handle json
        starts_dict, n_levels, seg = {}, len(dict_in['levels']), dict_in['duration']+dict_in['interval']
        lvl_str = ''
        for level in dict_in['levels']:
            lvl_str += str(int(level)) + '-'
        json_fn = 'perturb_{:.0f}_{}_{}_{:.0f}_{:.0f}.json'. \
            format(dict_in['start'], lvl_str[:-1], dict_in['n_repeat'],
            dict_in['duration'], dict_in['interval'])
        json_fn = os.path.join(para_dict['sim_dict']['data_path'], json_fn)
        # handle starts
        if os.path.isfile(json_fn):
            with open(json_fn, 'r') as jf:
                starts_dict = json.load(jf)
        else:
            for level in dict_in['levels']:
                starts_dict[str(level)] = []
            for i in range(dict_in['n_repeat']):
                tmps = np.array([dict_in['start'] + i*n_levels*seg + j*seg for j in range(n_levels)])
                np.random.shuffle(tmps)
                for k, tmp in enumerate(tmps):
                    starts_dict[str(dict_in['levels'][k])].append(tmp)
            with open(json_fn, 'w') as jf:
                json.dump(starts_dict, jf)
        para_dict['stim_dict']['perturbs']['starts'] = starts_dict
        # para_dict['stim_dict']['perturbs']['interval'] = intrv
        if dict_in['type'] == 'ac':
            para_dict['stim_dict']['perturbs']['ac']['amplitude'] = amp
            para_dict['stim_dict']['perturbs']['ac']['frequency'] = freq
        # return dict_in['n_repeat']*n_levels*seg

    # thalamic input
    def set_thalamic(self, para_dict, dict_in):
        if isinstance(dict_in['th_start'], np.ndarray) and len(dict_in['th_start']) > 0 \
            and type(dict_in['th_rate']) is float:
            if isinstance(dict_in['conn_probs_th'], np.ndarray) is not True:
                # Bruno, Simons, 2002; Oberlaender et al., 2011; Sermet et al., 2019; Constantinople, Bruno, 2013
                dict_in['conn_probs_th'] = np.array([0.062, 0.062, 0.0, 0.0, 0.4, 0.4, 0.0, 0.259, 0.259, 0.0, 0.09, 0.09, 0.0])
            para_dict['stim_dict']['thalamic_input'] = True
            para_dict['stim_dict'].update(dict_in)

    '''
    tools
    '''
    def print_all(self, all_dict):
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

# create pickles for scanning
if __name__ == "__main__":
    # get output names from system input
    outs = sys.argv[1:]

    # create ScanParams object
    scanparams = ScanParams()

    # parameters to be scanned
    for out in outs:
        scanparams.set_path(out)
        scanparams.save_pickle(out)
