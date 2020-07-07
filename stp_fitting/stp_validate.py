'''
for the final validation of fitted stp data
'''
import json
import pickle
import numpy as np
from stp_test import ConnTest

if __name__ == '__main__':
    spk_n = 10
    with open('stp_fitted_02.pickle', 'rb') as p:
        stp = pickle.load(p)
    with open('exp-data_02.json', 'r') as js:
        exp = json.load(js)
    for prek, prev in stp.items():
        for postk, postv in prev.items():
            pre_subtype, post_subtype = prek.replace('_', '-'), postk.replace('_', '-')
            spk_isi = pprs = peaks = None
            U, F, D = postv['U'], postv['tau_fac'], postv['tau_rec']
            syn_dict = {
                'model': 'tsodyks_synapse',
                'U': U,
                'tau_fac': F,
                'tau_rec': D
            }
            for k, v in exp.items():
                if v['pre_subtype'] == pre_subtype and v['post_subtype'] == post_subtype:
                    spk_isi = float(v['spk_isi'])
                    try:
                        pprs = np.array(v['pprs'].split('-')).astype(float)
                    except ValueError:
                        pprs = None
                    try:
                        peaks = np.array(v['peaks'].split('-')).astype(float)
                    except ValueError:
                        peaks = None
            conntest = ConnTest(syn_dict,
                                pre_subtype,
                                post_subtype,
                                pprs=pprs,
                                peaks=peaks,
                                spk_n=spk_n,
                                spk_isi=spk_isi,
                                verify=True)
            conntest.run_sim(spk_n*spk_isi*1.5)
            conntest.run_analysis()
