'''
get exp data to be fitted and dump to pickle
'''
import sys
import json
import pickle
import numpy as np

if __name__ == '__main__':
    in_file = sys.argv[1]
    outs = sys.argv[2:]
    with open(in_file, 'r') as j:
        exp_dicts = json.load(j)
    ks = []
    for i, k in enumerate(exp_dicts.keys()):
        ks.append(k)
    for i, out in enumerate(outs):
        dict_i = {}
        if i < len(ks):
            v = exp_dicts[ks[i]]
            try:
                pprs = np.array(v['pprs'].split('-')).astype(float)
            except ValueError:
                pprs = None
            try:
                peaks = np.array(v['peaks'].split('-')).astype(float)
            except ValueError:
                peaks = None
            dict_i = {
                'article': v['article'],
                'pre_subtype': v['pre_subtype'],
                'post_subtype': v['post_subtype'],
                'spk_isi': float(v['spk_isi']),
                'pprs': pprs,
                'peaks': peaks
            }
        with open(out, 'wb') as p:
            pickle.dump(dict_i, p)
