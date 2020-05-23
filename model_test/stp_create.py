import sys
import numpy as np
import pickle

if __name__ == '__main__':
    # scanned parameters
    Us = np.arange(0.05, 1.01, 0.05)
    Fs = np.arange(0.0, 200.1, 10.0)
    Ds = np.arange(0.0, 200.1, 10.0)
    # constant parameters
    spk_n = 10
    pre_subtype = sys.argv[1]
    post_subtype = sys.argv[2]
    spk_isi = float(sys.argv[3])
    try:
        pprs = np.array(sys.argv[4].split('-')).astype(float)
    except ValueError:
        pprs = None
    try:
        peaks = np.array(sys.argv[5].split('-')).astype(float)
    except ValueError:
        peaks = None
    out_list = sys.argv[6:]
    for out_file in out_list:
        params = out_file.split('_')[:-1]
        if len(params) > 2:
            out_dict = {
                'pre_subtype': pre_subtype,
                'post_subtype': post_subtype,
                'spk_n': spk_n,
                'spk_isi': spk_isi,
                'pprs': pprs,
                'peaks': peaks,
                'U': Us[int(params[0])],
                'tau_fac': Fs[int(params[1])],
                'tau_rec': Ds[int(params[2])]
            }
        with open(out_file, 'wb') as h:
            pickle.dump(out_dict, h)
