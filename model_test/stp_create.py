import sys
import numpy as np
import pickle

if __name__ == '__main__':
    out_list = sys.argv[1:]
    for out_file in out_list:
        params = out_file.split('_')[1:-1]
        try:
            pprs = np.array(params[4].split('-')).astype(float)
        except ValueError:
            pprs = None
        try:
            peaks = np.array(params[5].split('-')).astype(float)
        except ValueError:
            peaks = None
        if len(params) > 8:
            out_dict = {
                'pre_subtype': params[0],
                'post_subtype': params[1],
                'spk_n': int(params[2]),
                'spk_isi': float(params[3]),
                'pprs': pprs,
                'peaks': peaks,
                'U': float(params[6])*0.01,
                'tau_fac': float(params[7]),
                'tau_rec': float(params[8])
            }
        with open(out_file, 'wb') as h:
            pickle.dump(out_dict, h)
