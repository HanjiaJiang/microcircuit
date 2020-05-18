import os
import sys
import pickle
import numpy as np

if __name__ == '__main__':
    # list of .dat files
    in_list = sys.argv[1:-1]
    # article name
    article = sys.argv[-1]
    # fitness by each parameter set
    fits = []
    # minimum (best) of fitness
    min_fit = np.inf
    idx_min = 0
    # find minimum
    for i, in_file in enumerate(in_list):
        with open(in_file, 'r') as f:
            fit = float(f.read())
            if fit < min_fit:
                min_fit = fit
                idx_min = i
            fits.append(fit)
            f.close()
    best_set = in_list[idx_min].split('_')[-4:-1]
    best_str = '{},{},{}'.format(best_set[0], best_set[1], best_set[2])
    print('Best set: {}\nFitness = {}'.format(best_str, min_fit))

    params_list = in_list[0].split('_')

    # fitted STPs dump to pickle
    fname = 'stp_fitted.pickle'
    pretype = params_list[1].replace('-', '_')
    posttype = params_list[2].replace('-', '_')
    # check pickle exist
    if not os.path.isfile(fname):
        with open(fname, 'wb') as h:
            pickle.dump({}, h)
            h.close()
    # prepare dictionary
    stp_dicts = {}
    fitted_dict = {
        'model': 'tsodyks_synapse',
        'U': float(best_set[0])*0.1,
        'tau_fac': float(best_set[1]),
        'tau_rec': float(best_set[2]),
    }
    # read and add
    with open(fname, 'rb') as h:
        stp_dicts = pickle.load(h)
        if pretype not in stp_dicts:
            stp_dicts[pretype] = {}
        stp_dicts[pretype][posttype] = fitted_dict
        h.close()
    # write
    with open(fname, 'wb') as h:
        pickle.dump(stp_dicts, h)
        h.close()

    # create output directory and copy files
    out_dir = '_'.join(params_list[:7])
    os.system('mkdir -p {}'.format(out_dir))
    os.system('> {}/best-set={}'.format(out_dir, best_str))
    os.system('> {}/article={}'.format(out_dir, article))
    # os.system('mv STP*.dat {}/'.format(out_dir))
