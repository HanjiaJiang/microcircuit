'''
the fitting
'''
import os
import sys
import json
import pickle
import numpy as np
from stp_test import ConnTest
from stp_tools import using_multiindex, colormap
# from microcircuit.functions import verify_collect, verify_print
np.set_printoptions(precision=2, linewidth=500, suppress=True)

if __name__ == '__main__':
    # read experimental data
    in_file = sys.argv[1]
    with open(in_file, 'rb') as p:
        exp_dict = pickle.load(p)

    # run if pickle not empty
    if exp_dict == {}:
        pass
    else:
        # get input parameters
        Uin = np.array(sys.argv[2].split('-')).astype(float)
        Fin = np.array(sys.argv[3].split('-')).astype(float)
        Din = np.array(sys.argv[4].split('-')).astype(float)
        Us = np.arange(Uin[0], Uin[1]+1e-05, Uin[2])
        Fs = np.arange(Fin[0], Fin[1]+1e-05, Fin[2])
        Ds = np.arange(Din[0], Din[1]+1e-05, Din[2])
        Ds[Ds==0.0] = 0.01
        n_bestpoints = int(sys.argv[5])

        # spike number is constant
        spk_n = 10

        # initialize ConnTest object; syn_dict not meaningful
        syn_dict = {'model': 'tsodyks_synapse', 'U': 0.5, 'tau_fac': 100, 'tau_rec': 100}
        conntest = ConnTest(syn_dict, "Exc", "Exc")

        # assign experimental parameters and run
        pre_subtype = exp_dict['pre_subtype']
        post_subtype = exp_dict['post_subtype']
        prog_fn = 'prog.{}.txt'.format('-'.join([pre_subtype, post_subtype]))
        spk_isi = exp_dict['spk_isi']
        pprs = exp_dict['pprs']
        peaks = exp_dict['peaks']
        # caches to save fitness; U: z, D: y, F: x
        fits = np.full((len(Us), len(Ds), len(Fs)), np.nan)
        min_fit, max_fit = np.inf, -np.inf
        for i, U in enumerate(Us):
            with open(prog_fn, 'a') as f:
                f.write('{:.2f}\n'.format(U))
            for j, F in enumerate(Fs):
                for k, D in enumerate(Ds):
                    # STP dictionary
                    syn_dict = {
                        'model': 'tsodyks_synapse',
                        'U': U,
                        'tau_fac': F,
                        'tau_rec': D
                    }
                    # initiate and run
                    conntest.setup(syn_dict,
                                        pre_subtype,
                                        post_subtype,
                                        pprs=pprs,
                                        peaks=peaks,
                                        spk_n=spk_n,
                                        spk_isi=spk_isi,
                                        verify=False)
                    conntest.run_sim(spk_n*spk_isi*1.5)
                    fit = conntest.run_analysis()
                    fits[i, k, j] = fit
                    if fit < min_fit:
                        min_fit = fit
                    if fit > max_fit:
                        max_fit = fit

        # fitness result handling
        df = using_multiindex(fits, list('ZYX'), 'fitness')
        print('mem by df = {}'.format(df.memory_usage(index=True).sum()))
        df = df.sort_values(by=['fitness'])[:min(n_bestpoints, len(df))]
        best_xyzs = np.array([df['X'].to_numpy(), df['Y'].to_numpy(), df['Z'].to_numpy()]).T
        best_UFDs = []
        for xyz in best_xyzs:
            best_UFDs.append([Us[xyz[2]], Fs[xyz[0]], Ds[xyz[1]]])
        best_UFDs = np.array(best_UFDs)
        df.insert(4, 'U', best_UFDs[:, 0])
        df.insert(5, 'F', best_UFDs[:, 1])
        df.insert(6, 'D', best_UFDs[:, 2])

        # output
        conn_name = '_'.join([pre_subtype, post_subtype, exp_dict['article']])
        print('connection: {}'.format(conn_name))
        print(df)
        csv_fn = 'stp:' + conn_name + '.csv'
        df.to_csv(csv_fn)
        plot_fn = colormap('stp:' + conn_name, 'fitness(RMSE)', fits, Fs, Ds, Us, (min_fit, max_fit), xlbl='F', ylbl='D', points=best_UFDs)

        # fitted STPs dump to pickle
        pname = 'stp:{}.pickle'.format(conn_name)
        # prepare dictionary
        fitted_dict = {
            'pre_subtype': pre_subtype,
            'post_subtype': post_subtype,
            'syn_dict':{
                'model': 'tsodyks_synapse',
                'U': best_UFDs[0, 0],
                'tau_fac': best_UFDs[0, 1],
                'tau_rec': best_UFDs[0, 2],
            }
        }
        # write
        with open(pname, 'wb') as h:
            pickle.dump(fitted_dict, h)
            h.close()
        # os.system('mv {} stp-data/'.format(prog_fn))
        # os.system('mv {} {} stp-results/'.format(csv_fn, plot_fn))
