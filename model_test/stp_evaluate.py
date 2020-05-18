import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

def colormap(prefix, name, data, xs, ys, v_range, xlbl='x', ylbl='y', cmap='RdBu'):
    xs = np.array(xs)
    ys = np.array(ys)

    # set plotting variables
    fig, axs = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
    # axs = axs.ravel()
    plot_name = '{}_{}.png'.format(prefix, name)

    # color range
    vmin = v_range[0]
    vmax = v_range[1]

    plt.xlabel(xlbl)
    plt.ylabel(ylbl, va='center', rotation='vertical')
    # plt.xticks(xs, rotation=30)
    plt.xlim((xs[0], xs[-1]))
    plt.ylim((ys[0], ys[-1]))

    # define plot borders
    extent = [np.min(xs), np.max(xs), np.min(ys), np.max(ys)]

    # colormap plots
    cs = axs.imshow(data, interpolation='none', cmap=cmap, origin='lower',
                       extent=extent, vmax=vmax, vmin=vmin)

    # set off-limit colors
    cs.cmap.set_over("midnightblue")
    cs.cmap.set_under("firebrick")

    # set plot aspect ratio
    axs.set_aspect(float((xs[-1] - xs[0])/(ys[-1] - ys[0])))

    # set yticklabels
    # axs.set_yticklabels(ys, rotation=30)

    # colorbar
    cbar = fig.colorbar(cs, ax=axs, orientation='horizontal', shrink=0.6)

    fig.suptitle(name)
    fig.savefig(plot_name)
    plt.close()
    return plot_name

if __name__ == '__main__':
    # list of .dat files
    in_list = sys.argv[1:-1]
    # article name
    article = sys.argv[-1]
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
            f.close()

    # print the best set
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
    os.system('cp Snakefile {}/'.format(out_dir))
    os.system('cp stp*.py {}/'.format(out_dir))
    # os.system('mv STP*.dat {}/'.format(out_dir))
