import os
import sys
import pickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size'] = 20.0
np.set_printoptions(precision=4, linewidth=500, suppress=True)

def colormap(prefix, name, data, xs, ys, zs, v_range, xlbl='x', ylbl='y', cmap='Blues_r'):
    # set plotting variables
    xs = np.array(xs)
    ys = np.array(ys)
    n_col = 5
    n_row = 4
    fig, axs = plt.subplots(n_row, n_col, figsize=(12, 12), constrained_layout=True)

    # color range
    vmin = v_range[0]
    vmax = v_range[1]

    # plt.xlabel(xlbl)
    # plt.ylabel(ylbl, va='center', rotation='vertical')
    # plt.xticks(xs, rotation=30)
    plt.xlim((xs[0], xs[-1]))
    plt.ylim((ys[0], ys[-1]))

    # define plot borders
    extent = [np.min(xs), np.max(xs), np.min(ys), np.max(ys)]

    for r in range(n_row):
        for c in range(n_col):
            idx = r*n_col+c
            if idx >= len(data):
                axs[r, c].axis('off')
                continue
            # colormap plots
            cs = axs[r, c].imshow(data[idx], interpolation='none',
                                cmap=cmap, origin='lower',
                                extent=extent, vmax=vmax, vmin=vmin)

            # set off-limit colors
            cs.cmap.set_over("midnightblue")
            cs.cmap.set_under("firebrick")

            # set plot aspect ratio
            axs[r, c].set_aspect(float((xs[-1] - xs[0])/(ys[-1] - ys[0])))

            # title
            axs[r, c].set_title('U = {:.2f}'.format(float(zs[idx])*0.01))

    # x, y labels
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    # colorbar
    cbar = fig.colorbar(cs, ax=axs, orientation='vertical', shrink=0.6)

    fig.suptitle(name)
    plot_name = '{}_{}.png'.format(prefix, name)
    fig.savefig(plot_name)
    plt.close()
    return plot_name

if __name__ == '__main__':
    # list of .dat files
    in_list = sys.argv[1:-1]
    params = in_list[0].split('_')

    # article name
    article = sys.argv[-1]

    # get colormap parameters
    xs = []
    ys = []
    zs = []
    for in_file in in_list:
        params_tmp = in_file.split('_')
        xs.append(float(params_tmp[-3]))
        ys.append(float(params_tmp[-2]))
        zs.append(float(params_tmp[-4]))
    xs = sorted(list(set(xs)))
    ys = sorted(list(set(ys)))
    zs = sorted(list(set(zs)))
    data = np.full((len(zs), len(ys), len(xs)), np.nan)

    # get fitness data
    fits = []
    idx_min = 0
    for i, in_file in enumerate(in_list):
        params_tmp = in_file.split('_')
        x = float(params_tmp[-3])   # tau_fac
        y = float(params_tmp[-2])   # tau_rec
        z = float(params_tmp[-4])
        with open(in_file, 'r') as f:
            fit = float(f.read())
            fits.append(fit)
            data[zs.index(z), ys.index(y), xs.index(x)] = fit
            f.close()

    # veryfy
    print('xs = {},{}'.format(xs, type(xs)))
    print('ys = {},{}'.format(ys, type(ys)))
    print('zs = {},{}'.format(zs, type(zs)))
    print('data =')
    for d in data:
        print(d)

    # print the best set
    best_fname = in_list[np.array(fits).argmin()].replace('.dat', '')
    best_set = best_fname.split('_')[-4:-1]   # U, tau_fac, tau_rec
    U = float(best_set[0])*0.01
    tau_fac = int(best_set[1])
    tau_rec = int(best_set[2])
    best_str = 'U={:.2f},tau_fac={},tau_rec={}'.format(U, tau_fac, tau_rec)
    print('Best set: {}\nBest fitness = {}'.format(best_str, np.min(fits)))

    # colormap
    prefix = '_'.join(params[:7])
    name = 'fitness'
    plot_name = colormap(prefix, name, data, xs, ys, zs, (min(fits), max(fits)), xlbl='tau_fac', ylbl='tau_rec')

    # fitted STPs dump to pickle
    fname = '../stp_fitted.pickle'
    pretype = params[1].replace('-', '_')
    posttype = params[2].replace('-', '_')

    # check pickle exist
    if not os.path.isfile(fname):
        with open(fname, 'wb') as h:
            pickle.dump({}, h)
            h.close()

    # prepare dictionary
    stp_dicts = {}
    if tau_rec == 0:
        tau_rec = 0.01
    fitted_dict = {
        'model': 'tsodyks_synapse',
        'U': U,
        'tau_fac': tau_fac,
        'tau_rec': tau_rec,
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
    out_dir = '_'.join(params[:7])
    os.system('mkdir -p {}'.format(out_dir))
    os.system('> {}/best-set:{}'.format(out_dir, best_str))
    os.system('> {}/article:{}'.format(out_dir, article))
    os.system('cp Snakefile stp*.py {} {}* {}/'.format(plot_name, best_fname, out_dir))
    os.system('mkdir -p {}/png_backup/'.format(out_dir))
    os.system('mv STP*.png {}/png_backup/'.format(out_dir))
