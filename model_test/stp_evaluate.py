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
            axs[r, c].set_title('U = {:.2f}'.format(float(zs[idx])))

    # x, y labels
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    # colorbar
    cbar = fig.colorbar(cs, ax=axs, orientation='vertical', shrink=0.6)

    fig.suptitle(name)
    plot_name = '{}__{}.png'.format(prefix, name)
    fig.savefig(plot_name)
    plt.close()
    return plot_name

if __name__ == '__main__':
    # article name
    article = sys.argv[1]

    # list of .dat files
    in_list = sys.argv[2:]

    # get data dimensions
    idx_sets = {'U': [], 'F': [], 'D': []}
    for in_file in in_list:
        buff = in_file.split('_')
        if buff[0] not in idx_sets['U']:
            idx_sets['U'].append(buff[0])
        if buff[1] not in idx_sets['F']:
            idx_sets['F'].append(buff[1])
        if buff[2] not in idx_sets['D']:
            idx_sets['D'].append(buff[2])
    data = np.full((len(idx_sets['U']), len(idx_sets['D']), len(idx_sets['F'])), np.nan)
    # minU = np.min(np.array(idx_sets['U']).astype(int))
    # minF = np.min(np.array(idx_sets['F']).astype(int))
    # minD = np.min(np.array(idx_sets['D']).astype(int))
    print('index set: U, F, D =\n{},\n{},\n{}'.format(idx_sets['U'], idx_sets['F'], idx_sets['D']))

    # get data
    xs, ys, zs = [], [], []
    fits, min_fit, min_i = [], np.inf, 0
    U, F, D = None, None, None
    for i, in_file in enumerate(in_list):
        with open(in_file, 'r') as f:
            buff = f.readlines()
        # params: x:F, y:D, z:U
        x, y, z = float(buff[-3]), float(buff[-2]), float(buff[-4])
        if x not in xs:
            xs.append(x)
        if y not in ys:
            ys.append(y)
        if z not in zs:
            zs.append(z)
        # fitness
        fit = float(buff[-1])
        if fit < min_fit:
            min_fit, min_i = fit, i
            U, F, D = z, x, y
        fits.append(fit)
        idxs = in_file.split('_')
        # print(idxs)
        # idxU, idxF, idxD = int(idxs[0]), int(idxs[1]), int(idxs[2])
        try:
            data[idx_sets['U'].index(idxs[0]), idx_sets['D'].index(idxs[2]), idx_sets['F'].index(idxs[1])] = fit
        except IndexError:
            print('Error: idxs={}')
    xs = sorted(xs)
    ys = sorted(ys)
    zs = sorted(zs)

    # verify
    print('xs (F) = {}'.format(np.array(xs)))
    print('ys (D) = {}'.format(np.array(ys)))
    print('zs (U) = {}'.format(np.array(zs)))
    print('data =')
    for d in data:
        print(d)

    # output from the best set
    best_str = 'U={:.2f},F={},D={}'.format(U, int(F), int(D))
    print('Best set: {}\nBest fitness = {}'.format(best_str, min_fit))
    best_fname = in_list[np.array(fits).argmin()]
    with open(best_fname, 'r') as f:
        reading = f.read().splitlines()
    print('data of the best set=\n{}'.format(reading))
    constants = '_'.join(reading[0:6]).replace(' ', '').replace('[','').replace(']','')
    plot_name = colormap(constants, 'fitness:RMSE', data[:len(zs), :len(ys), :len(xs)], xs, ys, zs, (min(fits), max(fits)), xlbl='F', ylbl='D')

    # fitted STPs dump to pickle
    fname = '../stp_fitted.pickle'
    pretype = reading[0].replace('-', '_')
    posttype = reading[1].replace('-', '_')
    # check pickle exist
    if not os.path.isfile(fname):
        with open(fname, 'wb') as h:
            pickle.dump({}, h)
            h.close()
    # prepare dictionary
    stp_dicts = {}
    if D == 0:
        D = 0.01
    fitted_dict = {
        'model': 'tsodyks_synapse',
        'U': U,
        'F': F,
        'D': D,
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
    out_dir = constants
    os.system('mkdir -p {} {}/png/'.format(out_dir, out_dir))
    os.system('> {}/best-set:{}'.format(out_dir, best_str))
    os.system('> {}/article:{}'.format(out_dir, article))
    os.system('cp Snakefile *.json *.yml *.sh stp*.py {} *{}* {}/'.format(plot_name, best_fname.replace('.dat', ''), out_dir))
    os.system('mv stp_*.png {}/png/'.format(out_dir))
