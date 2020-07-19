import numpy as np
import string
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size'] = 20.0

def using_multiindex(A, columns, data_name):
    shape = A.shape
    index = pd.MultiIndex.from_product([range(s)for s in shape], names=columns)
    df = pd.DataFrame({data_name: A.flatten()}, index=index).reset_index()
    return df

def colormap(prefix, name, data, xs, ys, zs, v_range, xlbl='x', ylbl='y', cmap='Blues_r', points=None):
    # set plotting variables
    xs = np.array(xs)
    ys = np.array(ys)
    n_col = 5
    n_row = 4
    fig, axs = plt.subplots(n_row, n_col, figsize=(12, 12),  sharex=True, sharey=True, constrained_layout=True)

    # color range
    vmin = v_range[0]
    vmax = v_range[1]

    plt.xlim((xs[0], xs[-1]))
    plt.ylim((ys[0], ys[-1]))

    # define plot borders
    extent = [np.min(xs), np.max(xs), np.min(ys), np.max(ys)]

    for r in range(n_row):
        for c in range(n_col):
            idx_u = r*n_col+c
            if idx_u >= len(data):
                axs[r, c].axis('off')
                continue
            # colormap plots
            cs = axs[r, c].imshow(data[idx_u], interpolation='none',
                                cmap=cmap, origin='lower',
                                extent=extent, vmax=vmax, vmin=vmin)
            # best-fit points
            U = float(zs[idx_u])
            if isinstance(points, np.ndarray):
                best_xs = points[points[:, 0]==U][:, 1]
                best_ys = points[points[:, 0]==U][:, 2]
                axs[r, c].scatter(best_xs, best_ys, color='r')
                best_value = data[idx_u, ys.tolist().index(best_ys[0]), xs.tolist().index(best_xs[0])]
                axs[r, c].text(best_xs[0], best_ys[0], '{:.4f}'.format(best_value), color='magenta', fontsize=12, horizontalalignment='center')

            # set off-limit colors
            cs.cmap.set_over("midnightblue")
            cs.cmap.set_under("firebrick")

            # set plot aspect ratio
            axs[r, c].set_aspect(float((xs[-1] - xs[0])/(ys[-1] - ys[0])))

            # title
            axs[r, c].set_title('U = {:.2f}'.format(U))


    # x, y labels
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    # colorbar
    cbar = fig.colorbar(cs, ax=axs, orientation='vertical', shrink=0.6)

    fig.suptitle(name)
    plot_name = '{}:{}.png'.format(prefix, name)
    fig.savefig(plot_name)
    plt.close()
    return plot_name
