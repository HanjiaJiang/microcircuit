import os
import sys
import copy
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size'] = 15.0
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes

class ScanData:
    # directory list, dimension dictionary
    def __init__(self, inputs, dims=None, plotvars=None):
        self.lyrs = ['l23', 'l4', 'l5', 'l6']
        if plotvars is None:
            self.plotvars = ['fr-exc', 'corr', 'cvisi', 'fr-pv', 'fr-som', 'fr-vip']
        else:
            self.plotvars = plotvars
        self.mtxs = {}
        self.fits = {}
        self.vbounds = {'fr-exc': [0.0, 20.0],
                        'corr': [-0.02, 0.02],
                        'cvisi': [0.0, 1.5],
                        'fr-pv': [0.0, 40.0],
                        'fr-som': [0.0, 40.0],
                        'fr-vip': [0.0, 40.0]}
        self.criteria = {'fr-exc': [0.0, 10.0],
                            'corr': [0.0001, 0.008],
                            'cvisi': [0.76, 1.2]
                            }
        self.setup(inputs, dims)

    def setup(self, inputs, dims):
        # determine the order of plot dimenstions
        if isinstance(dims, dict) and 'x' in dims and 'y' in dims and 'za' in dims and 'zb' in dims:
            self.dims = copy.deepcopy(dims)
        else:
            print('\'dims\' input not a good dictionary, using:')
            self.dims = {
                'x': 'exc',
                'y': 'pv',
                'za': 'som',
                'zb': 'vip',
            }
            print(self.dims)
        # make DataFrame object
        data_list = []
        excs, pvs, soms, vips = [0, 4, 7, 10], [1, 5, 8, 11], [2, 6, 9, 12], [3, 3, 3, 3]
        for p in inputs:
            path_str = os.path.dirname(p)
            params_list = np.array(path_str.split('/')[-1].split('_')).astype(float)
            fr_arr = np.loadtxt(os.path.join(path_str, 'fr.dat'), delimiter=',')
            ai_arr = np.loadtxt(os.path.join(path_str, 'ai.dat'), delimiter=',')
            # ain_arr = np.loadtxt(os.path.join(path_str, 'ai_n.dat'), delimiter=',')
            for i, lyr in enumerate(self.lyrs):
                tmp_dict = {
                    self.dims['x']: params_list[0],
                    self.dims['y']: params_list[1],
                    self.dims['za']: params_list[2],
                    self.dims['zb']: params_list[3],
                    'lyr': lyr,
                    'fr-exc': fr_arr[excs[i], 0],
                    'fr-pv': fr_arr[pvs[i], 0],
                    'fr-som': fr_arr[soms[i], 0],
                    'fr-vip': fr_arr[vips[i], 0],
                    'corr': ai_arr[i, 0],
                    'cvisi': ai_arr[i, 1],
                    # 'corr.n': ain_arr[i, 0],
                    # 'cvisi.n': ain_arr[i, 1],
                    }
                data_list.append(tmp_dict)
        self.df = pd.DataFrame(data_list)
        self.update_data(dims)

    # set plot data: dimensions, coordinates, matrix
    def update_data(self, dims):
        self.update_dims(dims)
        xname, yname, zaname, zbname = self.dims['x'], self.dims['y'], self.dims['za'], self.dims['zb']
        xs, ys, zas, zbs = self.df[xname].tolist(), self.df[yname].tolist(), self.df[zaname].tolist(), self.df[zbname].tolist() # float
        lvls_s = self.get_lvls([xs, ys, zas, zbs])
        xlvls, ylvls, zalvls, zblvls = lvls_s[0], lvls_s[1], lvls_s[2], lvls_s[3]
        x_intrv, y_intrv = (xlvls[1] - xlvls[0]), (ylvls[1] - ylvls[0])
        self.extent = [xlvls[0] - x_intrv/2, xlvls[-1] + x_intrv/2, ylvls[0] - y_intrv/2, ylvls[-1] + y_intrv/2]
        self.x_lvls, self.y_lvls, self.za_lvls, self.zb_lvls = xlvls, ylvls, zalvls, zblvls
        # print('extend={}'.format(self.extent))
        print('levels=')
        print(self.x_lvls, self.y_lvls, self.za_lvls, self.zb_lvls)
        self.make_data()
        self.make_plots()

    # update dimensions
    def update_dims(self, dims):
        if isinstance(dims, dict):
            self.dims = deepcopy(dims)

    # get dimension levels
    def get_lvls(self, lists):
        lvls_s = []
        for l in lists:
            lvls_s.append(sorted(list(set(l))))
        return lvls_s

    # make plot matrix
    def make_data(self):
        # structure
        for zb in self.zb_lvls:
            self.mtxs[str(zb)] = {}
            self.fits[str(zb)] = {}
            for za in self.za_lvls:
                self.mtxs[str(zb)][str(za)] = {}
                self.fits[str(zb)][str(za)] = {}
                for plotvar in self.plotvars:
                    self.mtxs[str(zb)][str(za)][plotvar] = np.full((len(self.lyrs), len(self.y_lvls), len(self.x_lvls)), np.nan)
                    if plotvar in self.criteria:
                        self.fits[str(zb)][str(za)][plotvar] = np.full((len(self.lyrs), len(self.y_lvls), len(self.x_lvls)), 0)
            # print(self.mtxs[str(zb)])
        # data
        lyrs = self.df['lyr'].tolist()
        for plotvar in self.plotvars:
            for i, value in enumerate(self.df[plotvar].tolist()):
                x = self.df[self.dims['x']][i]
                y = self.df[self.dims['y']][i]
                za = self.df[self.dims['za']][i]
                zb = self.df[self.dims['zb']][i]
                value = self.df[plotvar][i]
                lyr = self.lyrs.index(lyrs[i])
                idx_x, idx_y = self.x_lvls.index(x), self.y_lvls.index(y)
                self.mtxs[str(zb)][str(za)][plotvar][lyr, idx_x, idx_y] = value
                if plotvar in self.criteria:
                    if self.criteria[plotvar][0] <= value <= self.criteria[plotvar][1]:
                        self.fits[str(zb)][str(za)][plotvar][lyr, idx_x, idx_y] = 1

    def make_plots(self):
        # plot
        for zb in self.zb_lvls:
            for za in self.za_lvls:
                self.colormap(za, zb)

    def colormap(self, za, zb):
        # set plotting variables
        fig, axs = plt.subplots(4, len(self.plotvars), figsize=(16, 8), sharex=True, sharey=True, constrained_layout=True)
        xs, ys = np.array(self.x_lvls), np.array(self.y_lvls)
        plt.xlim((xs[0], xs[-1]))
        plt.ylim((ys[0], ys[-1]))
        plt.xticks(xs, rotation=30)
        xlbl, ylbl = self.dims['x'], self.dims['y']
        # plot
        for c, plotvar in enumerate(self.plotvars):
            for r in range(4):
                ax = axs[r, c]
                # vip
                if plotvar == 'fr-vip' and r > 0:
                    ax.axis('off')
                    continue

                data = self.mtxs[str(zb)][str(za)][plotvar][r].T
                # print(data)
                cs = ax.imshow(data, interpolation='none',
                                        cmap='Blues',
                                        origin='lower',
                                        extent=self.extent,
                                        vmin=self.vbounds[plotvar][0],
                                        vmax=self.vbounds[plotvar][1])

                # fitness points
                if plotvar in self.criteria:
                    fits = np.ones(data.shape)
                    for k, v in self.criteria.items():
                        fits_k = self.fits[str(zb)][str(za)][k][r].T
                        fits = np.multiply(fits, fits_k)
                        # individual criteria
                        if plotvar == k:
                            idxs_y, idxs_x = np.where(fits_k == 1)
                            xs_fit, ys_fit = xs[idxs_x], ys[idxs_y]
                            ax.scatter(xs_fit, ys_fit, s=20, color='y', zorder=9)
                    # all criteria together
                    idxs_y, idxs_x = np.where(fits == 1)
                    xs_fit, ys_fit = xs[idxs_x], ys[idxs_y]
                    ax.scatter(xs_fit, ys_fit, s=20, color='r', zorder=10)

                # set off-limit colors
                cs.cmap.set_over('purple')
                cs.cmap.set_under('magenta')

                # set plot aspect ratio
                ax.set_aspect(float((xs[-1] - xs[0])/(ys[-1] - ys[0])))

                # title
                if r == 0:
                    ax.set_title(plotvar)
                    cbar = fig.colorbar(cs, ax=axs[:, c], orientation='horizontal', shrink=0.6, aspect=5)

                # xlabel
                if r == 3:
                    ax.set_xlabel(xlbl)

                # ylabel
                if c == 0:
                    ax.set_ylabel(ylbl)
        #
        plot_name = '{}={},{}={}'.format(self.dims['za'], str(za), self.dims['zb'], str(zb))
        fig.suptitle(plot_name)
        fig.savefig(plot_name + '.png')
        plt.close()

if __name__ == '__main__':
    inputs = sys.argv[1:]
    scandata = ScanData(inputs)
