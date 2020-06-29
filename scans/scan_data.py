import os
import sys
import copy
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size'] = 15.0
np.set_printoptions(precision=3, linewidth=500, suppress=True)
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes

class ScanData:
    # directory list, dimension dictionary
    def __init__(self, inputs, dims=None, plotvars=None, figsize=(20, 10)):
        self.lyrs = ['L2/3', 'L4', 'L5', 'L6']
        self.figsize = figsize
        if plotvars is None:
            self.plotvars = ['fr-exc', 'corr', 'cvisi', 'fr-pv', 'fr-som', 'fr-vip']
        else:
            self.plotvars = plotvars
        self.mtxs = {}
        self.fits = {}
        self.criteria = {'fr-exc': [0.0, 2.0],
                            'corr': [0.0001, 0.1],
                            'cvisi': [0.6, 0.76]
                            }
        self.cmaps = {'fr-exc': 'Blues',
                        'corr': 'RdBu',
                        'cvisi': 'Blues',
                        'fr-pv': 'Blues',
                        'fr-som': 'Blues',
                        'fr-vip': 'Blues'}
        self.fmt = {'fr-exc': '%1.1f',
                        'corr': '%1.4f',
                        'cvisi': '%1.2f',
                        'fr-pv': '%1.1f',
                        'fr-som': '%1.1f',
                        'fr-vip': '%1.1f'}
        self.units = {'fr-exc': '(spikes/s)',
                        'corr': '',
                        'cvisi': '',
                        'fr-pv': '(spikes/s)',
                        'fr-som': '(spikes/s)',
                        'fr-vip': '(spikes/s)'}
        self.setup(inputs, dims)

    def setup(self, inputs, dims):
        # determine the order of plot dimenstions
        if isinstance(dims, list) and len(dims) == 4:
            pass
        else:
            print('using default dimension items: exc, pv, som, vip')
            dims = ['exc', 'pv', 'som', 'vip']
        self.dims = {
            'x': dims[0],
            'y': dims[1],
            'za': dims[2],
            'zb': dims[3],
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
        fig, axs = plt.subplots(4, len(self.plotvars), figsize=self.figsize, sharex=True, sharey=True)
        xs, ys = np.array(self.x_lvls), np.array(self.y_lvls)
        plt.xlim((xs[0], xs[-1]))
        plt.ylim((ys[0], ys[-1]))
        plt.yticks(ys, rotation=30)
        plt.xticks(xs, rotation=30)
        xlbl, ylbl = self.dims['x'], self.dims['y']
        # plot
        for c, plotvar in enumerate(self.plotvars):
            vmin, vmax = 0.0, np.nanmax(np.abs(self.mtxs[str(zb)][str(za)][plotvar]))
            if plotvar == 'corr':
                vmin = -vmax
            for r in range(4):
                ax = axs[r, c]

                # plot data
                data = self.mtxs[str(zb)][str(za)][plotvar][r].T
                if plotvar == 'fr-vip' and r > 0:
                    ax.axis('off')
                    data = np.full(data.shape, np.nan)
                # cf = ax.contourf(data,
                #     cmap=self.cmaps[plotvar],
                #     origin='lower',
                #     extent=self.extent,
                #     vmin=vmin,
                #     vmax=vmax)
                cf = ax.imshow(data, interpolation='bilinear',
                    cmap=self.cmaps[plotvar],
                    origin='lower',
                    extent=self.extent,
                    vmin=vmin,
                    vmax=vmax)

                # single- and triple-fit
                if plotvar in self.criteria:
                    tri_fit = np.ones(data.shape)
                    for k, v in self.criteria.items():
                        fit = self.fits[str(zb)][str(za)][k][r].T
                        tri_fit = np.multiply(tri_fit, fit)
                    # data_fit = fitness
                    data_fit = np.zeros(data.shape)
                    # individual
                    data_fit[np.where((data>self.criteria[plotvar][0])&(data<self.criteria[plotvar][1]))] = 10.0
                    data_fit[np.where(tri_fit==1)] = 20.0
                    print('L{}:\n{}\n{}'.format(r, data[::-1], data_fit[::-1]))
                    cf_fit = ax.contourf(data_fit,
                        levels=[5.0, 15.0, 25.0],
                        origin='lower',
                        colors='gray',
                        extent=self.extent,
                        hatches=['//', '++', '//'],
                        alpha=0.0)
                    # ax.clabel(cf_fit, cf_fit.levels, fmt=self.fmt[plotvar], inline=True, fontsize=10)

                # many settings
                ax.set_aspect(float((xs[-1] - xs[0])/(ys[-1] - ys[0])))
                ax.set_xticklabels(xs, rotation=30)
                ax.set_yticklabels(ys, rotation=30)

                # title
                if r == 0:
                    ax.set_title(plotvar + self.units[plotvar] + '\n ')

                # xlabel
                if r == 3:
                    ax.set_xlabel(xlbl)
                    cbar = fig.colorbar(cf, ax=axs[:, c], orientation='horizontal', shrink=0.8, aspect=10)
                    # cbar_xlbls = cbar.ax.get_xticklabels()
                    # print('cbar_xlbls = {}'.format(cbar_xlbls))
                    # cbar.ax.set_xticklabels(cbar_xlbls, rotation=30)

                # ylabel
                if c == 0:
                    ax.set_ylabel(ylbl)
                    # ax.text(0.5, 0.5, self.lyrs[r], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        plot_name = '{}={},{}={}'.format(self.dims['za'], str(za), self.dims['zb'], str(zb))
        # fig.suptitle('ai state and firint rates')
        fig.savefig(plot_name + '.png', bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    inputs = sys.argv[5:]
    dims = sys.argv[1:5]
    scandata = ScanData(inputs, dims=dims)
    # scandata = ScanData(inputs, plotvars=['fr-exc', 'corr', 'cvisi'], figsize=(10, 8))
