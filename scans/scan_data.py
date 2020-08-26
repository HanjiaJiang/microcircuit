import os
import sys
import copy
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
matplotlib.rcParams['font.size'] = 15.0
np.set_printoptions(precision=3, linewidth=500, suppress=True)
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes

class ScanData:
    # directory list, dimension dictionary
    def __init__(self, inputs, dims=None, plotvars=None, criteria=None, figsize=(16, 12), mark_flg=False):
        self.set_params(plotvars, figsize, criteria, mark_flg)
        self.set_data(inputs, dims)

    def set_params(self, plotvars, figsize, criteria, mark_flg):
        self.lyrs = ['L2/3', 'L4', 'L5', 'L6']
        self.figsize = figsize
        if plotvars is None:
            self.plotvars = [r'$r_{Exc}$', 'pairwise\ncorrelation', 'CV(ISI)', r'$r_{PV}$', r'$r_{SOM}$', r'$r_{VIP}$']
        else:
            self.plotvars = plotvars
        self.mtxs, self.fits, self.rmse, self.rmse_n = {}, {}, {}, {}
        # ground state criteria
        if criteria is None:
            self.criteria = {r'$r_{Exc}$': [[0.0, 10.0], [0.0, 10.0], [0.0, 10.0], [0.0, 10.0]],
                                'pairwise\ncorrelation': [[0.0001, 0.008], [0.0001, 0.008], [0.0001,0.008], [0.0001, 0.008]],
                                'CV(ISI)': [[0.76, 1.2], [0.76, 1.2], [0.76, 1.2], [0.76, 1.2]]
                                }
        else:
            self.criteria = criteria
        # RMSE
        self.mark_flg = mark_flg
        self.rmse_criteria = {r'$r_{Exc}$': [2.7, 0.5, 6.8, 6.1],
                            r'$r_{PV}$': [13.8, 10.2, 7.5, 16.9],
                            r'$r_{SOM}$': [2.6, 2.6, 2.8, 3.9],
                            r'$r_{VIP}$': [14.6]}
        self.cmaps = {r'$r_{Exc}$': 'Blues',
                        'pairwise\ncorrelation': 'RdBu',
                        'CV(ISI)': 'Blues',
                        r'$r_{PV}$': 'Blues',
                        r'$r_{SOM}$': 'Blues',
                        r'$r_{VIP}$': 'Blues'}
        self.clabel_format = {r'$r_{Exc}$': '%1.0f',
                        'pairwise\ncorrelation': '%1.3f',
                        'CV(ISI)': '%1.2f',
                        r'$r_{PV}$': '%1.0f',
                        r'$r_{SOM}$': '%1.0f',
                        r'$r_{VIP}$': '%1.0f'}
        self.vmaxs = {r'$r_{Exc}$': 20.0,
                        'pairwise\ncorrelation': 0.1,
                        'CV(ISI)': 1.0,
                        r'$r_{PV}$': 50.0,
                        r'$r_{SOM}$': 50.0,
                        r'$r_{VIP}$': 50.0}

    def set_data(self, inputs, dims):
        # dimensions
        if isinstance(dims, list) and len(dims) == 4:
            pass
        else:
            dims = ['exc', 'pv', 'som', 'vip']
        self.dims = {'x': dims[0], 'y': dims[1], 'za': dims[2], 'zb': dims[3]}
        print('dimensions =\n{}'.format(self.dims))
        # make DataFrame object
        data_list = []
        excs, pvs, soms, vips = [0, 4, 7, 10], [1, 5, 8, 11], [2, 6, 9, 12], [3, 3, 3, 3]
        for p in inputs:
            path_str = os.path.dirname(p)
            params_list = np.array(path_str.split('/')[-1].split('_')).astype(float)
            fr_arr, ai_arr = np.full((4, 2), np.nan), np.full((4, 2), np.nan)
            if os.path.isfile(os.path.join(path_str, 'fr.dat')):
                fr_arr = np.loadtxt(os.path.join(path_str, 'fr.dat'), delimiter=',')
            if os.path.isfile(os.path.join(path_str, 'ai.dat')):
                ai_arr = np.loadtxt(os.path.join(path_str, 'ai.dat'), delimiter=',')
            for i, lyr in enumerate(self.lyrs):
                # include items in the dataframe
                tmp_dict = {
                    self.dims['x']: params_list[0],
                    self.dims['y']: params_list[1],
                    self.dims['za']: params_list[2],
                    self.dims['zb']: params_list[3],
                    'lyr': lyr,
                    r'$r_{Exc}$': fr_arr[excs[i], 0],
                    r'$r_{PV}$': fr_arr[pvs[i], 0],
                    r'$r_{SOM}$': fr_arr[soms[i], 0],
                    r'$r_{VIP}$': fr_arr[vips[i], 0],
                    'pairwise\ncorrelation': ai_arr[i, 0],
                    'CV(ISI)': ai_arr[i, 1]
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
        print('levels =')
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

    # make data matrix
    def make_data(self):
        # structure
        for zb in self.zb_lvls:
            self.mtxs[str(zb)] = {}
            self.fits[str(zb)] = {}
            self.rmse[str(zb)] = {}
            self.rmse_n[str(zb)] = {}
            for za in self.za_lvls:
                self.mtxs[str(zb)][str(za)] = {}
                self.fits[str(zb)][str(za)] = {}
                self.rmse[str(zb)][str(za)] = np.full((len(self.y_lvls), len(self.x_lvls)), 0)
                self.rmse_n[str(zb)][str(za)] = np.full((len(self.y_lvls), len(self.x_lvls)), 0)
                for plotvar in self.plotvars:
                    self.mtxs[str(zb)][str(za)][plotvar] = np.full((len(self.lyrs), len(self.y_lvls), len(self.x_lvls)), np.nan)
                    if plotvar in self.criteria:
                        self.fits[str(zb)][str(za)][plotvar] = np.full((len(self.lyrs), len(self.y_lvls), len(self.x_lvls)), 0)

        # data
        for i, lyr in enumerate(self.df['lyr'].tolist()):
            x = self.df[self.dims['x']][i]
            y = self.df[self.dims['y']][i]
            za = self.df[self.dims['za']][i]
            zb = self.df[self.dims['zb']][i]
            idx_lyr = self.lyrs.index(lyr)
            idx_x, idx_y = self.x_lvls.index(x), self.y_lvls.index(y)
            fr_exc = self.df[r'$r_{Exc}$'][i]
            fr_pv = self.df[r'$r_{PV}$'][i]
            fr_som = self.df[r'$r_{SOM}$'][i]
            fr_vip = self.df[r'$r_{VIP}$'][i]
            for plotvar in self.plotvars:
                value = self.df[plotvar][i]
                self.mtxs[str(zb)][str(za)][plotvar][idx_lyr, idx_x, idx_y] = value
                # fitness to ground state criteria
                if plotvar in self.criteria:
                    if self.criteria[plotvar][idx_lyr][0] <= value <= self.criteria[plotvar][idx_lyr][1]:
                        self.fits[str(zb)][str(za)][plotvar][idx_lyr, idx_x, idx_y] = 1
            if fr_exc != np.nan:
                self.rmse[str(zb)][str(za)][idx_x, idx_y] += (fr_exc - self.rmse_criteria[r'$r_{Exc}$'][idx_lyr])**2
                self.rmse_n[str(zb)][str(za)][idx_x, idx_y] += 1
            if fr_pv != np.nan:
                self.rmse[str(zb)][str(za)][idx_x, idx_y] += (fr_pv - self.rmse_criteria[r'$r_{PV}$'][idx_lyr])**2
                self.rmse_n[str(zb)][str(za)][idx_x, idx_y] += 1
            if fr_som != np.nan:
                self.rmse[str(zb)][str(za)][idx_x, idx_y] += (fr_som - self.rmse_criteria[r'$r_{SOM}$'][idx_lyr])**2
                self.rmse_n[str(zb)][str(za)][idx_x, idx_y] += 1
            if idx_lyr == 0 and fr_vip != np.nan:
                self.rmse[str(zb)][str(za)][idx_x, idx_y] += (fr_vip - self.rmse_criteria[r'$r_{VIP}$'][idx_lyr])**2
                self.rmse_n[str(zb)][str(za)][idx_x, idx_y] += 1

        # calculate RMSEs
        for k1, v1 in self.rmse.items():
            for k2, v2 in v1.items():
                self.rmse[k1][k2] = np.sqrt(np.divide(v2, self.rmse_n[k1][k2]))
                print('za:{}, zb:{}, rmse:'.format(k2, k1))
                print(self.rmse[k1][k2])
                print(self.rmse_n[k1][k2])

    def make_plots(self, afx=None):
        # plot
        for zb in self.zb_lvls:
            for za in self.za_lvls:
                self.colormap(za, zb, afx=afx)

    def colormap(self, za, zb, afx=None):
        # set plotting variables
        fig, axs = plt.subplots(4, len(self.plotvars), figsize=self.figsize, sharex=True, sharey=True)
        xs, ys = np.array(self.x_lvls), np.array(self.y_lvls)
        plt.xlim((xs[0], xs[-1]))
        plt.ylim((ys[0], ys[-1]))
        xlbl, ylbl = self.dims['x'], self.dims['y']
        ylbl = ylbl.replace('bg_rate', r'$r_{bg}$')
        # loop variables to plot
        for c, plotvar in enumerate(self.plotvars):
            vmin, vmax = 0.0, self.vmaxs[plotvar]
            if plotvar == 'pairwise\ncorrelation':
                vmin = -vmax
            # fit in all criteria and layers
            all_fit = np.ones((len(ys), len(xs)))
            # loop layers
            for r in range(4):
                ax = axs[r, c]
                # plot data
                data = self.mtxs[str(zb)][str(za)][plotvar][r].T    # (y, x)
                if plotvar == r'$r_{VIP}$' and r > 0:
                    ax.axis('off')
                    data = np.full(data.shape, np.nan)

                # simple grid (for colorbar)
                im = ax.imshow(data, interpolation='none',
                    cmap=self.cmaps[plotvar],
                    origin='lower', extent=self.extent,
                    vmin=vmin, vmax=vmax, zorder=1)

                # patch to cover grid (shitty)
                rect = patches.Rectangle((xs[0],ys[0]),(xs[-1]-xs[0]),(ys[-1]-ys[0]), edgecolor='w',facecolor='w', zorder=2)
                ax.add_patch(rect)

                # contour
                cf = ax.contourf(data,
                    levels=np.linspace(vmin, vmax, 11),
                    cmap=self.cmaps[plotvar],
                    origin='lower', extent=self.extent,
                    vmin=vmin, vmax=vmax, zorder=3,
                    extend='max')
                cf.cmap.set_over('darkblue')
                ct = ax.contour(data,
                    levels=np.linspace(vmin, vmax, 11),
                    origin='lower', extent=self.extent,
                    colors='gray', linewidths=0.5,
                    vmin=vmin, vmax=vmax, zorder=4)
                if self.mark_flg:
                    ax.clabel(ct, fmt=self.clabel_format[plotvar],
                    colors='k', inline=True, fontsize=10)

                # single- & triple-fit patches
                if plotvar in self.criteria:
                    # single fit
                    fits = self.fits[str(zb)][str(za)][plotvar][r].T
                    # triple-fit
                    tri_fits = np.ones(data.shape)
                    for k in self.criteria.keys():
                        tri_fits = np.multiply(tri_fits, self.fits[str(zb)][str(za)][k][r].T) # (y, x)
                        # save all-fit (across layers), for RMSE later
                        if c == 0:
                            for row in range(4):
                                all_fit = np.multiply(all_fit, self.fits[str(zb)][str(za)][k][row].T)
                    # matrix for single and triple-fit
                    fit_mtx = np.zeros(data.shape)
                    fit_mtx[np.where(fits == 1)] = 10.0 # single
                    fit_mtx[np.where(tri_fits==1)] = 20.0 # triple
                    # validate
                    print('{}, layer {}, criteria = {}'.format(plotvar, self.lyrs[r], self.criteria[plotvar][r]))
                    print('data:\n{}\n{}'.format(data[::-1], fit_mtx[::-1]))
                    # plot
                    cf_fit = ax.contourf(fit_mtx,
                        levels=[9., 19., 25.],
                        origin='lower',
                        extent=self.extent,
                        hatches=['//', '++', ''],
                        alpha=0.0,
                        linewidth=0.5,
                        zorder=5)
                    cf_fit = ax.contour(fit_mtx,
                        levels=[9., 19.],
                        origin='lower',
                        extent=self.extent,
                        linewidth=0.5,
                        colors='k', zorder=6)

                # RMSE
                if self.mark_flg:
                    rmse_mtx = self.rmse[str(zb)][str(za)].T  # (y, x)
                    # mark best RMSEs
                    if c == 0:
                        rmse_mtx = self.rmse[str(zb)][str(za)].T
                        best_rmse = np.min(rmse_mtx)
                        i_y, i_x = np.where(rmse_mtx==best_rmse)
                        ax.scatter(xs[i_x], ys[i_y], s=100, marker='o',
                        color='yellow', edgecolor='k', zorder=7)
                        # best RMSE in the triple-fit area
                        if len(rmse_mtx[np.where(all_fit==1)]) > 0:
                            best_rmse = np.min(rmse_mtx[np.where(all_fit==1)])
                            i_y, i_x = np.where(rmse_mtx==best_rmse)
                            ax.scatter(xs[i_x], ys[i_y], s=100, marker='*',
                            color='yellow', edgecolor='k', zorder=8)
                            ax.text(xs[i_x], ys[i_y],
                            '{:.1f}'.format(best_rmse),
                            color='r', fontsize=10, zorder=9)
                    # text RMSE values
                    if plotvar == r'$r_{PV}$':
                        for a, y in enumerate(ys):
                            for b, x in enumerate(xs):
                                clr = 'gray'
                                if all_fit[a, b] == 1:
                                    clr = 'magenta'
                                ax.text(x, y, '{:.1f}'.format(rmse_mtx[a, b]),
                                color=clr, fontsize=8,
                                horizontalalignment='center',
                                verticalalignment='center',
                                zorder=10)

                # many settings
                ax.set_aspect(float((xs[-1] - xs[0])/(ys[-1] - ys[0])))

                # titles
                if r == 0:
                    ax.set_title(plotvar + '\n ')

                # xlabels and colorbars
                if r == 3:
                    ax.set_xlabel(xlbl)
                    cbar = fig.colorbar(im, ax=axs[:, c], orientation='horizontal', shrink=0.8, aspect=10)

                # ylabel
                if c == 0:
                    ax.set_ylabel(ylbl)
                    ax.text(-0.75, 0.5, self.lyrs[r], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        # if self.mark_flg:
        #     plt.suptitle('minimum ' + r'$RMSE_{r}$' + '={:.2f}'.format(best_rmse))
        plot_name = '{}={},{}={}'.format(self.dims['za'], str(za), self.dims['zb'], str(zb))
        plot_name = os.getcwd().split('/')[-1] + '_' + plot_name
        if isinstance(afx, str):
            plot_name = afx + plot_name
        fig.savefig(plot_name + '.png', bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    inputs = sys.argv[5:]
    dims = sys.argv[1:5]
    # scandata = ScanData(inputs)
    scandata = ScanData(inputs, dims=dims)
    scandata.mark_flg = True
    scandata.make_plots(afx='mark_')
