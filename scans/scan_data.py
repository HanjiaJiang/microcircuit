import os
import sys
import copy
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.interpolate import interp1d, interp2d
matplotlib.rcParams['font.size'] = 25.0
np.set_printoptions(precision=3, linewidth=500, suppress=True)
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes

class ScanData:
    # directory list, dimension dictionary
    def __init__(self, inputs, dims=None, mark=False, star_loc=None, xybounds=None):
        self.set_plotvars()
        self.set_params(mark, star_loc)
        self.set_dataframe(inputs, dims)
        self.set_plot(dims, xybounds)

    def set_plotvars(self, plotvars=None, figsize=(16, 12)):
        self.plotvars = np.array([r'$r_{Exc}$', 'corr', 'CV(ISI)',
            r'$r_{PV}$', r'$r_{SOM}$', r'$r_{VIP}$'])
        self.figsize = figsize
        if isinstance(plotvars, list):
            self.plotvars = self.plotvars[np.array(plotvars)]
            self.figsize = (len(plotvars)*3, 12)

    def set_params(self, mark, star_loc):
        self.layer_labels = ['L2/3', 'L4', 'L5', 'L6']
        self.mtxs, self.fits, self.rmse, self.rmse_n = {}, {}, {}, {}
        # ground state criteria (Maksimov et al., 2018)
        self.criteria = {
            r'$r_{Exc}$': np.tile([0.0, 10.0], (4, 1)),
            'corr': np.tile([0.0001, 0.008], (4, 1)),
            'CV(ISI)': np.tile([0.76, 1.2],(4, 1))}
        self.vlims = {
            r'$r_{Exc}$': [0., 10.],
            'corr': [-0.05, 0.05],
            'CV(ISI)': [0.5, 1.5],
            r'$r_{PV}$': [0., 100],
            r'$r_{SOM}$': [0., 100],
            r'$r_{VIP}$': [0., 10]}
        # RMSE
        self.mark = mark
        self.star_loc = star_loc
        self.rmse_criteria = {
            r'$r_{Exc}$': [2.7, 0.5, 6.8, 6.1],
            r'$r_{PV}$': [13.8, 10.2, 7.5, 16.9],
            r'$r_{SOM}$': [2.6, 2.6, 2.8, 3.9],
            r'$r_{VIP}$': [14.6]}
        self.cmaps = {
            r'$r_{Exc}$': 'Blues',
            'corr': 'RdBu',
            'CV(ISI)': 'Blues',
            r'$r_{PV}$': 'Blues',
            r'$r_{SOM}$': 'Blues',
            r'$r_{VIP}$': 'Blues'}
        self.clabel_format = {
            r'$r_{Exc}$': '%1.0f',
            'corr': '%1.3f',
            'CV(ISI)': '%1.2f',
            r'$r_{PV}$': '%1.0f',
            r'$r_{SOM}$': '%1.0f',
            r'$r_{VIP}$': '%1.0f'}

    def set_dataframe(self, inputs, dims):
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
            fr_arr = np.loadtxt(os.path.join(path_str, 'fr.dat'), delimiter=',')
            ai_arr = np.loadtxt(os.path.join(path_str, 'ai.dat'), delimiter=',')
            if fr_arr.shape != (13, 2) or ai_arr.shape != (4, 2):
                print('{} failed'.format(path_str))
            for i, lyr in enumerate(self.layer_labels):
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
                    'corr': ai_arr[i, 0],
                    'CV(ISI)': ai_arr[i, 1]
                    }
                data_list.append(tmp_dict)
        self.df = pd.DataFrame(data_list)

    # set plot data: dimensions, coordinates, matrix
    def set_plot(self, dims, xybounds):
        self.update_dims(dims)
        xname, yname, zaname, zbname = self.dims['x'], self.dims['y'], self.dims['za'], self.dims['zb']
        xs, ys, zas, zbs = self.df[xname].tolist(), self.df[yname].tolist(), self.df[zaname].tolist(), self.df[zbname].tolist() # float
        lvls_s = self.get_lvls([xs, ys, zas, zbs])
        self.x_lvls, self.y_lvls, self.za_lvls, self.zb_lvls = lvls_s[0], lvls_s[1], lvls_s[2], lvls_s[3]
        print('levels =')
        print(self.x_lvls, self.y_lvls, self.za_lvls, self.zb_lvls)
        if isinstance(xybounds, list) and len(xybounds) == 4:
            self.xybounds = xybounds
        else:
            self.xybounds = [self.x_lvls[0], self.x_lvls[-1], self.y_lvls[0], self.y_lvls[-1]]
        self.make_data()
        # self.make_plots()

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
                    self.mtxs[str(zb)][str(za)][plotvar] = np.full((len(self.layer_labels), len(self.y_lvls), len(self.x_lvls)), np.nan)
                    if plotvar in self.criteria:
                        self.fits[str(zb)][str(za)][plotvar] = np.full((len(self.layer_labels), len(self.y_lvls), len(self.x_lvls)), 0)

        # data
        for i, lyr in enumerate(self.df['lyr'].tolist()):
            x = self.df[self.dims['x']][i]
            y = self.df[self.dims['y']][i]
            za = self.df[self.dims['za']][i]
            zb = self.df[self.dims['zb']][i]
            idx_lyr = self.layer_labels.index(lyr)
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

    def make_plots(self, afx=None, plotvars=None):
        self.set_plotvars(plotvars)
        # plot
        for zb in self.zb_lvls:
            for za in self.za_lvls:
                self.colormap(za, zb, afx=afx)

    def interpolate_frame(self, frame, resolution=(100,100), kind='linear'):
        dim_x, dim_y = frame.shape
        if np.iscomplexobj(frame):
            re_frame = interpolate_frame(np.real(frame),
                                         resolution=resolution,
                                         kind=kind)
            im_frame = interpolate_frame(np.imag(frame),
                                         resolution=resolution,
                                         kind=kind)
            return re_frame + 1j*im_frame
        else:
            f = interp2d(np.arange(dim_x), np.arange(dim_y), np.nan_to_num(frame), kind=kind)
            frame = f(np.linspace(0, dim_x-1, resolution[0]),
                      np.linspace(0, dim_y-1, resolution[1]))
        return frame

    def interpol(self, data, resolution=(100, 100), cut=0.5):
        data_interpol = copy.deepcopy(data) # to be interpolated
        data_mask = copy.deepcopy(data) # mask for nan values
        data_mask[np.isnan(data_mask) == False] = 1.
        data_mask[np.isnan(data_mask)] = 0.
        for i, row in enumerate(data_interpol):
            for j, item in enumerate(row):
                if np.isnan(item):
                    a = data[i, j+1]  if j < len(row)-1 else np.nan
                    b = data[i+1, j]  if i < len(data_interpol)-1 else np.nan
                    c = data[i, j-1]  if j > 0 else np.nan
                    d = data[i-1, j]  if i > 0 else np.nan
                    arr = np.array([a, b, c, d])
                    data_interpol[i, j] = np.nanmean(arr)
        data_interpol = self.interpolate_frame(data_interpol, resolution=resolution)
        data_mask = self.interpolate_frame(data_mask, resolution=resolution)
        data_interpol[data_mask<cut] = np.nan
        return data_interpol

    def get_multifit(self, za, zb, plotvar, r):
        fits = self.fits[str(zb)][str(za)][plotvar][r].T
        multi_fits = np.ones((len(self.y_lvls), len(self.x_lvls)))
        all_fits = np.ones((len(self.y_lvls), len(self.x_lvls)))
        for k in self.criteria.keys():
            # multi-fit (same layer, across criteria)
            multi_fits = np.multiply(multi_fits, self.fits[str(zb)][str(za)][k][r].T) # (y, x)
            # all-fit (across layers and criteria)
            for row in range(4):
                all_fits = np.multiply(all_fits, self.fits[str(zb)][str(za)][k][row].T)
        return fits, multi_fits, all_fits

    def plot_fitpatch(self, ax, fits_easy, fits_hard, extent):
        fit_mtx = np.zeros(fits_easy.shape)
        fit_mtx[np.where(fits_easy==1)] = 10.0
        fit_mtx[np.where(fits_hard==1)] = 20.0
        cf_fit = ax.contourf(fit_mtx,
            levels=[9., 19., 25.],
            origin='lower',
            extent=extent,
            hatches=['//', '++', ''],
            alpha=0.0,
            linewidth=0.25,
            zorder=5)
        cf_fit = ax.contour(fit_mtx,
            levels=[9., 19.],
            origin='lower',
            extent=extent,
            linewidth=0.25,
            zorder=6,
            colors='k')

    def plot_fitpoint(self, ax, fits_gray, fits_black, fits_green=None):
        iy, ix = np.where(fits_gray == 1)
        xs, ys = np.array(self.x_lvls)[ix], np.array(self.y_lvls)[iy]
        ax.scatter(xs, ys, color='gray', s=8, marker='o', zorder=5)
        iy, ix = np.where(fits_black == 1)
        xs, ys = np.array(self.x_lvls)[ix], np.array(self.y_lvls)[iy]
        ax.scatter(xs, ys, color='black', s=8, marker='o', zorder=6)
        if fits_green is not None:
            iy, ix = np.where(fits_green == 1)
            xs, ys = np.array(self.x_lvls)[ix], np.array(self.y_lvls)[iy]
            ax.scatter(xs, ys, color='green', s=8, marker='o', zorder=7)
        pass

    # not using
    def plot_contour(self, ax, data, extent, plotvar, vmin, vmax):
        # contour
        datamin, datamax = np.nanmin(data), np.nanmax(data)
        if datamin == datamax or np.isnan(datamin) or np.isnan(datamax):
            pass
        else:
            cf = ax.contourf(data,
                levels=np.linspace(np.nanmin(data), np.nanmax(data), 11),
                cmap=self.cmaps[plotvar],
                origin='lower', extent=extent,
                vmin=vmin, vmax=vmax, zorder=3,
                extend='max')
            cf.cmap.set_over('black')
        ct = ax.contour(data,
            levels=np.linspace(np.nanmin(data), np.nanmax(data), 11),
            origin='lower', extent=extent,
            colors='gray', linewidths=0.5,
            vmin=vmin, vmax=vmax, zorder=4)
        if self.mark:
            ax.clabel(ct, fmt=self.clabel_format[plotvar],
            colors='k', inline=True, fontsize=10)

    def contour_ai(self, ax, data, extent, plotvar, vmin, vmax):
        datamin, datamax = np.nanmin(data), np.nanmax(data)
        if datamin == datamax or np.isnan(datamin) or np.isnan(datamax):
            pass
        else:
            ct = ax.contour(data,
                levels=np.array([self.criteria[plotvar][0][0], self.criteria[plotvar][0][1]]),
                origin='lower', extent=extent,
                colors='k', linewidths=0.5, zorder=4)
            ax.clabel(ct, fmt=self.clabel_format[plotvar],
            colors='k', inline=True, fontsize=10)

    def mark_rmse(self, rmse_mtx, xs, ys, all_fits):
        for a, y in enumerate(ys):
            for b, x in enumerate(xs):
                clr = 'gray'
                if all_fits[a, b] == 1:
                    clr = 'magenta'
                ax.text(x, y, '{:.1f}'.format(rmse_mtx[a, b]),
                color=clr, fontsize=8,
                horizontalalignment='center',
                verticalalignment='center',
                zorder=10)

    def colormap(self, za, zb, afx=None):
        # set plotting variables
        fig, axs = plt.subplots(4, len(self.plotvars), figsize=self.figsize, sharey=True)
        xs, ys = np.array(self.x_lvls), np.array(self.y_lvls)
        extent1 = [xs[0], xs[-1], ys[0], ys[-1]]
        plt.setp(axs, xticks=xs[::2], yticks=ys[::2])
        xlbl, ylbl = self.dims['x'], self.dims['y']
        ylbl = ylbl.replace('bg', r'$r_{bg}$')
        # add invisible frame for x and y labels
        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        plt.xlabel(xlbl)
        plt.ylabel(ylbl)
        plt.gca().xaxis.set_label_coords(.5, 0.25)
        plt.gca().yaxis.set_label_coords(-0.01, 0.65)
        # loop variables to plot
        for c, plotvar in enumerate(self.plotvars):
            vmin, vmax = self.vlims[plotvar][0], self.vlims[plotvar][1]
            # loop layers
            for r in range(4):
                ax = axs[r, c]
                # plot interpolated data (y, x)
                data = self.interpol(self.mtxs[str(zb)][str(za)][plotvar][r].T)
                if plotvar == r'$r_{VIP}$' and r > 0:
                    ax.axis('off')
                    data = np.full(data.shape, np.nan)

                current_cmap = matplotlib.cm.get_cmap()
                current_cmap.set_bad(color='gray')

                # simple grid (for colorbar)
                im = ax.imshow(data, interpolation='none',
                    cmap=self.cmaps[plotvar],
                    origin='lower', extent=extent1,
                    vmin=vmin, vmax=vmax, zorder=1)
                im.cmap.set_over('midnightblue')
                if plotvar in ['corr', 'CV(ISI)']:
                    self.contour_ai(ax, data, extent1, plotvar, vmin, vmax)

                # single- & triple-fit (not interpolated)
                if c == 0:
                    fits, tri_fits, all_fits = self.get_multifit(za, zb, plotvar, r)
                    self.plot_fitpoint(ax, tri_fits, all_fits)
                    # self.plot_fitpatch(ax, fits, tri_fits, extent2)

                # RMSE
                if self.mark:
                    rmse_mtx = self.rmse[str(zb)][str(za)].T  # (y, x)
                    # mark best RMSEs in the triple-fit area
                    if c == 0:
                        rmse_mtx = self.rmse[str(zb)][str(za)].T
                        if len(rmse_mtx[np.where(all_fits==1)]) > 0:
                            best_rmse = np.min(rmse_mtx[np.where(all_fits==1)])
                            i_y, i_x = np.where(rmse_mtx==best_rmse)
                            ax.scatter(xs[i_x], ys[i_y], s=200, marker='*',
                            facecolor='none', edgecolor='k', zorder=8)
                            # best RMSE point
                            ax.text(xs[i_x], ys[i_y],
                            '{:.1f}'.format(best_rmse),
                            color='r', fontsize=10, zorder=9)
                    # text RMSE values
                    # if plotvar == r'$r_{PV}$':
                    #     self.mark_rmse(rmse_mtx, xs, ys, all_fits)

                # mark defined star
                if isinstance(self.star_loc, list):
                    if c==0:
                        ax.scatter(self.star_loc[0], self.star_loc[1],
                        s=150, marker='*', facecolor='none', edgecolor='k', zorder=8)

                # many settings
                ax.set_xlim(self.xybounds[0], self.xybounds[1])
                ax.set_ylim(self.xybounds[2], self.xybounds[3])
                ax.set_aspect(float((self.xybounds[1] - self.xybounds[0])/(self.xybounds[3] - self.xybounds[2])))
                ax.tick_params(axis='both', labelsize=matplotlib.rcParams['font.size']*.7)

                # titles
                if r == 0:
                    ax.set_title(plotvar, y=1.0, pad=50)
                    # xlabels and colorbars
                    # ax.set_xlabel(xlbl)
                    cbar = fig.colorbar(im, ax=axs[:, c], location='bottom', shrink=0.8, aspect=10)
                    cbar.ax.tick_params(labelsize=matplotlib.rcParams['font.size']*.8)

                # if plotvar == r'$r_{VIP}$' and r == 0:
                #     ax.set_xlabel(xlbl)

                if r < 3 and plotvar != r'$r_{VIP}$':
                    plt.setp(ax.get_xticklabels(), visible=False)

                # ylabel
                if c == 0:
                    # ax.set_ylabel(ylbl)
                    ax.text(-0.8, 0.5, self.layer_labels[r], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        plot_name = '{}={},{}={}'.format(self.dims['za'], str(za), self.dims['zb'], str(zb))
        plot_name = os.getcwd().split('/')[-1] + '_' + plot_name
        if isinstance(afx, str):
            plot_name = afx + '_' + plot_name
        fig.savefig(plot_name + '.png', bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    inputs = sys.argv[5:]
    dims = sys.argv[1:5]
    scandata = ScanData(inputs, dims=dims, xybounds=[3.8, 10.2, 3.8, 7.2])
    # scandata = ScanData(inputs, dims=dims, star_loc=[8., 4.5], xybounds=[3.8, 10.2, 3.8, 7.2])
    scandata.make_plots(afx='all')
    scandata.make_plots(afx='012', plotvars=[0, 1, 2])
    scandata.make_plots(afx='0345', plotvars=[0, 3, 4, 5])
