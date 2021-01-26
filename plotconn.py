import os
import copy
import numpy as np
import matplotlib.pyplot as plt
import microcircuit.network as network
import microcircuit.conn as conn

class PlotConn():
    def __init__(self, path):
        self.setup(path)

    def setup(self, path):
        self.network = network
        self.probs_path = path
        self.parameters()
        self.set_plt()

    def set_plt(self):
        plt.rcParams['font.size'] = 15
        plt.rcParams['figure.figsize'] = 10, 10
        plt.rcParams['xtick.top'] = True
        plt.rcParams['xtick.bottom'] = False
        plt.rcParams['xtick.labeltop'] = True
        plt.rcParams['xtick.labelbottom'] = False

    def run_plots(self):
        fns = os.listdir(self.probs_path)
        for fn in fns:
            if fn.startswith('conn') and fn.endswith('.csv'):
                plotconn.plot_by_subtype(os.path.join(self.probs_path, fn))
            if fn.startswith('raw') and fn.endswith('.csv'):
                plotconn.plot_by_subtype(os.path.join(self.probs_path, fn), dtype='raw')

    def parameters(self):
        self.pops = ['Exc', 'PV', 'SOM', 'VIP',
                    'Exc', 'PV', 'SOM',
                    'Exc', 'PV', 'SOM',
                    'Exc', 'PV', 'SOM']

        self.layers = ['L2/3', 'L4', 'L5', 'L6']
        self.positions = [2/13, 5.5/13, 8.5/13, 11.5/13]
        self.bbp_data = np.loadtxt(os.path.join(self.probs_path, 'raw_bbp.csv'), delimiter=',')

    def plot_by_subtype(self, fn, vmin=-2.0, vmax=2.0, cmap='RdBu', fontsizes=(12, 15), dtype='conn', thr=None):
        # get data
        data = np.loadtxt(fn, delimiter=',')
        # need to import bbp data, if it is raw data
        if dtype == 'raw':
            # raw = copy.deepcopy(data)
            for i, row in enumerate(data):
                for j, item in enumerate(row):
                    if data[i, j] == 0.:
                        data[i, j] = self.bbp_data[i, j]
        # else:
        #     raw_fn = os.path.join(os.path.dirname(fn), os.path.basename(fn).replace('conn', 'raw'))
        #     raw = np.loadtxt(raw_fn, delimiter=',') if os.path.isfile(raw_fn) else np.ones(data.shape)

        # color map
        data[:, [1,2,3,5,6,8,9,11,12]] *= -1.
        cx = plt.imshow(data, interpolation='none', cmap=cmap,
            extent=[0, 13, 0, 13], vmin=vmin, vmax=vmax)

        # values
        flg_fn = os.path.join(os.path.dirname(fn), os.path.basename(fn).replace('raw', 'flg').replace('conn', 'flg'))
        flg_mtx = np.loadtxt(flg_fn, delimiter=',')
        for i in range(13):
            for j in range(13):
                prob = np.abs(data[i, j])
                if thr is not None and prob < thr:
                    continue
                text, fcolor = '{:.0f}%'.format(prob*100), 'k'
                weight = 'heavy' if flg_mtx[i][j] == 1 else 'normal'
                fsize = fontsizes[1] if flg_mtx[i][j] == 1 else fontsizes[0]
                if flg_mtx[i][j] == 2:
                    text += '*'
                plt.text(j+0.5, 12-i+0.5, text,
                    fontsize=fsize, color=fcolor, horizontalalignment='center',
                    verticalalignment='center', weight=weight)
        ax = plt.gca()

        # source/target labels
        plt.text(6.5/13, 1.15, 'presynaptic', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20)
        plt.text(-0.15, 6.5/13, 'postsynaptic', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, rotation='90', fontsize=20)
        for i in range(4):
            plt.text(self.positions[i], 1.1, self.layers[i], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20)
            plt.text(-0.1, 1-self.positions[i], self.layers[i], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20)

        # population labels
        plt.xticks(np.arange(0.5, 13, 1.0), self.pops)
        plt.yticks(np.arange(0.5, 13, 1.0), self.pops[::-1])

        # others
        plt.vlines([4, 7, 10], 0, 13)
        plt.hlines([3, 6, 9], 0, 13)
        plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")
        plt.setp(ax.get_yticklabels(), ha="right")
        # plt.colorbar(orientation='vertical', shrink=0.5)
        # plt.title('connection probabilities')

        # save
        plt.savefig(fn.replace('.csv', '.png'), bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    probs_path = './microcircuit/conn_probs/'
    plotconn = PlotConn(probs_path)
    plotconn.run_plots()
