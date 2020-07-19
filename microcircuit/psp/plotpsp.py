import os
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return sp.ma.masked_array(sp.interp(value, x, y))

def plotpsp(fn, raw=False, vipconn=True, vlimits=(-1.5, 1.5)):
    pop_lbls = ['L2/3 Exc', 'L2/3 PV', 'L2/3 SOM', 'L2/3 VIP',
                'L4 Exc', 'L4 PV', 'L4 SOM',
                'L5 Exc', 'L5 PV', 'L5 SOM',
                'L6 Exc', 'L6 PV', 'L6 SOM']

    # get data
    data = np.loadtxt(fn, delimiter=',')
    if vipconn is True:
        data[[6, 9, 12], 3] = data[2, 3]
    data[:, [1,2,3,5,6,8,9,11,12]] /= 8.

    # set plot parameters
    # plt.rcParams['figure.constrained_layout.use'] = True
    plt.rcParams['font.size'] = 15
    plt.rcParams['figure.figsize'] = 10, 10
    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.bottom'] = False
    plt.rcParams['xtick.labeltop'] = True
    plt.rcParams['xtick.labelbottom'] = False

    # color map
    # norm = MidpointNormalize(vmin=vlimits[0], vmax=vlimits[1], midpoint=0)
    cx = plt.imshow(data,
    interpolation='none',
    cmap='RdBu',
    extent=[0, 13, 0, 13],
    vmin=vlimits[0],
    vmax=vlimits[1]
    # norm=norm
    )

    # values
    flg_mtx = ipsp_estimate_mtx()
    for i in range(13):
        for j in range(13):
            value, fsize, fcolor = data[i, j], 10, 'k'
            text = '{:.2f}'.format(value)
            if np.abs(value) > (vlimits[1] - vlimits[0])/3.:
                fcolor = 'w'                
            if 'ipsp' in fn and flg_mtx[i][j] == 1:
                # fcolor = 'gray'
                text += '*'
            plt.text(j+0.5, 12-i+0.5, text,
            fontsize=fsize,
            color=fcolor,
            horizontalalignment='center',
            verticalalignment='center')
    ax = plt.gca()

    # source/target labels
    # plt.xlabel('presynaptic')
    # plt.ylabel('postsynaptic')
    plt.text(6.5/13, 1.15, 'presynaptic', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=20)
    plt.text(-0.2, 6.5/13, 'postsynaptic', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, rotation='90', fontsize=20)

    # layer labels
    # for i in range(4):
    #     plt.text(-0.14, lyr_pos[i], lyr_lbls[3-i], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    #     plt.text(1 - lyr_pos[i], 1.12, lyr_lbls[3-i], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    # population labels
    plt.xticks(np.arange(0.5, 13, 1.0), pop_lbls)
    plt.yticks(np.arange(0.5, 13, 1.0), pop_lbls[::-1])

    # others
    plt.vlines([4, 7, 10], 0, 13)
    plt.hlines([3, 6, 9], 0, 13)
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")
    plt.setp(ax.get_yticklabels(), ha="right")
    plt.colorbar(orientation='vertical', shrink=0.5)
    # plt.title('connection probabilities')

    # save
    plt.savefig(fn.replace('.csv', '.png'), bbox_inches='tight')
    plt.close()

def ipsp_estimate_mtx():
    mtx = [[0,0,0,1,0,0,0,0,0,0,0,0,0],
           [0,0,0,1,0,0,0,0,0,0,0,0,0],
           [0,0,1,1,0,0,1,0,0,1,0,0,1],
           [0,1,1,1,0,1,1,0,1,1,0,1,1],
           [0,0,0,1,0,0,0,0,0,0,0,0,0],
           [0,0,0,1,0,0,0,0,0,0,0,0,0],
           [0,0,1,1,0,0,1,0,0,1,0,0,1],
           [0,0,0,1,0,0,0,0,0,0,0,0,0],
           [0,0,0,1,0,0,0,0,0,0,0,0,0],
           [0,0,1,1,0,0,1,0,0,1,0,0,1],
           [0,0,0,1,0,0,0,0,0,0,0,0,0],
           [0,0,0,1,0,0,0,0,0,0,0,0,0],
           [0,0,1,1,0,0,1,0,0,1,0,0,1],
           ]
    return mtx

if __name__ == '__main__':
    fns = os.listdir()
    for fn in fns:
        if fn.startswith('psp') and fn.endswith('.csv'):
            plotpsp(fn, vipconn=False)
