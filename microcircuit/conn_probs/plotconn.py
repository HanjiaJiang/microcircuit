import os
import numpy as np
import matplotlib.pyplot as plt

def plotconn(fn, raw=False, threshold=False, vipconn=True):
    # set labels
    # lyr_lbls = ['L2/3', 'L4', 'L5', 'L6']
    # lyr_pos = [1.5/13, 4.5/13, 7.5/13, 11/13]
    pop_lbls = ['L2/3 Exc', 'L2/3 PV', 'L2/3 SOM', 'L2/3 VIP',
                'L4 Exc', 'L4 PV', 'L4 SOM',
                'L5 Exc', 'L5 PV', 'L5 SOM',
                'L6 Exc', 'L6 PV', 'L6 SOM']
    # pop_lbls = ['Exc', 'PV', 'SOM', 'VIP',
    #             'Exc', 'PV', 'SOM',
    #             'Exc', 'PV', 'SOM',
    #             'Exc', 'PV', 'SOM']

    # get data
    data = np.loadtxt(fn, delimiter=',')
    if vipconn is True:
        data[[6, 9, 12], 3] = data[2, 3]

    # set plot parameters
    # plt.rcParams['figure.constrained_layout.use'] = True
    plt.rcParams['font.size'] = 15
    plt.rcParams['figure.figsize'] = 10, 10
    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.bottom'] = False
    plt.rcParams['xtick.labeltop'] = True
    plt.rcParams['xtick.labelbottom'] = False

    # color map
    cx = plt.imshow(data,
    interpolation='none',
    cmap='Blues',
    extent=[0, 13, 0, 13],
    vmin=0.0,
    vmax=1.0
    )

    # values
    mtx = conn_estimate_mtx()
    for i in range(13):
        for j in range(13):
            prob, fsize, fcolor = data[i, j], 10, 'k'
            text = '{:.2f}'.format(prob)
            if threshold is True and prob < 0.05:
                fcolor = 'gray'
            if '7-15' in fn and mtx[i][j] == 1:
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

def conn_estimate_mtx():
    mtx = [[0,1,0,0,0,0,0,0,0,0,0,0,0],
           [0,0,0,0,0,0,0,0,0,0,0,0,0],
           [0,0,0,0,0,0,0,0,0,0,0,0,0],
           [0,0,0,0,0,0,0,0,0,0,0,0,0],
           [0,0,0,0,0,0,0,0,0,0,0,0,0],
           [0,0,0,0,1,0,0,0,0,0,0,0,0],
           [0,0,0,0,0,0,0,0,0,0,0,0,0],
           [0,0,0,0,0,0,0,0,1,1,0,0,0],
           [0,0,0,0,0,0,0,1,1,1,0,0,0],
           [0,0,0,0,0,0,0,1,1,1,0,0,0],
           [0,0,0,0,0,0,0,0,0,0,0,1,1],
           [0,0,0,0,0,0,0,0,0,0,1,1,1],
           [0,0,0,0,0,0,0,0,0,0,1,1,1],
           ]
    return mtx

if __name__ == '__main__':
    fns = os.listdir()
    for fn in fns:
        if fn.startswith('conn') and fn.endswith('.csv'):
            plotconn(fn, threshold=True, vipconn=False)
