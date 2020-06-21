import os
import numpy as np
import matplotlib.pyplot as plt

def plotconn(fn, raw=False):
    pop_lbls = ['L2/3_Exc', 'L2/3_PV', 'L2/3_SOM', 'L2/3_VIP', 'L4_Exc', 'L4_PV', 'L4_SOM', 'L5_Exc', 'L5_PV', 'L5_SOM', 'L6_Exc', 'L6_PV', 'L6_SOM']
    data = np.loadtxt(fn, delimiter=',')
    data[[6, 9, 12], 3] = data[2, 3]
    plt.rcParams['figure.constrained_layout.use'] = True
    plt.rcParams['font.size'] = 15
    plt.rcParams['figure.figsize'] = 8, 8
    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.bottom'] = False
    plt.rcParams['xtick.labeltop'] = True
    plt.rcParams['xtick.labelbottom'] = False
    cx = plt.imshow(data,
    interpolation='none',
    cmap='Reds',
    extent=[0, 13, 0, 13],
    vmin=0.0,
    vmax=1.0
    )
    for i in range(13):
        for j in range(13):
            plt.text(j+0.5, 12-i+0.5, '{:.2f}'.format(data[i, j]), 
            fontsize=10, horizontalalignment='center', verticalalignment='center')
    plt.vlines([4, 7, 10], 0, 13)
    plt.hlines([3, 6, 9], 0, 13)
    plt.xticks(np.arange(0.5, 13, 1.0), pop_lbls, rotation='30')
    plt.yticks(np.arange(0.5, 13, 1.0), pop_lbls[::-1])
    plt.colorbar(orientation='horizontal', shrink=0.5)
    plt.title('connection probabilities')
    plt.savefig(fn.replace('.csv', '.png'))
    plt.close()

if __name__ == '__main__':
    fns = os.listdir()
    for fn in fns:
        if fn.startswith('conn') and fn.endswith('.csv'):
            plotconn(fn)
