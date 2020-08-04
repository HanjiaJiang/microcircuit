import os
import copy
import json
import pickle
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
matplotlib.rcParams['font.size'] = 17.0

# load fitted STPs
stp_fns, stps = [], {'no-stp': {}}
for fn in os.listdir('microcircuit/stp/'):
    if fn.startswith('stp_fitted') and fn.endswith('.pickle'):
        stp_fns.append(fn)
stp_fns = sorted(stp_fns)
for fn in stp_fns:
    with open(os.path.join('microcircuit/stp/', fn), 'rb') as p:
        stps[fn] = pickle.load(p)

'''
Doiron data
'''
doiron_e = {
    'model': 'tsodyks_synapse',
    'U': 0.75,
    'tau_fac': 0.0,
    'tau_rec': 800.0,
}

doiron_pv = {
    'model': 'tsodyks_synapse',
    'U': 0.9,
    'tau_fac': 0.0,
    'tau_rec': 800.0,
}

doiron_e_weak = {
    'model': 'tsodyks_synapse',
    'U': 0.75,
    'tau_fac': 0.0,
    'tau_rec': 100.0,
}

doiron_pv_weak = {
    'model': 'tsodyks_synapse',
    'U': 0.9,
    'tau_fac': 0.0,
    'tau_rec': 100.0,
}

doiron_e2som = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_fac': 200.0,
    'tau_rec': 0.01,
}

stps['doiron'] = {
    'Exc': {
        'Exc': doiron_e,
        'PV': doiron_e,
        'SOM': doiron_e2som
    },
    'PV': {
        'Exc': doiron_pv,
        'PV': doiron_pv,
        'SOM': doiron_pv,
        'VIP': doiron_pv
    }
}

stps['doiron-w'] = {
    'Exc': {
        'Exc': doiron_e_weak,
        'PV': doiron_e_weak,
        'SOM': doiron_e2som
    },
    'PV': {
        'Exc': doiron_pv_weak,
        'PV': doiron_pv_weak,
        'SOM': doiron_pv_weak,
        'VIP': doiron_pv_weak
    }
}

'''
BBP data
'''
e1 = {
    'model': 'tsodyks_synapse',
    'U': 0.09,
    'tau_rec': 138.0,
    'tau_fac': 670.0
}
e2 = {
    'model': 'tsodyks_synapse',
    'U': 0.5,
    'tau_rec': 671.0,
    'tau_fac': 17.0
}
e3 = {
    'model': 'tsodyks_synapse',
    'U': 0.29,
    'tau_rec': 329.0,
    'tau_fac': 326.0
}
i1 = {
    'model': 'tsodyks_synapse',
    'U': 0.016,
    'tau_rec': 45.0,
    'tau_fac': 376.0
}
i2 = {
    'model': 'tsodyks_synapse',
    'U': 0.25,
    'tau_rec': 706.0,
    'tau_fac': 21.0
}
i3 = {
    'model': 'tsodyks_synapse',
    'U': 0.32,
    'tau_rec': 144.0,
    'tau_fac': 62.0
}
stps['bbp'] = {
    'Exc': {
        'Exc': e2,
        'PV': e2,
        'SOM': e1,
        'VIP': e1
    },
    'PV': {
        'Exc': i3
    },
    'SOM': {
        'Exc': i2
    },
    'VIP': {
        'Exc': i2
    },
}

print('stps = {}'.format(list(stps)))

def func1(x):
    return 1000.*x

def func2(x):
    return x/1000.

# def get_paper(stp_fn):
#     data, json_fn = {}, os.path.join('microcircuit/stp', stp_fn.replace('.pickle', '.json'))
#     if os.path.isfile(json_fn):
#         with open(json_fn, 'r') as js:
#             data = js.load(js)
#     return data

if __name__ == "__main__":
    colors = ['#66CCEE', '#CCBB44', '#AA3377',]
    legends = ['U', 'F', 'D']
    EIset = ['Exc', 'PV', 'SOM', 'VIP']
    EEset = ['L23_Exc', 'L4_Exc', 'L5_Exc', 'L6_Exc']
    # for marking paper
    papers = [['Kapfer, 2007', 'Ma, 2012', 'Karnani, 2016'], ['Lefort, 2017']]
    paper_clrs = [['b', 'r', 'g'], ['magenta']]
    paper_rects = [[[[0.6, 2.4]], [[3.6, 6.4], [7.6, 9.4]], [[6.6, 7.4], [10.6, 11.4], [13.6, 15.4]]], [[[-0.4, 15.4]]]]
    for stp_fn, stp in stps.items():
        EIshift = False
        # E to I
        names_ei, Us_ei, Fs_ei, Ds_ei = [], [], [], []
        for i, pre in enumerate(EIset):
            for j, post in enumerate(EIset):
                name, U, F, D = pre + ' to ' + post, 0., 0., 0.
                try:
                    tmp = stp[pre][post]
                    U, F, D = tmp['U'], tmp['tau_fac']/1000., tmp['tau_rec']/1000.
                except KeyError:
                    # get rid of Exc-to-Exc if layer-specific
                    if pre == 'Exc' and post == 'Exc':
                        EIshift = True
                        continue
                names_ei.append(name.replace('_', '').replace('23', '2/3'))
                Us_ei.append(U)
                Fs_ei.append(F)
                Ds_ei.append(D)

        # E to E
        names_ee, Us_ee, Fs_ee, Ds_ee = [], [], [], []
        for i, pre in enumerate(EEset):
            for j, post in enumerate(EEset):
                name, U, F, D = pre + ' to ' + post, 0., 0., 0.
                try:
                    tmp = stp[pre][post]
                    U, F, D = tmp['U'], tmp['tau_fac']/1000., tmp['tau_rec']/1000.
                except KeyError:
                    pass
                names_ee.append(name.replace('_', '').replace('23', '2/3'))
                Us_ee.append(U)
                Fs_ee.append(F)
                Ds_ee.append(D)

        # combine
        if EIshift is True:
            xs = [np.arange(len(names_ei)) + 1, np.arange(len(names_ee))]
        else:
            xs = [np.arange(len(names_ei)), np.arange(len(names_ee))]
        names = [names_ei, names_ee]
        Us, Fs, Ds, w = [Us_ei, Us_ee], [Fs_ei, Fs_ee], [Ds_ei, Ds_ee], 0.15

        # plot
        fig, ax = plt.subplots(2, 1, figsize=(15, 10), sharey=True)
        for i in range(2):
            for j, (dx, params) in enumerate(zip([-w, 0., w], [Us, Fs, Ds])):
                ax[i].bar(xs[i]+dx, params[i], w, color=colors[j])
            for x in [3.5, 7.5, 11.5]:
                ax[i].axvline(x, 0, 1, color='k')
            ax[i].set_ylabel('probability')
            secax = ax[i].secondary_yaxis('right', functions=(func1, func2))
            secax.set_ylabel(r'$\tau$ (ms)')
            ax[i].set_xticks(xs[i])
            ax[i].set_xticklabels(names[i], rotation='45')
            ax[i].set_ylim(-0.05, 1.05)
            plt.setp(ax[i].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
            if stp_fn == 'stp_fitted_02.pickle':
                for a, paper in enumerate(papers[i]):
                    rect_sets = paper_rects[i][a]
                    x = rect_sets[0][0]
                    ax[i].text(x, 1.12, paper, color=paper_clrs[i][a], fontsize=15)
                    for b, rect_set in enumerate(rect_sets):
                        rect = Rectangle(
                        (rect_set[0],-0.1),
                        rect_set[1] - rect_set[0], 1.2,
                        linewidth=1,edgecolor=paper_clrs[i][a],
                        facecolor='none', linestyle='--',
                        clip_on=False)
                        ax[i].add_patch(rect)

        # mark papers

        # legends
        for i in range(3):
            ax[0].scatter([], [], s=100, marker='s', label=legends[i], color=colors[i])
        ax[0].legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.7))

        fig.tight_layout()
        plt.savefig('stps:{}.png'.format(stp_fn))
        plt.close()
