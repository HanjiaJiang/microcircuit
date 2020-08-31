import os
import copy
import pickle
import numpy as np
from random import sample
import scipy.stats as stats
from scipy import interpolate

from tkinter import Tk
from tkinter.filedialog import askopenfilename
from multiprocessing import Process
from multiprocessing import Manager

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size'] = 20.0

class Spikes:
    def __init__(self, path, name):
        self.setup(path, name)
        self.load()

    def setup(self, path, name):
        self.path = path
        self.name = name
        self.gids = []
        self.detectors = []
        self.data = []
        self.react_lines = []
        self.veri_dict = {}
        self.fr_result = []
        self.fr_musig = [(2.7, 3.7), (13.8, 8.9), (2.6, 3.6),
                         (14.6, 7.3),
                         (0.5, 0.8), (10.2, 7.2), (2.6, 3.2),
                         (6.8, 5.2), (7.5, 5.2), (2.8, 4.5),
                         (6.1, 6.9), (16.9, 14.3), (3.9, 4.9)]
        self.fr_qrt = [ (0.5, 0.6, 4.5), (7.5, 11.7, 23.3), (0.03, 0.4, 4.1),
                        (8.5, 11.1, 21.0),
                        (0.0, 0.1, 0.7), (4.3, 7.8, 14.7), (0.3, 0.6, 4.9),
                        (2.7, 5.2, 11.2), (4.3, 7.6, 8.7), (0.2, 0.8, 3.6),
                        (0.4, 2.6, 11.5), (4.6, 17.2, 22.0), (0.5, 1.7, 6.9)]
        if os.path.isdir(self.path):
            print('Spikes.__init__(): data directory already exists')
        else:
            os.mkdir(self.path)
            print('Spikes.__init__(): data directory created')
        self.set_labels()

    def read_name(self):
        # Import filenames
        for file in os.listdir(self.path):
            if file.endswith('.gdf') and file.startswith(self.name):
                temp = file.split('-')[0] + '-' + file.split('-')[1]
                if temp not in self.detectors:
                    self.detectors.append(temp)
        # Import GIDs
        gidfile = open(os.path.join(self.path, 'population_GIDs.dat'), 'r')
        for l in gidfile:
            a = l.split()
            self.gids.append([int(a[0]), int(a[1])])
        self.detectors = sorted(self.detectors)
        self.verify_collect('self.detectors={}\n'.format(self.detectors), 'reading')

    def load(self):
        self.read_name()
        if len(self.detectors) > 0 and len(self.gids) > 0:
            for i in list(range(len(self.detectors))):
                all_filenames = os.listdir(self.path)
                thread_filenames = [
                    all_filenames[x] for x in list(range(len(all_filenames)))
                    if all_filenames[x].endswith('gdf') and
                       all_filenames[x].startswith('spike') and
                       (all_filenames[x].split('-')[0]
                        + '-' + all_filenames[x].split('-')[1]) in
                       self.detectors[i]
                ]
                self.verify_collect('thread_filenames:\n', 'reading')
                data_temp = []
                for f in thread_filenames:
                    self.verify_collect('{}\n'.format(f), 'reading')
                    load_tmp = np.array([])
                    try:
                        load_tmp = np.loadtxt(os.path.join(self.path, f))
                    except ValueError:
                        print(os.path.join(self.path, f))
                    else:
                        pass
                    if len(load_tmp.shape) == 2:
                        data_temp.append(load_tmp)
                if len(data_temp) > 0:
                    data = np.concatenate(data_temp)
                    data = data[np.argsort(data[:, 1])]
                    self.data.append(data)
                    self.verify_collect('data_final:\n', 'reading')
                    self.verify_collect('{}...\n'.format(data[:100]), 'reading')
                    # for d in data_final:
                    #     self.verify_collect('{}\n'.format(d), 'reading')
                else:
                    self.data.append([])
                    self.verify_collect('data_final=[]\n', 'reading')
        else:
            print('Spikes.load_spikes_all(): detectors or gids not found')

    def get_data(self, begin, end):
        data_ret = []
        for data in self.data:
            if isinstance(data, np.ndarray):
                data_ret.append(data[(data[:, 1] > begin) & (data[:, 1] <= end)])
            else:
                data_ret.append([])
        return data_ret

    def get_ts_by_id(self, data, id):
        return_data = []
        if isinstance(data, np.ndarray) and data.ndim == 2:
            ids = data[:, 0]
            ts = data[:, 1]
            return_data = ts[ids == id]
        return return_data

    def verify_collect(self, in_str, tag):
        if tag in self.veri_dict.keys():
            self.veri_dict[tag] += in_str
        else:
            self.veri_dict[tag] = in_str

    def verify_print(self, path=None):
        if os.path.isdir(path):
            path_flg = True
        else:
            path_flg = False
        for key, value in self.veri_dict.items():
            if path_flg:
                fpath = os.path.join(path, 'verify-{}.txt'.format(key))
            else:
                fpath = 'verify-{}.txt'.format(key)
            with open(fpath, 'w') as f:
                f.write(value)
                f.close()


    def set_labels(self):
        self.populations = ['L2/3 Exc', 'L2/3 PV', 'L2/3 SOM', 'L2/3 VIP',
                       'L4 Exc', 'L4 PV', 'L4 SOM',
                       'L5 Exc', 'L5 PV', 'L5 SOM',
                       'L6 Exc', 'L6 PV', 'L6 SOM']

        self.subtypes = ['Exc', 'PV', 'SOM', 'VIP', 'Exc', 'PV', 'SOM', 'Exc', 'PV', 'SOM', 'Exc', 'PV', 'SOM']

        self.positions = [0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 0, 1, 2]

        self.layers = [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]

        self.colors = [(68/255,119/255,170/255), (238/255,102/255,119/255), (34/255,136/255,51/255), (204/255,187/255,68/255),
                        (68/255,119/255,170/255), (238/255,102/255,119/255), (34/255,136/255,51/255),
                        (68/255,119/255,170/255), (238/255,102/255,119/255), (34/255,136/255,51/255),
                        (68/255,119/255,170/255), (238/255,102/255,119/255), (34/255,136/255,51/255)]


'''
Calculation
'''
# correlation
def get_corr(list1, list2=None, proc_id=None, rtr_dict=None):
    # multiprocessing
    if proc_id is not None and rtr_dict is not None:
        print('get_corr() proc {} start'.format(proc_id))

    if type(list2) is not list:    # list2 no data, use the same list
        flg_samepop = True
    else:
        flg_samepop = False

    coef_list = []
    for i, hist1 in enumerate(list1):
        if flg_samepop:
            idxs = list(range(i + 1, len(list1)))
        else:
            idxs = list(range(len(list2)))
        for j in idxs:
            if flg_samepop:
                hist2 = list1[j]
            else:
                hist2 = list2[j]
            if np.sum(hist1) != 0 and np.sum(hist2) != 0:
                coef = np.corrcoef(hist1, hist2)[0, 1]
                coef_list.append(coef)
                # if int(proc_id) == 0:
                #     print('hist({})={}\nhist({})={}\n'.format(i, hist1, j, hist2))
            else:
                print('oops, no data in one/both histograms')

    # multiprocessing
    if proc_id is not None and rtr_dict is not None:
        rtr_dict[str(proc_id)] = coef_list
        print('get_corr() proc {} end'.format(proc_id))

    return coef_list


'''
System
'''
# let print() function print to file (use: exec(set2txt))
def set2txt(path):
    str_out = 'import sys\norig_stdout = sys.stdout\nf = open(\'{}\', \'w\')\nsys.stdout = f\n'.format(os.path.join(path, 'out.txt'))
    str_error = 'orig_stderr = sys.stderr\ng = open(\'{}\', \'w\')\nsys.stderr = g\n'.format(os.path.join(path, 'err.txt'))
    return str_out + str_error


def end2txt():
    re_str = 'sys.stdout = orig_stdout\nsys.stderr = orig_stderr\n' \
             'f.close()\ng.close()\n'
    return re_str


'''
Files and folders, etc.
'''
def openfile():
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    filename = askopenfilename(initialdir = os.getcwd()) # show an "Open" dialog box and return the path to the selected file
    return filename

# print pickle content
def print_pickle():
    pickle_path = openfile()
    with open(pickle_path, 'rb') as handle:
        para_dict = pickle.load(handle)
        print(para_dict)
    handle.close()

# find folders with a target string
def folders_with(target_str):
    folder_list = next(os.walk('.'))[1]
    return_list = []
    for folder in folder_list:
        if target_str in folder:
            return_list.append(folder)
    return return_list

# make a folder with specific name
def make_folder(folder_name):
    data_path = ''
    dir_exist = False
    if isinstance(folder_name, str):
        data_path = os.path.join(os.getcwd(), folder_name)
        if os.path.isdir(data_path):
            dir_exist = True
        else:
            os.mkdir(data_path)
    else:
        print('not a string')
    return data_path, dir_exist


'''
From the original helpers.py
'''
def fire_rate(spikes, begin, end):
    data = spikes.get_data(begin, end)
    spikes.verify_collect('data[2]=\n{}\n'.format(data[2]), 'fr')
    gids = spikes.gids
    rates_averaged_all = []
    rates_std_all = []
    for h in list(range(len(data))):
        if len(data[h]) > 0:
            n_fil = data[h][:, 0]
            n_fil = n_fil.astype(int)
            count_of_n = np.bincount(n_fil) # bin count from 0
            count_of_n_fil = count_of_n[gids[h][0]:gids[h][1]+1] # get data with corresponding indices
            rate_each_n = count_of_n_fil * 1000. / (end - begin)
            rate_averaged = np.mean(rate_each_n)
            rate_std = np.std(rate_each_n)
            rates_averaged_all.append(float('%.3f' % rate_averaged))
            rates_std_all.append(float('%.3f' % rate_std))
            np.save(os.path.join(spikes.path, ('rate' + str(h) + '.npy')), rate_each_n)
            if h == 2:
                spikes.verify_collect('gids[2]=\n{}\n'.format(gids[2]), 'fr')
                spikes.verify_collect('len(n_fil)={}, values=\n'.format(len(n_fil)), 'fr')
                for n in n_fil:
                    spikes.verify_collect('{} '.format(n), 'fr')
                spikes.verify_collect('\n', 'fr')
                spikes.verify_collect('count_of_n(tail)=\n{}\n'.format(count_of_n[-len(count_of_n_fil):]), 'fr')
                spikes.verify_collect('count_of_n_fil=\n{}\n'.format(count_of_n_fil), 'fr')
                spikes.verify_collect('np.sum(count_of_n_fil)=\n{}\n'.format(np.sum(count_of_n_fil)), 'fr')
                spikes.verify_collect('len(count_of_n)=\n{}\n'.format(len(count_of_n)), 'fr')
                spikes.verify_collect('len(count_of_n_fil)=\n{}\n'.format(len(count_of_n_fil)), 'fr')
                spikes.verify_collect('rate_each_n=\n{}\n'.format(rate_each_n), 'fr')
                spikes.verify_collect('rate_averaged=\n{:.2f}\n'.format(rate_averaged), 'fr')
        else:
            rates_averaged_all.append(0.0)
            rates_std_all.append(0.0)
            np.save(os.path.join(spikes.path, ('rate' + str(h) + '.npy')), [])
    print('Mean rates: %r Hz' % rates_averaged_all)
    print('Standard deviation of rates: %r Hz' % rates_std_all)

    f_rates = open(os.path.join(spikes.path, 'fr.dat'), 'w')
    for rate_mean, rate_std in zip(rates_averaged_all, rates_std_all):
        f_rates.write(str(rate_mean) + ', ' + str(rate_std) + '\n')
    f_rates.close()
    spikes.fr_result = np.array([rates_averaged_all, rates_std_all])
    return rates_averaged_all, rates_std_all


def plot_raster(spikes, begin, end):
    data = spikes.get_data(begin, end)
    gids = spikes.gids
    highest_gid = gids[-1][-1]
    gids_numpy = np.asarray(gids)
    gids_numpy_changed = abs(gids_numpy - highest_gid) + 1

    # set y label parameters
    L23_label_pos = (gids_numpy_changed[0][0] + gids_numpy_changed[3][1])/2
    L4_label_pos = (gids_numpy_changed[4][0] + gids_numpy_changed[6][1])/2
    L5_label_pos = (gids_numpy_changed[7][0] + gids_numpy_changed[9][1])/2
    L6_label_pos = (gids_numpy_changed[10][0] + gids_numpy_changed[12][1])/2
    ylabels = ['L2/3', 'L4', 'L5', 'L6']

    fig, ax = plt.subplots(figsize=(16, 12))
    for i in list(range(len(data))):
        if len(data[i]) > 0:
            times = data[i][:, 1]
            neurons = np.abs(data[i][:, 0] - highest_gid) + 1
            plt.plot(times, neurons, '.', color=spikes.colors[i])

    # legend
    for i in range(4):
        plt.scatter([], [], s=100, label=spikes.subtypes[i], color=spikes.colors[i])
    plt.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.1))

    # set top and right frames invisible
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # reset y limits to contain legend
    bottom, top = ax.get_ylim()

    plt.xlabel('time (ms)')
    plt.xticks(np.arange(begin, end + 1.0, (end - begin)/4.0))
    plt.yticks(
        [L23_label_pos, L4_label_pos, L5_label_pos, L6_label_pos],
        ylabels, rotation=10
        )
    fig.tight_layout()
    plt.savefig(os.path.join(spikes.path, 'raster_plot.png'), dpi=300)
    plt.close()


def do_boxplot(data, cri, path, title, colors, ylbls, xlbl, xlims=None):
    layers = ['L2/3', 'L4', 'L5', 'L6']
    label_pos = list(range(1, len(data)+1))
    medianprops = dict(linestyle='-', linewidth=5, color='red')
    fig, ax = plt.subplots(figsize=(12, 12))
    bp = plt.boxplot(data, 0, 'k+', 0, medianprops=medianprops, patch_artist=True,)
    plt.setp(bp['boxes'], color='k')
    plt.setp(bp['whiskers'], color='k')
    if isinstance(xlims, tuple):
        plt.xlim(xlims)
    plt.ylim(label_pos[0]-0.5, label_pos[-1]+0.5)
    # box filling color
    for i, box in enumerate(bp['boxes']):
        box.set_facecolor(colors[i])
        # print(label_pos[i], cri[i])
        ax.scatter(list(cri[i])[1], label_pos[i], s=400, marker='*', edgecolor='k', color='yellow', zorder=10)
    for i in [3, 6, 9]:
        ax.hlines(label_pos[i]-0.5, plt.xlim()[0], plt.xlim()[1], linestyles='solid')
    for i in [0, 3, 6, 9]:
        ax.text(plt.xlim()[1], label_pos[i] + 1.0, layers[3 - int(i/3)], horizontalalignment='center')
    for i in range(12, 8, -1):
        plt.scatter([], [], s=100, marker='s', label=ylbls[i], color=colors[i])
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_visible(False)
    legend = plt.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.1))
    plt.xlabel(xlbl)
    fig.tight_layout()
    plt.savefig(os.path.join(path, 'boxplot_' + title + '.png'), dpi=300)
    plt.close()

def do_bars(data, cri, path, title, colors, ylbl, figsize=(15, 10)):
    layers = ['L2/3', 'L4', 'L5', 'L6']
    legends = ['Exc', 'PV', 'SOM', 'VIP']

    # bars
    x = np.arange(13)  # the label locations
    w = 0.3  # the width of the bars
    fig, ax = plt.subplots(figsize=figsize)
    rects1 = ax.bar(x - w/2, data[0, :], w, yerr=data[1, :], color=colors, edgecolor=colors)
    rects2 = ax.bar(x + w/2, cri[0, :], w, yerr=cri[1, :], fill=False, edgecolor=colors, hatch='///')

    # colors
    # for i, (r1, r2) in enumerate(zip(rects1, rects2)):
    #     r1.set_color(colors[i])
    #     r2.set_color(colors[i])
    #     r1.set_hatch('//')
    #     r2.set_hatch('xx')

    # legends
    for i in range(4):
        plt.scatter([], [], s=100, marker='s', label=legends[i], color=colors[i])
    legend = plt.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.15))

    # ticks
    ax.set_ylabel(ylbl)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_xlim((0-1.5*w, 12+1.5*w))

    # vlines
    for x in [3.5, 6.5, 9.5]:
        ax.axvline(x, 0, 1, color='k')

    for i, x in enumerate([1.5, 5, 8, 11]):
        ax.text(x, ax.get_ylim()[1], layers[i], horizontalalignment='center', verticalalignment='bottom')

    fig.tight_layout()
    plt.savefig(os.path.join(path, 'bars_' + title + '.png'))
    plt.close()

def fr_plot(spikes):
    rates = []
    for i in range(len(spikes.populations)):
        fpath = os.path.join(spikes.path, ('rate' + str(i) + '.npy'))
        if os.path.isfile(fpath):
            rates.append(np.load(fpath))
    # do_boxplot(rates[::-1], spikes.fr_qrt[::-1], spikes.path, 'fr', spikes.colors[::-1], spikes.subtypes[::-1], 'firing rate (spike/s)', xlims=(-1.0, 60.0))
    do_bars(spikes.fr_result, np.array(spikes.fr_musig).T, spikes.path, 'fr', spikes.colors, 'spikes/s')

'''
Other analysis
'''
def filter_by_spike_n(data, ids, n_spk=4):
    rdata = []
    for i in range(len(data)):
        d = data[i]
        if type(d) == np.ndarray and d.ndim == 2:
            cache = None
            for id in range(ids[i][0], ids[i][1]+1):
                tmp = d[d[:, 0] == id]
                if len(tmp) >= n_spk:
                    if cache is None:
                        cache = tmp
                    else:
                        cache = np.concatenate((cache, tmp))
            rdata.append(cache)
        else:
            rdata.append([])
    return rdata

def sample_by_layer(data, ids, layers, n_sample=140):
    rdata = []  # return data
    set_selected_by_lyr = []    # selected id sets by layer
    set_leftover_by_lyr = []    # unselected id sets by layer
    cnt_by_lyr = []             # count of selected sets by layer

    # sample n criteria: must >= 90% of desired number
    sample_cri = 0.9
    validity = np.ones(len(layers))     # validity with this criteria

    for i, layer in enumerate(layers):
        len_lyr = ids[layer[-1]][1] - ids[layer[0]][0] + 1
        cnt_lyr = 0
        selected = np.array([])
        leftover = np.array([])
        flg_exc = False
        for j, g in enumerate(layer):
            len_grp = ids[g][1] - ids[g][0] + 1
            sample_ratio = float(len_grp)/len_lyr
            d = data[g]
            if type(d) == np.ndarray and d.ndim == 2:
                set_0 = set(d[:, 0])    # original
                set_1 = sample(set_0, min(len(set_0), round(n_sample*sample_ratio))) # the set of this group that is selected
                rdata.append(d[np.in1d(d[:, 0], set_1)])    # the part of data that is in set_1
                selected = np.concatenate((selected, set_1))    # for this layer concatenate the sets that is selected
                leftover = np.concatenate((leftover, list(set_0.difference(set(set_1)))))   # for this layer concatenate the sets that is left over
                cnt_lyr += len(set_1)
                print('sample_by_layer(): group {} collected/desired n = {}/{}'.format(g, len(set_1), round(n_sample*sample_ratio)))
                if j == 0 and len(set_1) < n_sample*sample_ratio*sample_cri:
                    flg_exc = True
            else:
                rdata.append([])
        # layer sample n must > 90% desired
        if len(selected) < n_sample * sample_cri or flg_exc is True:
            validity[i] = 0
            print('layer of {} n < 90% desired; abandon this layer'.format(i))
        set_selected_by_lyr.append(selected)
        set_leftover_by_lyr.append(leftover)
        cnt_by_lyr.append(cnt_lyr)

    # makeup so that collected = desired (n_layer)
    for i, layer in enumerate(layers):
        if validity[i] == 1:    # do it only if the layer is valid
            n_diff = n_sample - cnt_by_lyr[i]    # difference of desired vs. collected
            if n_diff > 0:
                print('sample_by_layer(): layer of {} leftover n = {}'.format(i, len(set_leftover_by_lyr[i])))
                n_leftover = min(len(set_leftover_by_lyr[i]), n_diff)
                set_makeup = sample(list(set_leftover_by_lyr[i]), n_leftover)
                print('sample_by_layer(): layer of {} added {} samples'.format(i, n_leftover))
                for g in layer:
                    d = data[g]
                    if type(d) == np.ndarray and d.ndim == 2 and type(rdata[g]) == np.ndarray and rdata[g].ndim == 2:
                        rdata[g] = np.concatenate((rdata[g], d[np.in1d(d[:, 0], set_makeup)]))
            elif n_diff < 0:
                set_preserve = sample(list(set_selected_by_lyr[i]), n_sample)
                for g in layer:
                    if type(rdata[g]) == np.ndarray and rdata[g].ndim == 2:
                        rdata[g] = rdata[g][np.in1d(rdata[g][:, 0], set_preserve)]
    return rdata, validity

# Asynchronous irregular state calculation
def ai_score(spikes, begin, end, bw=10, seg_len=5000.0, layers=None, n_sample=140, n_spk=4):
    if layers is None:
        layers = [[0, 1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]
    data_all = spikes.get_data(begin, end)
    gids = spikes.gids
    seg_list = np.arange(begin, end, seg_len)

    # multiprocessing
    return_dict = Manager().dict()
    procs = []

    # ai calculation
    ai = open(os.path.join(spikes.path, 'ai.dat'), 'w')
    ai_n = open(os.path.join(spikes.path, 'ai_n.dat'), 'w')
    cvs_by_seg_by_layer = []
    validity_by_seg = []
    grp = 2
    lyr = 0
    # by segment
    for i, seg_head in enumerate(seg_list):
        seg_end = seg_head + seg_len
        data_seg = []

        # filter and sample data:
        # get segment data
        for j in range(len(data_all)):
            data_group = data_all[j]
            if type(data_group) == np.ndarray and data_group.ndim == 2:
                data_seg.append(data_group[(data_group[:, 1] >= seg_head) & (data_group[:, 1] < seg_end)])
            else:
                data_seg.append([])
        #
        if i == 0:
            spikes.verify_collect('data[grp]=\n{}\n'.format(data_all[grp]), 'gs')
            spikes.verify_collect('data_seg[grp]=\n{}\n'.format(data_seg[grp]), 'gs')

        # filter
        data_seg = filter_by_spike_n(data_seg, gids, n_spk=n_spk)
        # sample
        data_seg, valids = sample_by_layer(data_seg, gids, layers, n_sample=n_sample)
        # cache for layer data validity
        validity_by_seg.append(valids)

        #
        if i == 0:
            spikes.verify_collect('data_seg[grp](sampled)=\n{}\n'.format(data_seg[grp]), 'gs')

        # calculation:
        cvs_by_layer = []
        # by layer
        for j, layer in enumerate(layers):
            layer_ids = np.array([])
            layer_ts = np.array([])
            hists = []
            cvs = []
            cnt_corr = cnt_cv = 0
            # do it only if layer data is valid
            if valids[j] == 1:
                # get data
                for k in layer:
                    data = data_seg[k]
                    if type(data) == np.ndarray and data.ndim == 2:
                        layer_ids = np.concatenate((layer_ids, data[:, 0]))
                        layer_ts = np.concatenate((layer_ts, data[:, 1]))

                #
                if i==0 and j == lyr:
                    spikes.verify_collect('layer_ids=\n', 'gs')
                    for id in layer_ids:
                        spikes.verify_collect('{} '.format(id), 'gs')
                        spikes.verify_collect('\n', 'gs')
                    spikes.verify_collect('layer_ts=\n', 'gs')
                    for t in layer_ts:
                        spikes.verify_collect('{:.1f} '.format(t), 'gs')
                        spikes.verify_collect('\n', 'gs')

                # obtain histogram and calculate pair-corr and cv-isi
                for gid in list(set(layer_ids)):   # each id in the selected
                    ts = layer_ts[layer_ids == gid]
                    if len(ts) > 0:
                        hist_bins = np.arange(seg_head, seg_end + bw, bw)
                        hist = np.histogram(ts, bins=hist_bins)[0]
                        if not np.all(hist == hist[0]):
                            hists.append(hist)
                            cnt_corr += 1
                        if len(ts) > 3:
                            isi = np.diff(ts)
                            cvs.append(np.std(isi) / np.mean(isi))
                            cnt_cv += 1
                            #
                            if i==0 and j == lyr and gid == layer_ids[0]:
                                spikes.verify_collect('hist_bins={}\n'.format(hist_bins), 'gs')
                                spikes.verify_collect('hist={}\n'.format(hist), 'gs')
                                spikes.verify_collect('isi={}\n'.format(isi), 'gs')
                # conclude this layer
                print('seg {}, layer of {}, n of (corr, cv) = ({}, {})'.format(seg_head, j, cnt_corr, cnt_cv))
                ai_n.write('{}, {}\n'.format(cnt_corr, cnt_cv))
                # set multiprocessing
                proc = Process(target=get_corr, args=(hists, None, int(i * (len(layers)) + j), return_dict))
                procs.append(proc)
                proc.start()
            cvs_by_layer.append(cvs)
        ai_n.write('\n')
        cvs_by_seg_by_layer.append(cvs_by_layer)

    # join multiprocessing for correlation
    for proc in procs:
        proc.join()

    # collect data and save
    for j in range(len(layers)):
        corr_means_by_seg = []
        cv_means_by_seg = []
        for i in range(len(seg_list)):
            n_corrs = np.nan
            n_cvs = np.nan
            if validity_by_seg[i][j] == 1:
                corrs = return_dict[str(int(i * (len(layers)) + j))]
                corr_means_by_seg.append(np.mean(corrs))
                cvs = cvs_by_seg_by_layer[i][j]
                cv_means_by_seg.append(np.mean(cvs))
                n_corrs = len(corrs)
                n_cvs = len(cvs)
                for a, corr in enumerate(corrs):
                    if np.isnan(corr):
                        print('corr no. {} is nan'.format(a))
        ai.write('{}, {}\n'.format(str(np.mean(corr_means_by_seg)), str(np.mean(cv_means_by_seg))))

    ai.close()
    ai_n.close()

def selectivity(spikes, stims, duration, bin_w=10.0, n_bin=10, raw=False):
    n_stim = len(stims)
    if len(stims) > 1:
        interval = stims[1] - stims[0]
    data_all = spikes.get_data(stims[0], stims[-1] + interval/2)
    SIs_by_popbin = []
    xs = np.arange(0.0, bin_w*n_bin, bin_w)
    fig, axs = plt.subplots(4, 1, figsize=(8, 16), sharex=True, sharey=True, constrained_layout=True)
    ylims = np.array([-3.0, 3.0])
    plt.ylim(ylims[0], ylims[1])
    plt.xlabel('time (ms)')
    plt.ylabel('selectivity')
    axs = axs.ravel()
    safe_flg = True
    for i in range(len(data_all)):
        data_pop = data_all[i]
        if type(data_pop) == np.ndarray and data_pop.ndim == 2:
            # calculate SI by neuron
            pop_len = spikes.gids[i][1] - spikes.gids[i][0] + 1
            SI_by_neuronbin = np.full((pop_len, n_bin), np.nan)
            SI_by_neuronbin = np.full((pop_len, n_bin), np.nan)
            for j, gid in enumerate(range(spikes.gids[i][0], spikes.gids[i][1]+1)):
                # get ts by id
                ts = spikes.get_ts_by_id(data=data_pop, id=gid)
                for b in range(n_bin):
                    # get spike counts and stds by stims
                    cnts_by_stim = np.array([len(np.where((ts > stim + b*bin_w) & (ts <= stim + (b + 1)*bin_w))[0]) for stim in stims])
                    pooled_std = np.sqrt((np.var(cnts_by_stim[0::2], ddof=1) + np.var(cnts_by_stim[1::2], ddof=1))/2)
                    if pooled_std != 0: # set some more criteria, e.g. spike number must > 0 ...?
                        if raw is True:
                            SI_by_neuronbin[j, b] = (np.mean(cnts_by_stim[0::2]) - np.mean(cnts_by_stim[1::2]))/pooled_std
                        else:
                            SI_by_neuronbin[j, b] = 2*(np.abs(np.mean(cnts_by_stim[0::2]) - np.mean(cnts_by_stim[1::2])))/pooled_std    # multiplied by 2 to compare!

            # Calculate population selectivity and plot
            len_q = int(len(SI_by_neuronbin)/4)
            idx_lyr = spikes.layers[i]
            SIs_s1 = SI_by_neuronbin[len_q:3*len_q]
            SIs_s2 = np.concatenate((SI_by_neuronbin[:len_q], SI_by_neuronbin[3*len_q:]), axis=0)
            if raw is True:
                SIs_by_popbin.append(np.nanmean(SIs_s1, axis=0) - np.nanmean(SIs_s2, axis=0))
            else:
                SIs_by_popbin.append(np.nanmean(SI_by_neuronbin, axis=0))
            # main lines
            axs[idx_lyr].plot(xs, SIs_by_popbin[-1], color=spikes.colors[i], linewidth=2, label=spikes.subtypes[i])
            # boxplot
            if raw is True:
                for k, SIs in enumerate(SIs_s1.T):
                    axs[idx_lyr].boxplot(SIs[~np.isnan(SIs)], widths=[bin_w/10.0],
                    positions=[xs[k] + 2*spikes.positions[i]*bin_w/10.0],
                    boxprops=dict(color=spikes.colors[i], linewidth=2),
                    whiskerprops=dict(color=spikes.colors[i], linewidth=2),
                    flierprops=dict(markeredgecolor=spikes.colors[i], marker='.', markersize=2),
                    medianprops=dict(color=spikes.colors[i], linewidth=2),
                    capprops=dict(color=spikes.colors[i], linewidth=2))
                for k, SIs in enumerate(SIs_s2.T):
                    axs[idx_lyr].boxplot(SIs[~np.isnan(SIs)], widths=[bin_w/10.0],
                    positions=[xs[k] + (2*spikes.positions[i] + 1)*bin_w/10.0],
                    boxprops=dict(color=spikes.colors[i], linewidth=1),
                    whiskerprops=dict(color=spikes.colors[i], linewidth=1),
                    flierprops=dict(markeredgecolor=spikes.colors[i], marker='.', markersize=1),
                    medianprops=dict(color=spikes.colors[i], linewidth=1),
                    capprops=dict(color=spikes.colors[i], linewidth=1))
            else:
                for k, SIs in enumerate(SI_by_neuronbin.T):
                    axs[idx_lyr].boxplot(SIs[~np.isnan(SIs)], widths=[bin_w/10.0],
                    positions=[xs[k] + 2*spikes.positions[i]*bin_w/10.0],
                    boxprops=dict(color=spikes.colors[i], linewidth=1.5),
                    whiskerprops=dict(color=spikes.colors[i], linewidth=1.5),
                    flierprops=dict(markeredgecolor=spikes.colors[i], marker='.', markersize=1.5),
                    medianprops=dict(color=spikes.colors[i], linewidth=1.5),
                    capprops=dict(color=spikes.colors[i], linewidth=1.5))
            axs[idx_lyr].spines['right'].set_visible(False)
            axs[idx_lyr].spines['top'].set_visible(False)
            if i == 0:
                axs[idx_lyr].plot([0.0, duration], [ylims[1], ylims[1]], color='k', linewidth=10)
            if i < 4:
                axs[idx_lyr].legend(loc='upper right', fontsize=16, ncol=2)
            axs[idx_lyr].text(-25.0, np.mean(ylims), spikes.layers[idx_lyr])
    for ax in axs:
        ax.plot([0.0, bin_w*n_bin], [0.0, 0.0], color='k', linestyle='--')
    plt.xticks(xs[::2], labels=xs[::2].astype(str))
    plt.savefig(os.path.join(spikes.path, 'selectivity-raw={}.png'.format(raw)))
    plt.close()


# responses to transient/thalamic input
def response(spikes, begin, stims, window, bw=1.0, exportplot=False, interpol=False):
    n_stim = len(stims)
    if len(stims) > 1:
        interval = stims[1] - stims[0]
    data_all = spikes.get_data(begin, begin+n_stim*interval)
    # response data to file
    f = open(os.path.join(spikes.path, 'sf.dat'), 'w')

    # plot settings
    fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
    colors = ['hotpink', 'dodgerblue', 'black', 'black']
    linestyles = ['solid', 'solid', 'solid', 'dashed']
    xy_by_layer, empty_flg = [], False

    # calculate response spread and amplitude
    exc_lyrs = [0, 4, 7, 10]
    for a, i in enumerate(exc_lyrs):
        data = data_all[i]
        # skip if no data
        if type(data) != np.ndarray or data.ndim != 2:
            f.write('{}, {}\n'.format(np.nan, np.nan))
            empty_flg = True
            continue
        data = data[np.argsort(data[:, 1])] # sort by time
        rsp_amp, rsp_spread = [], [] # response amplitude and spread
        ts, ids = data[:, 1], data[:, 0]
        ltcs_by_stim = []

        # cache for neuron: [[id, sum, n], ...]
        id_set = np.arange(spikes.gids[i][0], spikes.gids[i][1]+1)
        neuron_ltc_cache = np.array([id_set, np.zeros(id_set.shape), np.zeros(id_set.shape)]).T

        # loop by stimulation
        for j in range(n_stim):
            # define start and end of window
            win_begin = stims[j]
            win_end = stims[j] + window
            # spread and amplitude
            ts_stim = ts[(ts > win_begin) & (ts <= win_end)]
            rsp_amp.append(len(ts_stim))
            if len(ts_stim) >= 3:
                std = np.std(ts_stim)
                rsp_spread.append(std)
            # latency by neurons
            ids_stim = ids[(ts > win_begin) & (ts <= win_end)]
            neuron_ltcs = []
            for k, id in enumerate(id_set):
                if id not in ids_stim:
                    continue
                idx = ids_stim.tolist().index(id)
                neuron_ltc_cache[k, 1] += ts_stim[idx] - win_begin   # latency
                neuron_ltc_cache[k, 2] += 1  # sample n

        # calculate:
        # include neurons with n of responses > n_stim/2
        neuron_ltc_cache = neuron_ltc_cache[np.where(neuron_ltc_cache[:, 2]>=n_stim/2)]
        # neuron mean latencies
        mean_ltcs = np.sort(np.divide(neuron_ltc_cache[:, 1], neuron_ltc_cache[:, 2]))
        # layer average latency
        lyr_avg_ltc = np.mean(mean_ltcs)
        # for raster plot vline marking (not using)
        spikes.react_lines.append(lyr_avg_ltc + win_begin)

        # plot:
        hist, bins = np.histogram(mean_ltcs, bins=np.arange(0.0, window, bw))
        xs, ys = bins[:-1], hist/np.sum(hist)
        xy_by_layer.append([xs, ys])
        f_cubic = interpolate.interp1d(xs, ys, kind='quadratic')
        # interpolate
        if interpol:
            xs = np.linspace(min(xs), max(xs), len(xs)*5)
            ys = f_cubic(xs)
        ax.plot(xs, ys,
            linestyle=linestyles[a],
            linewidth=3,
            label=spikes.populations[i],
            color=colors[a])
        ax.set_xlabel('mean spike latency (ms)')
        ax.set_ylabel('fraction')
        ax.set_ylim(top=0.41)

        # save response spread, amplitude, layer average latency
        f.write('{:.2f}, {:.2f}, {:.2f}\n'.format(np.mean(rsp_spread), np.mean(rsp_amp), lyr_avg_ltc))
    f.close()
    plt.legend()
    # export latency plot
    plt.savefig(os.path.join(spikes.path, 'hist-ltc.png'))
    if exportplot:
        os.system('cp -p {} {}'.format(
        os.path.join(spikes.path, 'hist-ltc.png'),
        spikes.path.split('/')[-1] + '_hist-ltc.png'))
    # save neuron latency data
    if empty_flg:
        xy_by_layer = []
    np.save(os.path.join(spikes.path, 'lts_distr.npy'), xy_by_layer)


def paradox_plot(spikes, xs, frs_all, normalize=False, ylims=(0.0, 5.0), shrink=8, frs_ei=None):
    fig, axs = plt.subplots(4, 1, figsize=(6, 15), sharex=True, sharey=True)
    purple = (170/255, 51/255, 119/255)
    plt.xlabel('I (pA)')
    affix = 'normalized' if normalize else ''
    plt.ylabel('{} r (spike/s)'.format(affix))
    for i, gids in enumerate(spikes.gids):
        # axs[spikes.layers[i]].hlines(1., params_dict['offsets'][0], params_dict['offsets'][-1], linestyles=':', color='k')
        if normalize is False and spikes.subtypes[i] in ['PV', 'VIP']:
            frs = np.array(frs_all[i])/shrink
        else:
            frs = np.array(frs_all[i])
        axs[spikes.layers[i]].plot(xs, frs, color=spikes.colors[i], linewidth=2, marker='.', markersize=10)
    if frs_ei is not None:
        for i, frs in enumerate(frs_ei):
            axs[i].plot(xs, frs[0], color=spikes.colors[0], linewidth=4, marker='.', markersize=10)
            if normalize:
                axs[i].plot(xs, frs[1], color=purple, linewidth=4, marker='.', markersize=10)
            else:
                axs[i].plot(xs, frs[1]/shrink, color=purple, linewidth=4, marker='.', markersize=10)
    # legend
    for i in range(4):
        axs[0].plot([], [], label=spikes.subtypes[i], color=spikes.colors[i])
    axs[0].legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.5))
    if normalize is False and isinstance(ylims, tuple):
        plt.ylim(ylims)
    else:
        plt.ylim((0., 1.5))
    plt.savefig(os.path.join(spikes.path, 'paradox_{}.png'.format(affix)), bbox_inches='tight')
    plt.close()

def paradox_calc(spikes, prdx_dict):
    if prdx_dict['n'] <= 0:
        return
    # f = open(os.path.join(spikes.path, 'paradox.dat'), 'w')
    n_trial, duration = prdx_dict['n'], prdx_dict['duration']
    frs_all, frs_all_norm = [], []
    spksum_by_pop, ns_by_pop = [], []
    for i, gids in enumerate(spikes.gids):
        spksum_by_stim = []
        frs, n_pop = [], gids[1] - gids[0] + 1
        # by stim. level
        for j, (offset, starts) in enumerate(prdx_dict['starts'].items()):
            spk_sum, fr = 0, np.nan
            # by trial
            for k, start in enumerate(starts):
                data = spikes.get_data(start+duration/2, start+duration)[i]
                if type(data) != np.ndarray or data.ndim != 2:
                    spk_sum = np.nan
                    break
                spk_sum += len(data)
            spksum_by_stim.append(spk_sum)
            fr = 1000*(spk_sum/n_pop)/(n_trial*(duration/2))
            frs.append(fr)
            # f.write('{:.3f}, '.format(fr))
        spksum_by_pop.append(spksum_by_stim)
        frs_all.append(frs)
        norm = np.array(frs)/frs[0] if frs[0] != 0 else np.full(len(prdx_dict['offsets']), np.nan)
        frs_all_norm.append(norm)
        # f.write('\n')
        ns_by_pop.append(n_pop)
    spksum_by_pop = np.array(spksum_by_pop)
    # print(frs_all)
    # frs_all = (1000/(duration/2))*spksum_by_pop/n_trial
    # print(frs_all)
    # frs_all_norm = np.divide(frs_all, np.tile(frs_all[:, 0], (len(frs_all[0]), 1)).T)
    ns_by_pop = np.array(ns_by_pop)
    spk_sum_by_ei, n_sum_by_ei = [], []
    for i in range(4):
        idxs = np.where(np.array(spikes.layers) == i)[0]
        spk_sum_by_ei.append([spksum_by_pop[idxs[0], :]/ns_by_pop[idxs[0]],
        np.sum(spksum_by_pop[idxs[1:], :]/np.sum(ns_by_pop[idxs[1:]]), axis=0)])
    frs_by_ei = (1000/(duration/2))*np.array(spk_sum_by_ei)/n_trial
    frs_by_ei_norm = copy.deepcopy(frs_by_ei)
    for i, frs_lyr in enumerate(frs_by_ei_norm):
        for j, frs in enumerate(frs_lyr):
            frs_by_ei_norm[i, j] /= frs_by_ei_norm[i, j, 0]
    # f.close()
    paradox_plot(spikes, prdx_dict['offsets'], frs_all, normalize=False, frs_ei=frs_by_ei)
    paradox_plot(spikes, prdx_dict['offsets'], frs_all_norm, normalize=True, frs_ei=frs_by_ei_norm)
