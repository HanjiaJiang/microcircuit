import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib
matplotlib.rcParams['font.size'] = 20.0
import numpy as np
from random import sample
from multiprocessing import Process
from multiprocessing import Manager
import pickle
from tkinter import Tk
from tkinter.filedialog import askopenfilename

'''
Notes:
way to get file path by GUI:
import easygui
path = easygui.fileopenbox()
'''

populations = ['L2/3 Exc', 'L2/3 PV', 'L2/3 SOM', 'L2/3 VIP',
               'L4 Exc', 'L4 PV', 'L4 SOM',
               'L5 Exc', 'L5 PV', 'L5 SOM',
               'L6 Exc', 'L6 PV', 'L6 SOM']

subtype_label = ['Exc', 'PV', 'SOM', 'VIP']

use_box_xlim = True
box_xlim = 51.0
cv_isi_min_n = 4

'''
Plots
'''
def interaction_barplot(arr, y_bottom, y_top, labels=None, ylabel=None):
    arr_len = len(arr)
    arr_shape = np.array(arr).shape
    if len(arr_shape) is not 2:
        print('array not a 2-D array')
        return
    elif arr_shape[0] != arr_shape[1]:
        print('array not a square array')
        return

    if isinstance(labels, list) is not True:
        print('labels not a list')
        labels = list(range(arr_len))
    elif len(labels) > arr_len:
        print('labels too long, use the first {}'.format(arr_len))
        labels = labels[:len(arr)]
    elif len(labels) < arr_len:
        print('labels too short')
        labels = list(range(arr_len))

    labels = np.array(labels).astype(str)
    x = np.arange(arr_len)
    barwidth = 0.75/arr_len
    fig, ax = plt.subplots(figsize=(12, 12))
    for i in range(len(arr)):
        ax.bar(x + barwidth * i, arr[i, :], barwidth, label=labels[i])
    ax.legend(bbox_to_anchor=(0., 1.2, 1.0, 0.1),
              ncol=int(len(labels)/2), mode="expand", borderaxespad=0.)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    # ax.legend()
    plt.ylim((y_bottom, y_top))
    fig.tight_layout()
    plt.savefig('interaction_barplot.png')
    plt.show()

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
def set2txt():
    re_str = 'import sys\n' \
             'orig_stdout = sys.stdout\n' \
             'f = open(\'out.txt\', \'w\')\n' \
             'sys.stdout = f'
    return re_str


def end2txt():
    re_str = 'sys.stdout = orig_stdout\n' \
             'f.close()'
    return re_str


'''
Files and folders, etc.
'''
def openfile():
    Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
    filename = askopenfilename(initialdir = "~/Documents/") # show an "Open" dialog box and return the path to the selected file
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
def read_name(path, name):
    # Import filenames
    files = []
    for file in os.listdir(path):
        if file.endswith('.gdf') and file.startswith(name):
            temp = file.split('-')[0] + '-' + file.split('-')[1]
            if temp not in files:
                files.append(temp)

    # Import GIDs
    gidfile = open(os.path.join(path, 'population_GIDs.dat'), 'r')
    gids = []
    for l in gidfile:
        a = l.split()
        gids.append([int(a[0]), int(a[1])])
    files = sorted(files)
    return files, gids


def load_spike_times(path, name, begin, end):
    detector_names, gids = read_name(path, name)    # names: populations
    data = {}
    if len(detector_names) > 0 and len(gids) > 0:
        for i in list(range(len(detector_names))):
            all_filenames = os.listdir(path)
            thread_filenames = [
                all_filenames[x] for x in list(range(len(all_filenames)))
                if all_filenames[x].endswith('gdf') and
                   all_filenames[x].startswith('spike') and
                   (all_filenames[x].split('-')[0]
                    + '-' + all_filenames[x].split('-')[1]) in
                   detector_names[i]
                ]
            # ids_temp = np.array([])
            # times_temp = np.array([])
            data_temp = []
            for f in thread_filenames:
                load_tmp = np.array([])
                try:
                    load_tmp = np.loadtxt(os.path.join(path, f))
                except ValueError:
                    print(os.path.join(path, f))
                else:
                    pass
                if len(load_tmp.shape) == 2:
                    data_temp.append(load_tmp)
                    # ids_temp = np.append(ids_temp, load_tmp[:, 0])
                    # times_temp = np.append(times_temp, load_tmp[:, 1])
            # data_concatenated = np.array([ids_temp, times_temp]).T  # transpose
            if len(data_temp) > 0:
                data_concatenated = np.concatenate(data_temp)
                data_raw = \
                    data_concatenated[np.argsort(data_concatenated[:, 1])]
                idx = ((data_raw[:, 1] > begin) * (data_raw[:, 1] < end))
                data[i] = data_raw[idx]
            else:
                data[i] = []
    return data, gids   # an extra return: gids


def fire_rate(path, name, begin, end):
    files, gids = read_name(path, name)
    data_all, gids = load_spike_times(path, name, begin, end)
    rates_averaged_all = []
    rates_std_all = []
    for h in list(range(len(files))):
        if len(data_all[h]) > 0:
            n_fil = data_all[h][:, 0]
            n_fil = n_fil.astype(int)
            count_of_n = np.bincount(n_fil)
            count_of_n_fil = count_of_n[gids[h][0]-1:gids[h][1]]
            rate_each_n = count_of_n_fil * 1000. / (end - begin)
            rate_averaged = np.mean(rate_each_n)
            rate_std = np.std(rate_each_n)
            rates_averaged_all.append(float('%.3f' % rate_averaged))
            rates_std_all.append(float('%.3f' % rate_std))
            np.save(os.path.join(path, ('rate' + str(h) + '.npy')),
                    rate_each_n)
        else:
            #190606
            rates_averaged_all.append(0.0)
            rates_std_all.append(0.0)
            np.save(os.path.join(path, ('rate' + str(h) + '.npy')), [])
    print('Mean rates: %r Hz' % rates_averaged_all)
    print('Standard deviation of rates: %r Hz' % rates_std_all)

    f_rates = open(os.path.join(path, 'fr.dat'), 'w')
    for rate_mean, rate_std in zip(rates_averaged_all, rates_std_all):
        f_rates.write(str(rate_mean) + ', ' + str(rate_std) + '\n')
    f_rates.close()
    return rates_averaged_all, rates_std_all


def plot_raster(path, name, begin, end):
    files, gids = read_name(path, name)
    data_all, gids = load_spike_times(path, name, begin, end)
    highest_gid = gids[-1][-1]
    gids_numpy = np.asarray(gids)
    gids_numpy_changed = abs(gids_numpy - highest_gid) + 1

    # set y label parameters
    L23_label_pos = (gids_numpy_changed[0][0] + gids_numpy_changed[3][1])/2
    L4_label_pos = (gids_numpy_changed[4][0] + gids_numpy_changed[6][1])/2
    L5_label_pos = (gids_numpy_changed[7][0] + gids_numpy_changed[9][1])/2
    L6_label_pos = (gids_numpy_changed[10][0] + gids_numpy_changed[12][1])/2
    ylabels = ['L2/3', 'L4', 'L5', 'L6']

    color_list = [
        'blue', 'red', 'orange', 'green', 'blue', 'red', 'orange',
        'blue', 'red', 'orange', 'blue', 'red', 'orange'
    ]

    fig, ax = plt.subplots(figsize=(16, 12))
    for i in list(range(len(files))):
        if len(data_all[i]) > 0:
            times = data_all[i][:, 1]
            neurons = np.abs(data_all[i][:, 0] - highest_gid) + 1
            # 190819
            if i < 4:
                plt.plot(times, neurons, '.', color=color_list[i],
                         label=subtype_label[i])
            else:
                plt.plot(times, neurons, '.', color=color_list[i])

    # legend handling
    legend = plt.legend(loc='upper center', ncol=4)
    for legend_handle in legend.legendHandles:
        legend_handle._legmarker.set_markersize(30)

    # set top and right frames invisible
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # reset y limits to contain legend
    bottom, top = ax.get_ylim()
    ax.set_ylim((bottom, top + 1500.0))

    # plt.title('(A) raster plot', fontsize=40.0)
    plt.xlabel('time (ms)')
    plt.xticks(np.arange(begin, end + 1.0, (end - begin)/4.0))
    plt.yticks(
        [L23_label_pos, L4_label_pos, L5_label_pos, L6_label_pos],
        ylabels, rotation=10
        )
    fig.tight_layout()
    plt.savefig(os.path.join(path, 'raster_plot.png'), dpi=300)
    plt.close()


def do_boxplot(data, path, title, xlbls, ylbls, clr_list, xlims):
    label_pos = list(range(len(data), 0, -1))

    clr_list = clr_list[::-1]
    medianprops = dict(linestyle='-', linewidth=2.5, color='black')
    fig, ax = plt.subplots(figsize=(12, 12))
    bp = plt.boxplot(data, 0, 'k+', 0, medianprops=medianprops)
    plt.xlim(xlims)
    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    for h in range(len(data)):
        boxX = []
        boxY = []
        box = bp['boxes'][h]
        for j in list(range(5)):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = list(zip(boxX, boxY))
        boxPolygon = Polygon(boxCoords, facecolor=clr_list[h])
        ax.add_patch(boxPolygon)
    # set top and right frames invisible
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.xlabel(xlbls)
    plt.yticks(label_pos, ylbls)
    fig.tight_layout()
    plt.savefig(os.path.join(path, 'box_plot_' + title + '.png'), dpi=300)
    plt.close()

def fr_boxplot(net_dict, path):
    pops = net_dict['N_full']
    # list of population length e.g. [0, 1, ..., 12]
    reversed_order_list = list(range(len(pops) - 1, -1, -1))
    list_rates_rev = []
    for h in reversed_order_list:
        rate_filepath = os.path.join(path, ('rate' + str(h) + '.npy'))
        if os.path.isfile(rate_filepath):
            list_rates_rev.append(
                np.load(rate_filepath)
                )
    color_list = [
        'blue', 'red', 'orange', 'green', 'blue', 'red', 'orange',
        'blue', 'red', 'orange', 'blue', 'red', 'orange'
    ]
    do_boxplot(list_rates_rev, path, 'fr', 'firing rate (spikes/s)', populations, color_list, (-1.0, 50.0))


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
    for i, layer in enumerate(layers):
        len_lyr = ids[layer[-1]][1] - ids[layer[0]][0] + 1
        cnt_lyr = 0
        selected = np.array([])
        leftover = np.array([])
        for g in layer:
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
            else:
                rdata.append([])
        set_selected_by_lyr.append(selected)
        set_leftover_by_lyr.append(leftover)
        cnt_by_lyr.append(cnt_lyr)

    # makeup so that collected = desired (n_layer)
    for i, layer in enumerate(layers):
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
    return rdata

# Asynchronous irregular state calculation
def ai_score(path, name, begin, end, bw=10, seg_len=5000.0, layers=None, n_sample=140, n_spk=4):
    if layers is None:
        layers = [[0, 1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]
    data_all, gids = load_spike_times(path, name, begin, end)
    seg_list = np.arange(begin, end, seg_len)

    # multiprocessing
    return_dict = Manager().dict()
    procs = []

    # calculation and save
    ai = open(os.path.join(path, 'ai.dat'), 'w')
    cvs_by_seg_by_layer = []
    for i, seg_head in enumerate(seg_list):
        seg_end = seg_head + seg_len
        data_seg = []
        for j in range(len(data_all)):
            data_group = data_all[j]
            if type(data_group) == np.ndarray and data_group.ndim == 2:
                data_seg.append(data_group[(data_group[:, 1] >= seg_head) & (data_group[:, 1] < seg_end)])
            else:
                data_seg.append([])
        data_seg = filter_by_spike_n(data_seg, gids, n_spk=n_spk)
        data_seg = sample_by_layer(data_seg, gids, layers, n_sample=n_sample)

        cvs_by_layer = []
        for j, layer in enumerate(layers):
            layer_ids = np.array([])
            layer_ts = np.array([])
            hists = []
            cvs = []
            cnt_corr = cnt_cv = 0
            for k in layer:
                data = data_seg[k]
                if type(data) == np.ndarray and data.ndim == 2:
                    layer_ids = np.concatenate((layer_ids, data[:, 0]))
                    layer_ts = np.concatenate((layer_ts, data[:, 1]))
            for k in layer:
                for gid in range(gids[k][0], gids[k][1] + 1):   # each neuron id of this group
                    ts = layer_ts[layer_ids == gid]
                    if len(ts) > 0:
                        hists.append(
                            np.histogram(ts, bins=np.arange(seg_head, seg_end + bw, bw))[0])
                        cnt_corr += 1
                        if len(ts) > 3:
                            isi = np.diff(ts)
                            cvs.append(np.std(isi) / np.mean(isi))
                            cnt_cv += 1
            print('seg {}, layer of {}, n of (corr, cv) = ({}, {})'.format(seg_head, j, cnt_corr, cnt_cv))
            proc = Process(target=get_corr,
                           args=(hists, None, int(i * (len(layers)) + j), return_dict))
            procs.append(proc)
            proc.start()
            cvs_by_layer.append(cvs)
        cvs_by_seg_by_layer.append(cvs_by_layer)

    for proc in procs:
        proc.join()
    for j in range(len(layers)):
        corr_means_by_seg = []
        cv_means_by_seg = []
        for i in range(len(seg_list)):
            corrs = return_dict[str(int(i * (len(layers)) + j))]
            corr_means_by_seg.append(np.mean(corrs))
            cvs = cvs_by_seg_by_layer[i][j]
            cv_means_by_seg.append(np.mean(cvs))
        ai.write(str(np.mean(corr_means_by_seg)) + ', ' + str(np.mean(cv_means_by_seg)) + '\n')

    ai.close()
    # do_boxplot(corrs_by_layer, path, 'pair-corr', 'pairwise correlation',
    #            ['L2/3', 'L4', 'L5', 'L6'], ['gray', 'gray', 'gray', 'gray'], (-1.0, 1.0))
    # do_boxplot(cvs_by_layer, path, 'cv-isi', 'CV of ISI',
    #            ['L2/3', 'L4', 'L5', 'L6'], ['gray', 'gray', 'gray', 'gray'], (-0.1, 2.0))


def ai_score_200(path, name, begin, end,
             limit=200, bw=10, filter_1hz=False, seg_len = 5000.0):
    data_all, gids = load_spike_times(path, name, begin, end)
    seg_list = np.arange(begin, end, seg_len)

    # selected neurons with fr > 1 Hz
    pass_ids_all = []

    # group number in each layer
    layers = [[0, 1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]
    layer_head = [0, 0, 0, 0, 4, 4, 4, 7, 7, 7, 10, 10, 10]
    layer_tail = [3, 3, 3, 3, 6, 6, 6, 9, 9, 9, 12, 12, 12]

    # sampling for each group
    for i in range(len(data_all)):
        pass_ids_group = []
        if len(data_all[i]) > 0:
            ids = data_all[i][:, 0]
            ts = data_all[i][:, 1]
            n_group = gids[i][1] - gids[i][0] + 1
            n_layer = gids[layer_tail[i]][1] - gids[layer_head[i]][0] + 1
            n_desired = int(limit * n_group / n_layer)
            # sampling for each segment
            for seg_head in seg_list:
                seg_ids = ids[(ts >= seg_head) & (ts < seg_head + seg_len)]
                pass_ids_seg = []
                for gid in range(gids[i][0], gids[i][-1]+1):
                    if filter_1hz:
                        if len(seg_ids[seg_ids == gid]) >= seg_len / 1000.0:
                            pass_ids_seg.append(gid)
                    else:
                        if len(seg_ids[seg_ids == gid]) > 3:
                            pass_ids_seg.append(gid)
                if len(pass_ids_seg) > 0:
                    pass_ids_seg = sample(pass_ids_seg, min(len(pass_ids_seg), n_desired))
                print('group {}, seg {}, desired and obtained n = {}, {}'.format(i, seg_head, n_desired, len(pass_ids_seg)))
                pass_ids_group.append(pass_ids_seg)
                # print(pass_ids_group)
        else:
            for seg_head in seg_list:
                pass_ids_group.append([])
        pass_ids_all.append(pass_ids_group)
        # print(pass_ids_group[:10])


    # calculation and save
    ai = open(os.path.join(path, 'ai.dat'), 'w')
    corrs_by_layer = []
    cvs_by_layer = []
    for i, layer in enumerate(layers):
        layer_ids = np.array([])
        layer_ts = np.array([])
        for j in layer:
            if len(data_all[j]) > 0:
                layer_ids = np.concatenate((layer_ids, data_all[j][:, 0].astype(int)))
                layer_ts = np.concatenate((layer_ts, data_all[j][:, 1]))

        # for layer
        corr_means_lyr = []
        cv_means_lyr = []
        corrs_lyr = np.array([])
        cvs_lyr =  np.array([])

        # correlation and irregularity
        for j, seg_head in enumerate(seg_list):
            seg_end = seg_head + seg_len
            hists = []
            cvs = []
            sample_cnt = 0
            seg_ts = layer_ts[(layer_ts >= seg_head) & (layer_ts < seg_end)]
            seg_ids = layer_ids[(layer_ts >= seg_head) & (layer_ts < seg_end)]
            for k in layer:
                for gid in pass_ids_all[k][j]:
                    ts = seg_ts[seg_ids == gid]
                    if len(ts) > 3:
                        hists.append(
                            np.histogram(ts, bins=np.arange(seg_head, seg_end + bw, bw))[0])
                        isi = np.diff(ts)
                        cvs.append(np.std(isi) / np.mean(isi))
                        sample_cnt += 1
            print('layer {}, seg {}, n = {}'.format(i, seg_head, sample_cnt))
            corrs = get_corr(hists)
            corrs_lyr = np.concatenate((corrs_lyr, corrs))
            corr_means_lyr.append(np.mean(corrs))
            cvs_lyr = np.concatenate((cvs_lyr, cvs))
            cv_means_lyr.append(np.mean(cvs))
        ai.write(str(np.mean(corr_means_lyr)) + ', ' + str(np.mean(cv_means_lyr)) + '\n')
        corrs_by_layer.append(corrs_lyr)
        cvs_by_layer.append(cvs_lyr)
    ai.close()
    do_boxplot(corrs_by_layer, path, 'pair-corr', 'pairwise correlation',
               ['L2/3', 'L4', 'L5', 'L6'], ['gray', 'gray', 'gray', 'gray'], (-1.0, 1.0))
    do_boxplot(cvs_by_layer, path, 'cv-isi', 'CV of ISI',
               ['L2/3', 'L4', 'L5', 'L6'], ['gray', 'gray', 'gray', 'gray'], (-0.1, 2.0))


# response spread and amplitude
def response(path, name, begin, window, n_stim=20, interval=1000.0):
    data_all, gids = load_spike_times(path, name, begin, begin+n_stim*interval)
    data_save = np.full((13, n_stim, 2), np.nan)
    f = open(os.path.join(path, 'sf.dat'), 'w')
    for i in range(len(data_all)):
        if 'Exc' in populations[i]:
            data = data_all[i]
            t_stds = []
            n_spikes_list = []
            if len(data) > 0:
                ts = data[:, 1]
                for j in range(n_stim):
                    ts_sf = ts[(ts > begin + j*interval) & (ts <= begin + j*interval+window)]
                    n_spikes_list.append(len(ts_sf))
                    if len(ts_sf) >= 3:
                        std = np.std(ts_sf)
                        t_stds.append(std)
                    else:
                        std = np.nan
                    data_save[i, j, 0] = std
                    data_save[i, j, 1] = len(ts_sf)
            f.write('{:.2f}, {:.2f}\n'.format(np.mean(t_stds), np.mean(n_spikes_list)))
    f.close()
    np.save(os.path.join(path, 'sf.npy'), data_save)


# to be improved ..
def plot_psth(path, name, begin, end):
    files, gids = read_name(path, name)
    data_all, gids = load_spike_times(path, name, begin, end)
    fig, axs = plt.subplots(4, 1, figsize=(12, 12),
                            sharex=True, sharey=True, constrained_layout=True)
    colors = ['b', 'r', 'orange', 'g',
              'b', 'r', 'orange',
              'b', 'r', 'orange',
              'b', 'r', 'orange']
    bin_width = 1.0
    for i in list(range(len(data_all))):
        # times = np.array([])
        if len(data_all[i]) > 0:
            times = data_all[i][:, 1]
            # times = times.astype(int)

            # indexing for plots
            if i < 4:
                a = 0
                # b = i % 3.0
                # if i == 3:
                #     b = 3.0
            else:
                a = int((i - 1) / 3)
                # b = (i - 1) % 3.0

            # bar_width = window / 5.0
            # axs[a].bar(t_plot+b*bar_width, si_abs_arr[:, i], width=bar_width,
            #            label=populations[i] + ' all cells (abs)', color=colors[i], align='center')
            axs[a].hist(times, np.arange(begin, end + bin_width, bin_width), color=colors[i], label=populations[i])
            axs[a].legend(loc='upper right')
            # axs[a].set_ylim(-0.7, 0.7)
            # axs[a].set_xlim(t_plot[0] - 10.0, t_plot[-1] + 50.0)

    plt.xlabel('t (ms)')
    plt.ylabel('spikes')
    plt.savefig(os.path.join(path, 'psth.png'))
    plt.close()
