import os
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib
matplotlib.rcParams['font.size'] = 30.0
import numpy as np

'''
Notes:
way to get file path by GUI:
import easygui
path = easygui.fileopenbox() 
'''

populations = ['L2/3 E', 'L23 PV', 'L23 SOM', 'L23 VIP',
               'L4 E', 'L4 PV', 'L4 SOM',
               'L5 E', 'L5 PV', 'L5 SOM',
               'L6 E', 'L6 PV', 'L6 SOM']

subtype_label = ['E', 'PV', 'SOM', 'VIP']

use_box_xlim = True
box_xlim = 51.0

# Plots
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

# Calculation
# correlation
def get_mean_corr(list_1, list_2=None):
    coef_list = []
    for i, hist1 in enumerate(list_1):
        if list_2 is None:  # same population
            list_2_tmp = list_1[i + 1:]
        else:
            list_2_tmp = list_2
        for j, hist2 in enumerate(list_2_tmp):
            if np.sum(hist1) != 0 and np.sum(hist2) != 0:
                coef = np.corrcoef(hist1, hist2)[0, 1]
                coef_list.append(coef)
            # else:
            #     print('{} in list1 or {} in list2 no spikes'.format(i, j))
    return np.mean(coef_list)

# System
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


# File & folder, etc.
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


# Potjans_2014 helpers
def read_name(path, name):
    # Import filenames$
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
    # if check_data(path, name):
    #     return data_dict['data'], data_dict['gids']
    # else:
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
    # data_dict['path'] = path
    # data_dict['name'] = name
    # data_dict['gids'] = gids
    # data_dict['data'] = data
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
    f_rates.write('Mean rates: %r Hz\n' % rates_averaged_all)
    f_rates.write('Standard deviation of rates: %r Hz\n' % rates_std_all)
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
            # plt.plot(times, neurons, '.', color=color_list[i])
            # 190819
            if i < 4:
                plt.plot(times, neurons, '.', color=color_list[i],
                         label=subtype_label[i])
            else:
                plt.plot(times, neurons, '.', color=color_list[i])

    # legend handling
    #legend = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),
    #       ncol=4, mode="expand", borderaxespad=0.)
    legend = plt.legend(loc='upper center', ncol=4)
    for legend_handle in legend.legendHandles:
        legend_handle._legmarker.set_markersize(30)

    # set top and right frames invisible
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # reset y limits to contain legend
    bottom, top = ax.get_ylim()
    ax.set_ylim((bottom, top + 1500.0))

    plt.title('(A) raster plot', fontsize=40.0)
    plt.xlabel('time (ms)')
    plt.xticks(np.arange(begin, end + 1.0, (end - begin)/4.0))
    plt.yticks(
        [L23_label_pos, L4_label_pos, L5_label_pos, L6_label_pos],
        ylabels, rotation=10
        )
    fig.tight_layout()
    plt.savefig(os.path.join(path, 'raster_plot.png'), dpi=300)
    plt.close()




def boxplot(net_dict, path):
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
    pop_names = populations
    label_pos = list(range(len(pops), 0, -1))
    color_list = [
        'blue', 'red', 'orange', 'green', 'blue', 'red', 'orange',
        'blue', 'red', 'orange', 'blue', 'red', 'orange'
    ]
    color_list = color_list[::-1]
    medianprops = dict(linestyle='-', linewidth=2.5, color='black')
    fig, ax1 = plt.subplots(figsize=(12, 12))
    bp = plt.boxplot(list_rates_rev, 0, 'k+', 0, medianprops=medianprops)

    plt.title('(B) firing rate', fontsize=40.0)
    if use_box_xlim:
        plt.xlim((-1, box_xlim))

    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    #plt.setp(bp['fliers'], color='black', marker='+')
    for h in list(range(len(pops))):
        boxX = []
        boxY = []
        box = bp['boxes'][h]
        for j in list(range(5)):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = list(zip(boxX, boxY))
        boxPolygon = Polygon(boxCoords, facecolor=color_list[h])
        ax1.add_patch(boxPolygon)

    # set top and right frames invisible
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    plt.xlabel('firing rate (Hz)')
    plt.yticks(label_pos, pop_names)
    fig.tight_layout()
    plt.savefig(os.path.join(path, 'box_plot.png'), dpi=300)
    plt.close()


# Other analysis
# response spread and amplitude
def response(path, name, begin, window, n_stim=20, interval=1000.0):
    # end = begin + window
    data_all, gids = load_spike_times(path, name, begin, begin+n_stim*interval)
    data_save = np.full((13, n_stim, 2), np.nan)
    f = open(os.path.join(path, 'sf-chain.txt'), 'w')
    f.write('window = {:.2f} ms\n'.format(window))
    for i in range(len(data_all)):
        data = data_all[i]
        if 'E' in populations[i]:
            # print(populations[i]+'\n')
            f.write(populations[i]+'\n')
        if len(data) > 0:
            for j in range(n_stim):
                times_sf = data[(data[:, 1] > begin + j*interval) &
                                (data[:, 1] <= begin + j*interval+window), 1]
                if len(times_sf) > 0:
                    std_sf = np.std(times_sf)
                else:
                    std_sf = np.nan
                if 'E' in populations[i]:
                    tmp_str = '{:.2f},{:d}\n'.format(std_sf, len(times_sf))
                    # print(tmp_str)
                    f.write(tmp_str)
                data_save[i, j, 0] = std_sf
                data_save[i, j, 1] = len(times_sf)
    f.close()
    np.save(os.path.join(path, 'sf-chain.npy'), data_save)


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