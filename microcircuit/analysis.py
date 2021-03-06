import os
import copy
import time
import pickle
import pandas as pd
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
import matplotlib.ticker as ticker
matplotlib.rcParams['font.size'] = 20.0
matplotlib.rcParams['hatch.linewidth'] = 2.0

class Spikes:
    def __init__(self, path, net_dict):
        self.setup(path, net_dict)
        t0 = time.time()
        self.read_name()
        t1 = time.time()
        self.load()
        t2 = time.time()
        print('read_name(): {:.3f}'.format(t1-t0))
        print('load(): {:.3f}'.format(t2-t1))

    def fake_sample(self, mu, sig, n, estimate=True):
        if estimate:
            sig = sig*(np.sqrt(n-1)/np.sqrt(n))
        fsample = np.full(n, mu)
        for i in range(int(n/2)):
            fsample[i] += sig
            fsample[-1-i] -= sig
        if n%2 == 1:
            fsample[0] += sig/2
            fsample[-1] -= sig/2
        # print('{:.4f}({:.4f}) vs. {:.4f}({:.4f})'.format(mu, sig, np.mean(fsample), np.std(fsample)))
        return fsample

    def setup(self, path, net_dict):
        # data
        self.path = path
        self.device_types, self.devices = net_dict['rec_dev'], {}
        self.gids, self.data = [], {}
        # net_dict
        self.net_dict = net_dict
        # firing rate results
        self.fr_result = {}
        # firing rate criteria
        self.fr_cri = {
            'musig': [
                (2.7, 3.7), (13.8, 8.9), (2.6, 3.6), (14.6, 7.3),
                (0.5, 0.8), (10.2, 7.2), (2.6, 3.2),
                (6.8, 5.2), (7.5, 5.2), (2.8, 4.5),
                (6.1, 6.9), (16.9, 14.3), (3.9, 4.9)],
            'quarters': [
                (0.5, 0.6, 4.5), (7.5, 11.7, 23.3), (0.03, 0.4, 4.1), (8.5, 11.1, 21.0),
                (0.0, 0.1, 0.7), (4.3, 7.8, 14.7), (0.3, 0.6, 4.9),
                (2.7, 5.2, 11.2), (4.3, 7.6, 8.7), (0.2, 0.8, 3.6),
                (0.4, 2.6, 11.5), (4.6, 17.2, 22.0), (0.5, 1.7, 6.9)],
            'n': [5, 8, 9, 9, 95, 43, 27, 23, 7, 18, 30, 15, 26],
        }
        cri_sample = [self.fake_sample(mu, sig, n)
            for i, (mu, sig, n) in enumerate(zip(
            np.array(self.fr_cri['musig'])[:, 0],
            np.array(self.fr_cri['musig'])[:, 1],
            np.array(self.fr_cri['n'])))]
        self.fr_cri['fake_sample'] = cri_sample
        # AI state data
        self.corrs = []
        self.cvs = []
        # others
        self.verify_dict = {}
        if os.path.isdir(self.path):
            print('Spikes.__init__(): data directory already exists')
        else:
            os.mkdir(self.path)
            print('Spikes.__init__(): data directory created')
        self.set_labels()

    def set_labels(self):
        self.populations = ['L2/3 Exc', 'L2/3 PV', 'L2/3 SOM', 'L2/3 VIP',
                       'L4 Exc', 'L4 PV', 'L4 SOM',
                       'L5 Exc', 'L5 PV', 'L5 SOM',
                       'L6 Exc', 'L6 PV', 'L6 SOM']

        self.subtypes = ['Exc', 'PV', 'SOM', 'VIP', 'Exc', 'PV', 'SOM', 'Exc', 'PV', 'SOM', 'Exc', 'PV', 'SOM']

        self.positions = [0, 1, 2, 3, 0, 1, 2, 0, 1, 2, 0, 1, 2]

        self.layer_by_pops = [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]

        self.layer_labels = ['L2/3', 'L4', 'L5', 'L6']

        self.pops_by_layer = [[0, 1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]

        self.colors = [(68/255,119/255,170/255), (238/255,102/255,119/255), (34/255,136/255,51/255), (204/255,187/255,68/255),
                        (68/255,119/255,170/255), (238/255,102/255,119/255), (34/255,136/255,51/255),
                        (68/255,119/255,170/255), (238/255,102/255,119/255), (34/255,136/255,51/255),
                        (68/255,119/255,170/255), (238/255,102/255,119/255), (34/255,136/255,51/255)]

        self.colors_10 = [(51/255,34/255,136/255), (136/255,204/255,238/255), (68/255,170/255,153/255),
                          (17/255,119/255,51/255), (153/255,153/255,51/255), (221/255,204/255,119/255),
                          (204/255,102/255,119/255), (136/255,34/255,85/255), (170/255,68/255,153/255),
                          (187/255,187/255,187/255)]

    def read_name(self):
        # Import filenames
        for file in os.listdir(self.path):
            # print(self.device_types)
            for i, name in enumerate(self.device_types):
                if name not in self.devices:
                    self.devices[name] = []
                # print(name)
                if file.startswith(name):
                    # print(file)
                    device = file.split('-')[0] + '-' + file.split('-')[1]
                    if device not in self.devices[name]:
                        self.devices[name].append(device)
        for key, value in self.devices.items():
            self.devices[key] = sorted(self.devices[key])
            # print(key, value)

        # Import GIDs
        gidfile = open(os.path.join(self.path, 'population_GIDs.dat'), 'r')
        for l in gidfile:
            a = l.split()
            self.gids.append([int(a[0]), int(a[1])])

        self.verify_collect('self.devices=\n{}\n'.format(self.devices), 'reading')

    def get_threadfns(self, device_type, device):
        all_fns = os.listdir(self.path)
        thread_fns = [
            all_fns[x] for x in list(range(len(all_fns)))
            if all_fns[x].startswith(device_type) and
               (all_fns[x].split('-')[0] + '-' + all_fns[x].split('-')[1]) == device
        ]
        return thread_fns

    def load(self):
        self.df = {}
        self.df_columns={
            'spike_detector': ['id', 'time', 'population', 'layer'],
            'weight_recorder': ['source', 'target', 'time', 'weight', 'population', 'layer']
        }
        # loop device type
        for i, dev_type in enumerate(self.device_types):
            if dev_type not in self.devices \
                or len(self.devices[dev_type]) == 0 \
                or len(self.gids) == 0:
                print('load(): {} devices or gids not found'.format(i))
                continue
            if dev_type not in self.df:
                self.df[dev_type] = pd.DataFrame(columns=self.df_columns[dev_type])
            # loop device/population
            for j, device in enumerate(self.devices[dev_type]):
                data, thread_fns = [], self.get_threadfns(dev_type, device)
                for k, fn in enumerate(thread_fns):
                    data_thread, abs_fn = [], os.path.join(self.path, fn)
                    if os.path.isfile(abs_fn) and os.stat(abs_fn).st_size > 0:
                        data_thread = np.loadtxt(abs_fn)
                        # if data_thread.ndim < 2:
                        #     print(abs_fn, data_thread.shape)
                    if isinstance(data_thread, np.ndarray) and len(data_thread) > 0:
                        if data_thread.ndim == 1:
                            data.append(np.array([data_thread]))
                        elif data_thread.ndim == 2:
                            data.append(data_thread)
                        else:
                            print('bug: data_thread.ndim != 1 or 2')
                if len(data) > 0:
                    data = np.concatenate(data)
                if isinstance(data, np.ndarray) and data.ndim == 2:
                    # data = data[np.argsort(data[:, 1])]  # time consuming
                    lyr = self.layer_by_pops[j%13]
                    data = np.concatenate((data, np.full((len(data), 1), j),
                        np.full((len(data), 1), lyr)), axis=1)
                    df = pd.DataFrame(data=data, columns=self.df_columns[dev_type])
                    self.df[dev_type] = self.df[dev_type].append(df, ignore_index=True)

    def plot_weight(self, dev_type='weight_recorder', src_pop='L2/3 Exc', trg_pop='L2/3 Exc', n_pre=40, n_post=10, bw=10.):
        t0 = time.time()
        src = self.populations.index(src_pop)
        trg = self.populations.index(trg_pop)
        df = self.df[dev_type]
        src_ids = np.array(list(set(df['source'])))
        src_ids = src_ids[(self.gids[src][0]<=src_ids)&(src_ids<=self.gids[src][-1])].tolist()
        if len(src_ids) == 0:
            return
        src_ids = sample(src_ids, min(n_pre, len(src_ids)))
        # print('src_ids = {}'.format(list(map(int, src_ids))))
        fig, axs = plt.subplots(2, 1, figsize=(16, 16), sharex=True)
        axs[0].set_ylabel('weight (pA)')
        axs[1].set_ylabel('average weight (pA)')
        # samples
        for i, src_id in enumerate(src_ids):
            trg_ids = np.array(list(set(df[df['source']==src_id].target.values)))
            trg_ids = trg_ids[(self.gids[trg][0]<=trg_ids)&(trg_ids<=self.gids[trg][-1])].tolist()
            if len(trg_ids) == 0:
                continue
            trg_ids = sample(trg_ids, min(n_post, len(trg_ids)))
            # if i < 10:
            #     print('trg_ids={}'.format(list(map(int, trg_ids))))
            for j, trg_id in enumerate(trg_ids):
                tmp = df[(df['source']==src_id)&(df['target']==trg_id)]
                data = np.array([tmp.time.values, tmp.weight.values]).T
                data = data[np.argsort(data[:, 0])]
                axs[0].scatter(data[:, 0], data[:, 1], s=4, color='black')
                # if i < 10 and j == 0:
                #     ax.plot(data[:, 0], data[:, 1], color=self.colors_10[i])
        # mean
        mws, mws_pop, bins = [], [], np.arange(0., max(df.time.values), bw)
        df_roi = df[(self.gids[src][0]<=df.source)&(df.source<=self.gids[src][-1]) \
            &(self.gids[trg][0]<=df.target)&(df.target<=self.gids[trg][-1])]
        # print(df_roi)
        for i, bhead in enumerate(bins):
            ws = df_roi[(df_roi.time>=bhead)&(df_roi.time<bhead+bw)].weight.values
            # print('len(ws) = {}'.format(len(ws)))
            mws.append(np.mean(ws))
            mws_pop.append(np.sum(ws)/(self.gids[src][-1]-self.gids[src][0]+1))
        axs[1].plot(bins+bw/2, mws, color='g')
        axs[1].plot(bins+bw/2, mws_pop, color='b')
        t1 = time.time()
        plt.xlabel('time (ms)')
        plt.tight_layout()
        plt.savefig('weight_{}->{}.png'.format(src_pop.replace(' ','').replace('/',''), trg_pop.replace(' ', '').replace('/','')))
        plt.close()
        # print('time plot_weight(): {:.4f}, {:.4f}'.format(t1-t0, time.time()-t1))

    def seg_musig(self, ids, df, shead, stail, bw, id_str, return_dict):
        # multiprocessing
        if type(id_str) != str:
            print('seg_musig(): id_str must be str!')
            return
        print('seg_musig() proc {} start'.format(id_str))
        t0 = time.time()
        # get segment data
        data_seg = df[(shead<=df[:, 1])&(df[:, 1]<stail)]
        # data_seg = df[(shead<=df.time)&(df.time<stail)][['target', 'time', 'weight']].values
        # bins
        bheads, btails = np.arange(shead, stail, bw), np.arange(shead, stail, bw) + bw
        # caches
        means, vars = [], []
        # by neuron
        for j, id in enumerate(ids):
            data = data_seg[data_seg[:, 0]==id]
            ts, ws = data[:, 1], data[:, 2]
            # sums of weights*taus of nueron j in this segment
            wsums = np.full(len(bheads), 0.)
            # by bin
            for k, (bhead, btail) in enumerate(zip(bheads, btails)):
                # weights in bin k
                ws_bin = ws[(bhead<=ts)&(ts<btail)]
                # sum of weights*taus in bin k
                if len(ws_bin) > 0:
                    taus = np.full(ws_bin.shape, 0.)
                    taus[ws_bin>=0.] = self.net_dict['neuron_params']['tau_syn_ex']
                    taus[ws_bin<0.] = self.net_dict['neuron_params']['tau_syn_in']
                    wsums[k] = np.sum(np.multiply(taus, ws_bin))
                    # if j < 5:
                    #     print('taus={}'.format(taus))
                    #     print('ws_bin={}'.format(ws_bin))
            means.append(np.mean(wsums))
            vars.append(np.var(wsums))
            # if j < 100:
            #     self.verify_collect('{:04d}: {:2d},{:2d},{:2d}\n'.format(id, t3-t2,t4-t3,t5-t4), 'musig')
        return_dict[id_str] = (means, vars)
        print('seg_musig() proc {} end, total time = {} ms'.format(id_str, round((time.time()-t0)*1000)))

    def stationary_musig(self, begin, endin, sw=100., bw=1.0, pop_name='L2/3 Exc', verify=False, mp=False):
        t0 = time.time()
        # index of population
        pop = self.populations.index(pop_name)
        # ids of this population
        ids = np.arange(self.gids[pop][0], self.gids[pop][-1]+1)
        # segments
        sheads, stails = np.arange(begin, endin, sw), np.arange(begin, endin, sw) + sw
        # get data
        dev_type='weight_recorder'
        df = self.df[dev_type][(begin<=self.df[dev_type].time)&(self.df[dev_type].time<endin)][['target', 'time', 'weight']].values
        # df = self.df[dev_type][(begin<=self.df[dev_type].time)&(self.df[dev_type].time<endin)]
        # self.verify_collect('sources = \n{}'.format(list(set(self.df[dev_type].source.values))), 'musig')
        # caches
        means_all, vars_all, pvals_mean, pvals_var = [], [], [], []
        # multiprocessing or not
        if mp:
            re_dict, procs = Manager().dict(), []
        else:
            re_dict = {}
        # by segment
        for i, (shead, stail) in enumerate(zip(sheads, stails)):
            if mp:
                proc = Process(target=self.seg_musig, args=(ids, df, shead, stail, bw, str(shead), re_dict))
                procs.append(proc)
                proc.start()
            else:
                self.seg_musig(ids, df, shead, stail, bw, str(shead), re_dict)

        # join multiprocessing for correlation
        if mp:
            for proc in procs:
                proc.join()

        # collect data
        for i, (shead, stail) in enumerate(zip(sheads, stails)):
            means, vars = re_dict[str(shead)]
            # t test
            if i > 0:
                stat, pval_mean = stats.ttest_ind(means_all[-1], means)
                stat, pval_var = stats.ttest_ind(vars_all[-1], vars)
                pvals_mean.append(pval_mean)
                pvals_var.append(pval_var)
            # collect
            means_all.append(means)
            vars_all.append(vars)

        # plot
        fig, axs = plt.subplots(2, 1, figsize=(16, 16))
        plt.xlabel('time (ms)')
        xs = sheads + sw/2
        flierprops = dict(marker='+', markerfacecolor='k', markersize=5,
                  linestyle='none', markeredgecolor='k')
        axs[0].boxplot(means_all, positions=xs, widths=sw/4, sym='.', flierprops=flierprops)
        axs[1].boxplot(vars_all, positions=xs, widths=sw/4, sym='.', flierprops=flierprops)

        # plot significant p-values
        mbot, mtop = axs[0].get_ylim()
        vbot, vtop = axs[1].get_ylim()
        for i, stail in enumerate(stails[:-1]):
            mean_seg, var_seg = np.mean(means_all[i]), np.mean(vars_all[i])
            axs[0].text(stail-sw/2, mbot+(0.5+0.01*(i%10))*(mtop-mbot), \
                '{:.3f}'.format(mean_seg), horizontalalignment='center', fontsize=12, color='b')
            axs[1].text(stail-sw/2, vbot+(0.5+0.01*(i%10))*(vtop-vbot), \
                '{:.0f}'.format(var_seg), horizontalalignment='center', fontsize=12, color='b')
            if pvals_mean[i] <= 0.05:
                axs[0].text(stail, mbot+(0.9+0.01*(i%10))*(mtop-mbot), \
                    '{:.3f}'.format(pvals_mean[i]), horizontalalignment='center', fontsize=12, color='g')
            if pvals_var[i] <= 0.05:
                axs[1].text(stail, vbot+(0.9+0.01*(i%10))*(vtop-vbot), \
                    '{:.3f}'.format(pvals_var[i]), horizontalalignment='center', fontsize=12, color='g')

        xticks = np.append(sheads, endin).astype(int)
        xticklabels, tick_interval = [], int(len(sheads)/10)
        tick_interval = 1 if tick_interval == 0 else tick_interval
        for i, tick in enumerate(xticks):
            lbl = str(int(tick)) if i%(tick_interval) == 0 else ''
            xticklabels.append(lbl)
        # xticklabels = xticks.astype(str)
        axs[0].set_xticks(xticks)
        axs[0].set_xticklabels(xticklabels, rotation=30)
        axs[1].set_xticks(xticks)
        axs[1].set_xticklabels(xticklabels, rotation=30)

        axs[0].set_ylabel(r'weight mean ($ms\cdot pA$)')
        axs[1].set_ylabel(r'weight variance ($ms^{2}\cdot pA^{2}$)')
        # plt.tight_layout()
        plt.savefig(os.path.join(self.path, 'stationary_musig.png'))
        plt.close()
        print('stationary_musig() running time = {:.3f}'.format(time.time()-t0))

    def get_data(self, begin, end, dev_type='spike_detector'):
        t0 = time.time()
        data_ret = []
        for i, pop in enumerate(self.populations):
            df = self.df[dev_type]
            columns = df.columns.tolist()
            data = df[(df.population==i)&(df.time>=begin)&(df.time<end)][columns[:-2]].values
            if 'time' in columns:
                data = data[np.argsort(data[:, columns.index('time')])]
            data_ret.append(data)
        # print('get_data(): {}'.format(time.time()-t0))
        return data_ret

    def get_ts_by_id(self, data, id):
        return_data = []
        if isinstance(data, np.ndarray) and data.ndim == 2:
            ids = data[:, 0]
            ts = data[:, 1]
            return_data = ts[ids == id]
        return return_data

    def verify_collect(self, in_str, tag):
        if tag in self.verify_dict.keys():
            self.verify_dict[tag] += in_str
        else:
            self.verify_dict[tag] = in_str

    def verify_print(self, path=None):
        if os.path.isdir(path):
            path_flg = True
        else:
            path_flg = False
        for key, value in self.verify_dict.items():
            if path_flg:
                fpath = os.path.join(path, 'verify-{}.txt'.format(key))
            else:
                fpath = 'verify-{}.txt'.format(key)
            with open(fpath, 'w') as f:
                f.write(value)
                f.close()

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
    t0 = time.time()
    data = spikes.get_data(begin, end)
    gids = spikes.gids
    rates_averaged_all = []
    rates_std_all = []
    for h in list(range(len(data))):
        if len(data[h]) > 0:
            ids = data[h][:, 0]
            ids = ids.astype(int)
            count_of_n = np.bincount(ids) # bin count from 0
            count_of_ids = count_of_n[gids[h][0]:gids[h][1]+1] # get data with corresponding indices
            rate_each_n = count_of_ids * 1000. / (end - begin)
            rate_averaged = np.mean(rate_each_n)
            rate_std = np.std(rate_each_n)
            rates_averaged_all.append(float('%.3f' % rate_averaged))
            rates_std_all.append(float('%.3f' % rate_std))
            np.save(os.path.join(spikes.path, ('rate' + str(h) + '.npy')), rate_each_n)
            if h == 2:
                spikes.verify_collect('data[{}]=\n{}\n'.format(h, data[h]), 'fr')
                spikes.verify_collect('gids[{}]=\n{}\n'.format(h, gids[h]), 'fr')
                spikes.verify_collect('len(ids)={}, values=\n'.format(len(ids)), 'fr')
                for n in ids:
                    spikes.verify_collect('{} '.format(n), 'fr')
                spikes.verify_collect('\n', 'fr')
                spikes.verify_collect('count_of_n(tail)=\n{}\n'.format(count_of_n[-len(count_of_ids):]), 'fr')
                spikes.verify_collect('count_of_ids=\n{}\n'.format(count_of_ids), 'fr')
                spikes.verify_collect('np.sum(count_of_ids)=\n{}\n'.format(np.sum(count_of_ids)), 'fr')
                spikes.verify_collect('len(count_of_n)=\n{}\n'.format(len(count_of_n)), 'fr')
                spikes.verify_collect('len(count_of_ids)=\n{}\n'.format(len(count_of_ids)), 'fr')
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
    spikes.fr_result['musig'] = np.array([rates_averaged_all, rates_std_all])
    fake_sample = [spikes.fake_sample(mu, sig, gids[i][1] - gids[i][0] + 1, estimate=False)
            for i, (mu, sig) in enumerate(zip(
            rates_averaged_all,
            rates_std_all
            ))]
    spikes.fr_result['pvalues'] = \
        [stats.ttest_ind(fake_sample[i], spikes.fr_cri['fake_sample'][i], equal_var=False)[1]
            for i in range(len(fake_sample))]
    print('fire_rate() running time = {:.3f} s'.format(time.time()-t0))
    return rates_averaged_all, rates_std_all


def plot_raster(spikes, begin, end, figsize=(12, 9), fontsize=30.):

    t0 = time.time()
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

    fig, ax = plt.subplots(figsize=figsize)
    for i in list(range(len(data))):
        if len(data[i]) > 0:
            times = data[i][:, 1]
            neurons = np.abs(data[i][:, 0] - highest_gid) + 1
            plt.plot(times, neurons, '.', color=spikes.colors[i])

    # legend
    for i in range(4):
        plt.scatter([], [], s=100, label=spikes.subtypes[i], color=spikes.colors[i])
    plt.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.2), fontsize=fontsize)

    # set top and right frames invisible
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # reset y limits to contain legend
    bottom, top = ax.get_ylim()

    plt.xlabel('time (ms)', fontsize=fontsize)
    plt.xticks(np.arange(begin, end + 1.0, (end - begin)/4.0), fontsize=fontsize)
    plt.yticks(
        [L23_label_pos, L4_label_pos, L5_label_pos, L6_label_pos],
        ylabels, rotation=10, fontsize=fontsize
        )
    fig.tight_layout()
    plt.savefig(os.path.join(spikes.path, 'raster_plot.png'), dpi=300)
    plt.close()
    print('plot_raster() running time = {:.3f} s'.format(time.time()-t0))

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

def label_diff(ax, xs, ys, pvalues, bw):
    props = {'connectionstyle':'bar','arrowstyle':'-',\
                 'shrinkA':20,'shrinkB':20,'linewidth':2}
    for i, p in enumerate(pvalues):
        if p <= 0.05:
            ax.annotate('*', xy=(xs[i], min(30., ys[i]+7)),
                horizontalalignment='center', verticalalignment='center', fontsize=20)
            # ax.annotate('*', xy=(xs[i], 1.1*ys[i]), zorder=10, horizontalalignment='center', fontsize=20)
            ax.hlines(min(29.5, ys[i]+6.5), xs[i]-bw/2, xs[i]+bw/2, colors='k')
            ax.vlines([xs[i]-bw/2, xs[i]+bw/2], min(29., ys[i]+6), min(29.5, ys[i]+6.5), colors='k')
            # ax.annotate('', xy=(xs[i], 1.1*ys[i]), xytext=(xs[i], 1.1*ys[i]), arrowprops=props)

def do_bars(data, cri, data_n, cri_n, path, title, colors, ylbl, use_SE=True,
        # figsize=(16, 12), ylims=(-1., 31.), fontsize=30):
        figsize=(12, 9), ylims=(-1., 31.), fontsize=30, stattest=True, pvalues=None):
    layers = ['L2/3', 'L4', 'L5', 'L6']
    legends = ['Exc', 'PV', 'SOM', 'VIP']
    err_data = data[1, :]/np.sqrt(data_n) if use_SE else data[1, :]
    err_cri = cri[0, :]/np.sqrt(cri_n) if use_SE else cri[0, :]

    # bars
    xs = np.arange(13)  # the label locations
    bw = 0.3  # the width of the bars
    fig, ax = plt.subplots(figsize=figsize)
    rects1 = ax.bar(xs - bw/2, data[0, :], bw, yerr=err_data, color=colors, edgecolor=colors, linewidth=2)
    rects2 = ax.bar(xs + bw/2, cri[0, :], bw, yerr=err_cri, fill=False, edgecolor=colors, linewidth=2, hatch='//')

    # legends
    for i in range(4):
        plt.scatter([], [], s=100, marker='s', label=legends[i], color=colors[i])
    legend = plt.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.2), fontsize=fontsize)

    # plot settings
    ax.set_ylabel(ylbl, fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_xlim(0-1.5*bw, 12+1.5*bw)
    ax.set_ylim(ylims[0], ylims[1])

    # vline for layer boundary
    for x in [3.5, 6.5, 9.5]:
        ax.axvline(x, 0, 1, color='k')

    # mark layer
    for i, x in enumerate([1.5, 5, 8, 11]):
        ax.text(x, ax.get_ylim()[0]-.5, layers[i], horizontalalignment='center',
            verticalalignment='top', fontsize=fontsize)

    # mark statistical compareing
    if isinstance(pvalues, list):
        ys = [max(data[0, i], cri[0, i]) for i in range(len(data[0, :]))]
        label_diff(ax, xs, ys, pvalues, bw)

    fig.tight_layout()
    plt.savefig(os.path.join(path, 'bars_' + title + '.png'))
    plt.close()

def fr_plot(spikes):
    t0 = time.time()
    rates, n_pops = [], []
    for i in range(len(spikes.populations)):
        fpath = os.path.join(spikes.path, ('rate' + str(i) + '.npy'))
        n_pops.append(spikes.gids[i][1]-spikes.gids[i][0]+1)
        if os.path.isfile(fpath):
            rates.append(np.load(fpath))
    # do_boxplot(rates[::-1], spikes.fr_cri_qaurters[::-1], spikes.path, 'fr', spikes.colors[::-1], spikes.subtypes[::-1], 'firing rate (spike/s)', xlims=(-1.0, 60.0))
    do_bars(spikes.fr_result['musig'], np.array(spikes.fr_cri['musig']).T, n_pops,
        spikes.fr_cri['n'], spikes.path, 'fr', spikes.colors, 'spikes/s')
        # do_bars(spikes.fr_result['musig'], np.array(spikes.fr_cri['musig']).T, n_pops,
        #     spikes.fr_cri['n'], spikes.path, 'fr', spikes.colors, 'spikes/s', pvalues=spikes.fr_result['pvalues'])
    print('fr_plot() running time = {:.3f} s'.format(time.time()-t0))

'''
Other analysis
'''
# Ground state calculation
def gs_analysis(spikes, begin, end, bw=10, seg_len=5000.0, n_sample=140, n_spk_thr=4):
    t0 = time.time()
    # setup
    segs = list(zip(np.arange(begin, end, seg_len), np.arange(begin, end, seg_len) + seg_len))
    ai = np.full((4, 2), np.nan)
    ai_n = ''
    ai_xcpt = ''

    # multiprocessing
    return_dict, procs = Manager().dict(), []

    # prepare calculation
    df = copy.deepcopy(spikes.df['spike_detector'])
    # df = copy.deepcopy(spikes.df[])
    df_corr = pd.DataFrame(columns=['corr', 'segment', 'layer'])
    df_cv = pd.DataFrame(columns=['cv', 'segment', 'layer'])

    # filter by spike n
    set_keep = set(df.id)
    for i, (seg_head, seg_tail) in enumerate(segs):
        df_seg = df[(df.time >= seg_head) & (df.time < seg_tail)]
        val_cnts = df_seg.id.value_counts()
        set_keep = set_keep & set(val_cnts[val_cnts >= n_spk_thr].index)
    df = df[df.id.isin(set_keep)]

    # sample
    df_sp = pd.DataFrame(columns=['id', 'time', 'population', 'layer'])
    for k in range(4):
        set_total = set(df[df.layer==k].id)
        if len(set_total) >= n_sample:
            set_sampled = sample(set_total, n_sample)
            df_sp = df_sp.append(df[df.id.isin(set_sampled)], ignore_index=True)

    # loop by segment
    for i, (seg_head, seg_tail) in enumerate(segs):
        df_seg = df_sp[(df_sp.time >= seg_head) & (df_sp.time < seg_tail)]
        # loop by layer
        for k in range(4):
            hists, cvs = [], []
            df_lyr = df_seg[df_seg.layer==k]
            for id in list(set(df_lyr.id)):
                ts = df_lyr[df_lyr.id == id].time.values
                # print(ts)
                if len(ts) > 0:
                    hist_bins = np.arange(seg_head, seg_tail + bw, bw)
                    hist = np.histogram(ts, bins=hist_bins)[0]
                    # correlation
                    if np.all(hist == hist[0]):
                        # report monotone histogram
                        ai_xcpt += 'seg {}, layer {}, id {}: monotone histogram\n'.format(i, k, id)
                    else:
                        hists.append(hist)
                    # CV ISI
                    if len(ts) < n_spk_thr:
                        # report ISI n < threshould
                        ai_xcpt += 'seg {}, layer {}, id {}: CV-ISI n<threshold\n'.format(i, k, id)
                    else:
                        isi = np.diff(ts)
                        cvs.append(np.std(isi) / np.mean(isi))
                else:
                    # report n = 0
                    ai_xcpt += 'seg {}, layer {}, id {}: CV-ISI n=0\n'.format(i, k, id)

            # set multiprocessing for correlation
            proc = Process(target=get_corr, args=(hists, None, int(i*4 + k), return_dict))
            procs.append(proc)
            proc.start()

            # get cv
            df_cv = df_cv.append({'cv': np.mean(cvs), 'segment': i, 'layer': k}, ignore_index=True)
            ai_n += '{}, {}\n'.format(len(hists), len(cvs))
        ai_n += '\n'

    # join multiprocessing for correlation
    for proc in procs:
        proc.join()

    # collect data
    for k in range(4):
        for i, seg_head in enumerate(segs):
            # get segment-layer correlation
            try:
                corr_proc = np.mean(return_dict[str(int(i*4 + k))])
            except KeyError:
                corr_proc = np.nan
            df_corr = df_corr.append({'corr': corr_proc, 'segment': i, 'layer': k}, ignore_index=True)
        corrs_lyr = df_corr[df_corr['layer']==k]['corr'].values
        cvs_lyr = df_cv[df_cv['layer']==k]['cv'].values
        corr_avg = np.mean(corrs_lyr)
        cv_avg = np.mean(cvs_lyr)
        spikes.corrs.append(corrs_lyr)
        spikes.cvs.append(cvs_lyr)
        ai[k, :] = [corr_avg, cv_avg]

    np.savetxt(os.path.join(spikes.path, 'ai.dat'), ai, fmt='%.10f', delimiter=',')
    f1 = open(os.path.join(spikes.path, 'ai_n.dat'), 'w')
    f1.write(ai_n)
    f1.close()
    if len(ai_xcpt) > 0:
        f2 = open(os.path.join(spikes.path, 'ai_xcpt.dat'), 'w')
        f2.write(ai_xcpt)
        f2.close()
    print('gs_analysis() running time = {:.3f} s'.format(time.time()-t0))

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
                    cnts_by_stim = np.array([len(np.where((ts >= stim + b*bin_w) & (ts < stim + (b + 1)*bin_w))[0]) for stim in stims])
                    pooled_std = np.sqrt((np.var(cnts_by_stim[0::2], ddof=1) + np.var(cnts_by_stim[1::2], ddof=1))/2)
                    if pooled_std != 0: # set some more criteria, e.g. spike number must > 0 ...?
                        if raw is True:
                            SI_by_neuronbin[j, b] = (np.mean(cnts_by_stim[0::2]) - np.mean(cnts_by_stim[1::2]))/pooled_std
                        else:
                            SI_by_neuronbin[j, b] = 2*(np.abs(np.mean(cnts_by_stim[0::2]) - np.mean(cnts_by_stim[1::2])))/pooled_std    # multiplied by 2 to compare!

            # Calculate population selectivity and plot
            len_q = int(len(SI_by_neuronbin)/4)
            idx_lyr = spikes.layer_by_pops[i]
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
            axs[idx_lyr].text(-25.0, np.mean(ylims), spikes.layer_by_pops[idx_lyr])
    for ax in axs:
        ax.plot([0.0, bin_w*n_bin], [0.0, 0.0], color='k', linestyle='--')
    plt.xticks(xs[::2], labels=xs[::2].astype(str))
    plt.savefig(os.path.join(spikes.path, 'selectivity-raw={}.png'.format(raw)))
    plt.close()


# responses to transient/thalamic input
def response(spikes, begin, stims, window, bw=1.0, exportplot=False,
        interpol='quadratic', ylims=(-0.05, 0.45), ptitle='STPs',
        figsize=(8, 4), font_puff=5/4):
    matplotlib.rcParams['font.size'] = matplotlib.rcParams['font.size']*font_puff
    t0 = time.time()
    n_stim = len(stims)
    if len(stims) > 1:
        interval = stims[1] - stims[0]
    data_all = spikes.get_data(begin, begin+n_stim*interval)
    # response data to file
    f = open(os.path.join(spikes.path, 'sf.dat'), 'w')

    # plot settings
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    plt.title(ptitle)
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
            ts_stim = ts[(ts >= win_begin) & (ts < win_end)]
            rsp_amp.append(len(ts_stim))
            if len(ts_stim) >= 3:
                std = np.std(ts_stim)
                rsp_spread.append(std)
            # latency by neurons
            ids_stim = ids[(ts >= win_begin) & (ts < win_end)]
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

        # plot:
        hist, bins = np.histogram(mean_ltcs, bins=np.arange(0.0, window, bw))
        xs, ys = bins[:-1], hist/np.sum(hist)
        xy_by_layer.append([xs, ys])
        # interpolate
        if isinstance(interpol, str):
            f_cubic = interpolate.interp1d(xs, ys, kind=interpol)
            xs = np.linspace(min(xs), max(xs), len(xs)*5)
            ys = f_cubic(xs)
            ys[ys<0] = 0
        ax.plot(xs, ys,
            linestyle=linestyles[a],
            linewidth=3,
            label=spikes.populations[i],
            color=colors[a])
        ax.set_xlabel('mean spike latency (ms)')
        ax.set_ylabel('fraction')
        ax.set_ylim(ylims[0], ylims[1])

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
    print('response() running time = {:.3f} s'.format(time.time()-t0))
    matplotlib.rcParams['font.size'] = matplotlib.rcParams['font.size']/font_puff

def perturb_plot_single_layer(spikes, xs, frs_all, stim_type, normalize=False,
        frs_ei=None, pops=None, targets=None, ylims=None, figsize=(4, 10)):
    fontsize_puff = 5/4
    matplotlib.rcParams['font.size'] = matplotlib.rcParams['font.size']*fontsize_puff
    # print(frs_all)
    if pops is None:
        pops = [0, 1, 2, 3]
    if ylims is None:
        frs_means = np.mean(np.array([np.array(frs_all).max(axis=1), np.array(frs_all).min(axis=1)]), axis=0)
        ylims = {
            'Exc': (frs_means[0]-0.3, frs_means[0]+0.3),
            'PV': (frs_means[1]-4., frs_means[1]+4.),
            'SOM': (frs_means[2]-2., frs_means[2]+2.),
            'VIP': (0., 1.)
        }

    # set labels
    xlbl = r'$r_{stim}$ (spikes/s)' if stim_type == 'poisson' else r'$I_{stim}$ (pA)'
    ylbl = r'normalized $r$' if normalize else r'$r$ (spikes/s)'

    # set figure
    if normalize:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        target_str = ''
        for target in targets:
            target_str += '{}, '.format(spikes.subtypes[target])
        target_str = target_str[:-2] + ' activated'
        ax.set_title(target_str, y=1.0, pad=100)
        # ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{:.1f}'))
        # yticks = ax.get_yticks()[ax.get_yticks()>=0.]
        # ax.set_yticks(yticks)
        # ax.set_yticklabels(yticks)
    else:
        fig, axs = plt.subplots(len(pops), 1, sharex=True, figsize=figsize)
        # fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=figsize)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)

    # plot
    thin, thick = 2, 4
    for i, pop in enumerate(pops):
        if pop >= len(frs_all):
            break
        # mark stimulation target
        if isinstance(targets, list):
            linewidth = thick if pop in targets else thin
            linestyle = '--' if pop in targets else '-'
        else:
            linewidth = thin
            linestyle = '-'
        # plot
        if normalize:
            ax.plot(xs, frs_all[pop], color=spikes.colors[pop], label=spikes.subtypes[pop],
                linewidth=linewidth, marker='.', markersize=10, linestyle=linestyle)
        else:
            axs[i].plot(xs, frs_all[pop], color=spikes.colors[pop],
                linewidth=linewidth, marker='.', markersize=10, linestyle=linestyle)
            if isinstance(ylims, dict):
                axs[i].set_ylim(ylims[spikes.subtypes[pop]])
            axs[0].plot([], [], label=spikes.subtypes[i], color=spikes.colors[i])
            # ax1.plot(xs, frs_all[i], color=spikes.colors[i], label=spikes.subtypes[i],
            #     linewidth=linewidth, marker='.', markersize=15, linestyle=linestyle)
            # ax2.plot(xs, frs_all[i], color=spikes.colors[i], label=spikes.subtypes[i],
            #     linewidth=linewidth, marker='.', markersize=15, linestyle=linestyle)

    # special
    if normalize:
        # plt.ylim((0., 2.))
        span = max(plt.ylim()[1] - 1., 1. - plt.ylim()[0])
        plt.ylim((1. - span, 1 + span))
        ax.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.16))
    else:
        axs[0].legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.7))
        # # spacing
        # ax_ylims = ax1.get_ylim()
        # yspan = min((ax_ylims[1] - ax_ylims[0])/2, yspan)
        # ax1.set_ylim((ax_ylims[1]-yspan, ax_ylims[1]))
        # ax2.set_ylim((ax_ylims[0], ax_ylims[0]+yspan))
        # # ticks
        # ax1.spines['bottom'].set_visible(False)
        # ax2.spines['top'].set_visible(False)
        # ax1.xaxis.tick_top()
        # ax1.tick_params(labeltop=False)
        # ax2.xaxis.tick_bottom()
        # # slanted lines
        # d = .5
        # kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
        #     linestyle="none", color='k', mec='k', mew=1, clip_on=False)
        # ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
        # ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
        # ax1.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.5))

    plt.savefig(os.path.join(spikes.path, 'perturb_single_normalize={}.png'.format(normalize)), bbox_inches='tight')
    plt.close()
    matplotlib.rcParams['font.size'] = matplotlib.rcParams['font.size']/fontsize_puff

def perturb_plot_all_layer(spikes, xs, frs_all, stim_type, normalize=False, shrink=1, frs_ei=None, targets=None):
    fig, axs = plt.subplots(4, 1, figsize=(6, 15), sharex=True, sharey=True)
    purple = (170/255, 51/255, 119/255)
    xlbl = r'$r_{stim}$ (spikes/s)' if stim_type == 'poisson' else r'$I_{stim}$ (pA)'
    plt.xlabel(xlbl)
    if normalize:
        ylbl = r'normalized $r$'
    else:
        ylbl = r'$r$ (spikes/s)'
    plt.ylabel(ylbl)
    thin, thick = 2, 4
    # plot 4-group data
    for i, gids in enumerate(spikes.gids):
        if isinstance(targets, list):
            linewidth = thick if i in targets else thin
            linestyle = '--' if i in targets else '-'
        else:
            linewidth = thin
            linestyle = '-'
        if normalize is False and spikes.subtypes[i] in ['PV', 'VIP']:
            frs = np.array(frs_all[i])/shrink
        else:
            frs = np.array(frs_all[i])
        axs[spikes.layer_by_pops[i]].plot(xs, frs, color=spikes.colors[i],
            linewidth=linewidth, marker='.', markersize=10, linestyle=linestyle)
    # plot e-i group data
    if frs_ei is not None:
        for i, frs in enumerate(frs_ei):
            axs[i].plot(xs, frs[0], color=spikes.colors[0], linewidth=thick, marker='.', markersize=10)
            if normalize:
                axs[i].plot(xs, frs[1], color=purple, linewidth=thick, marker='.', markersize=10)
            else:
                axs[i].plot(xs, frs[1]/shrink, color=purple, linewidth=thick, marker='.', markersize=10)
    # legend
    for i in range(4):
        axs[0].plot([], [], label=spikes.subtypes[i], color=spikes.colors[i])
    axs[0].legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.5))
    # layer labels
    for i in range(4):
        axs[i].text(-0.3, 0.5, spikes.layer_labels[i], horizontalalignment='center', verticalalignment='center', transform=axs[i].transAxes)
    # set ylims
    if normalize:
        plt.ylim((0.5, 1.5))
        # bottom, top = plt.ylim()
        # print(bottom, top)
        # plt.ylim((max(bottom, 0.), min(top, 2.)))
    plt.savefig(os.path.join(spikes.path, 'perturb_normalize={}.png'.format(normalize)), bbox_inches='tight')
    plt.close()

def perturb_plots(spikes, levels, frs_all, frs_all_norm, stim_type, targets):
    # perturb_plot_all_layer(spikes, levels, frs_all, stim_type, normalize=False, targets=targets)
    # perturb_plot_all_layer(spikes, levels, frs_all_norm, stim_type, normalize=True, targets=targets)
    perturb_plot_single_layer(spikes, levels, frs_all, stim_type, normalize=False, targets=targets)
    perturb_plot_single_layer(spikes, levels, frs_all_norm, stim_type, normalize=True, targets=targets)

def perturb_calc(spikes, perturb_dict):
    t0 = time.time()
    if perturb_dict['n_repeat'] <= 0:
        return
    n_repeat, duration = perturb_dict['n_repeat'], perturb_dict['duration']
    frs_all, frs_all_norm = [], []
    # spksum_by_pop, ns_by_pop = [], []
    # fr calculation
    for i, gids in enumerate(spikes.gids):
        # spksum_by_stim = []
        frs, n_pop = [], gids[1] - gids[0] + 1
        # by stim. level
        for j, (offset, starts) in enumerate(perturb_dict['starts'].items()):
            spksum, fr = 0, np.nan
            # by trial
            for k, start in enumerate(starts):
                data = spikes.get_data(start+duration/2, start+duration)[i]
                if type(data) != np.ndarray or data.ndim != 2:
                    spksum = np.nan
                    break
                spksum += len(data)
            # spksum_by_stim.append(spksum)
            fr = 1000*(spksum/n_pop)/(n_repeat*(duration/2))
            frs.append(fr)
        frs_all.append(frs)
        norm = np.array(frs)/frs[0] if frs[0] != 0 else np.full(len(perturb_dict['levels']), np.nan)
        frs_all_norm.append(norm)
        # spksum_by_pop.append(spksum_by_stim)
        # ns_by_pop.append(n_pop)
    # print(perturb_dict)
    # print(frs_all)
    # print(frs_all_norm)
    perturb_plots(spikes, perturb_dict['levels'], frs_all, frs_all_norm, perturb_dict['type'], perturb_dict['targets'])
    # spksum_by_pop = np.array(spksum_by_pop)
    # ns_by_pop = np.array(ns_by_pop)
    # spk_sum_by_ei, n_sum_by_ei = [], []
    # for i in range(4):
    #     idxs = np.where(np.array(spikes.layer_by_pops) == i)[0]
    #     spk_sum_by_ei.append([spksum_by_pop[idxs[0], :]/ns_by_pop[idxs[0]],
    #     np.sum(spksum_by_pop[idxs[1:], :]/np.sum(ns_by_pop[idxs[1:]]), axis=0)])
    # e-i group data
    # frs_by_ei = (1000/(duration/2))*np.array(spk_sum_by_ei)/n_repeat
    # frs_by_ei_norm = copy.deepcopy(frs_by_ei)
    # for i, frs_lyr in enumerate(frs_by_ei_norm):
    #     for j, frs in enumerate(frs_lyr):
    #         frs_by_ei_norm[i, j] /= frs_by_ei_norm[i, j, 0]
    # perturb_plot(spikes, perturb_dict['levels'], frs_all, stim_type, normalize=False, frs_ei=frs_by_ei)
    # perturb_plot(spikes, perturb_dict['levels'], frs_all_norm, stim_type, normalize=True, frs_ei=frs_by_ei_norm)
    print('perturb_calc() running time = {} s'.format(time.time()-t0))
