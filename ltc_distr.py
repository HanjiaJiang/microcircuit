import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib
matplotlib.rcParams['font.size'] = 20.0

if __name__ == "__main__":
    def plot_ltc(ltc_avg, mem):
        # to-do: plotting to be modulized
        labels = ['L2/3 Exc', 'L4 Exc', 'L5 Exc', 'L6 Exc']
        colors = ['hotpink', 'dodgerblue', 'black', 'black']
        linestyles = ['solid', 'solid', 'solid', 'dashed']
        fig = plt.figure(figsize=(8, 4), constrained_layout=True)
        ax = fig.add_subplot(111)
        for i, lyr in enumerate(ltc_avg):
            xs = lyr[0]
            ys = lyr[1]
            f_cubic = interpolate.interp1d(xs, ys, kind='cubic')
            xs = np.linspace(min(xs), max(xs), len(xs)*5)
            ys = f_cubic(xs)
            ax.plot(xs, ys,
                linestyle=linestyles[i],
                linewidth=4,
                label=labels[i],
                color=colors[i])
            ax.set_xlabel('mean spike latency (ms)')
            ax.set_ylabel('fraction')
        plt.legend()
        plt.savefig('{}_ltc-distr.png'.format(mem))

    ltc_files = sorted(sys.argv[1:])

    # tags for parameter values, e.g. 0, 1, 2, ...
    tags = []
    for file in ltc_files:
        tags.append(file.split('/')[1].split('_')[0])

    ltc_sum = None  # sum for averaging
    cnt = 0 # count for averaging
    for i, file in enumerate(ltc_files):
        print(file)
        npload = np.load(file)
        if npload.ndim != 3:
            continue
        if ltc_sum is None:
            ltc_sum = npload
        else:
            ltc_sum = np.add(ltc_sum, npload)
        cnt += 1
        # get next value
        if i+1 < len(ltc_files):
            next_tag = tags[i+1]
        else:
            next_tag = None
        # plot and reset when next value is different
        if tags[i] != next_tag:
            plot_ltc(ltc_sum/cnt, tags[i])
            ltc_sum = None
            cnt = 0
