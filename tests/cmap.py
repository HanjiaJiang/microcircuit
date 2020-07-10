import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size'] = 15.0

fig, axs = plt.subplots(4, 6, figsize=(16, 12), sharex=True, sharey=True)
xs, ys = np.arange(0, 5), np.arange(0, 5)
# plt.xlim((xs[0], xs[-1]))
# plt.ylim((ys[0], ys[-1]))
xlbl, ylbl = 'x', 'y'
# plot
for c, plotvar in enumerate(['a', 'b', 'c', 'd', 'e', 'f']):
    vmin, vmax = 0.0, 1.0
    for r in range(4):
        ax = axs[r, c]
        # plot data
        data = np.random.random((5, 5))
        # im = ax.contourf(data,
        #     cmap=self.cmaps[plotvar],
        #     origin='lower',
        #     extent=self.extent,
        #     vmin=vmin,
        #     vmax=vmax)
        im = ax.imshow(data, interpolation='none',
            cmap='RdBu',
            origin='lower',
            extent=[-0.5, 4.5, -0.5, 4.5],
            vmin=vmin,
            vmax=vmax)

        # many settings
        ax.set_aspect(float((xs[-1] - xs[0])/(ys[-1] - ys[0])))

        # title
        if r == 0:
            ax.set_title(plotvar + '\n ')

        # xlabel
        if r == 3:
            ax.set_xlabel(xlbl)
            cbar = fig.colorbar(im, ax=axs[:, c], orientation='horizontal', shrink=0.8, aspect=10)

        # ylabel
        if c == 0:
            ax.set_ylabel(ylbl)
            ax.text(-0.75, 0.5, 'L', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

fig.savefig('cmap.png', bbox_inches='tight')
plt.close()
