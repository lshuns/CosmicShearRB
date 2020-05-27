# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-05-23 11:02:03
# @Last Modified by:   lshuns
# @Last Modified time: 2020-05-27 21:55:15

"""
make histogram plot of total weight distribution from different samples.
    plot in 5 tomographic bins
"""

import numpy as np
import pandas as pd
import feather
import matplotlib as mpl
import matplotlib.pyplot as plt


# +++++++++++++++ plot information
# input directory
inpathF_split = "/disks/shear15/ssli/KV450/split/all_tomo"
inpathP_less3 = "_T_B_less3.feather"
inpathP_greater3 = "_T_B_greater3.feather"
#
inpathFs = [inpathF_split, inpathF_split]
inpathPs = [inpathP_less3, inpathP_greater3]

# half-half cut
halfline = [3.091, 3.182, 3.273, 3.273, 3.182]

# output directory
outfile_png = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/hist_TB.png"
outfile_pdf = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/hist_TB.pdf"


# plot related
LABELS = [r'$T_{\rm B} \leq 3$', r'$T_{\rm B} > 3$']
COLORS = ['red', 'blue']
NBINS = [10, 10]
XLABEL = r"$T_{\rm B}$"        
YLABEL = "Weighted counts"
XLIM = [0.5, 6.5]
YLIM = [0, 1.5]
DENSITY = True
HISTTYPE = 'step'

# +++++++++++++++ Main code
# general settings for plot
mpl.use('Agg')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
plt.rc('font', size=12)

fig, axes = plt.subplots(nrows=2, ncols=3, sharey=True)
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)

# index for tomo bin
ibins = 1
# for legend
handles = []
for i in range(2):
    for j in range(3):
        if ibins <= 5:
            wgratio = 0
            for k in range(len(inpathFs)):
                inpath = inpathFs[k] + str(ibins) + inpathPs[k]
                data = feather.read_dataframe(inpath)

                TB = data['T_B'].values
                wg = data['recal_weight'].values
                print("Data loaded from", inpath)

                # plot
                CR = COLORS[k]
                NB = NBINS[k]
                LAB = LABELS[k]

                _, _, p_tmp = axes[i, j].hist(x=TB, bins=NB, density=DENSITY, weights=wg, color=CR, label=LAB, histtype=HISTTYPE)
                axes[i, j].axvline(x=halfline[int(ibins-1)], ymin=YLIM[0], ymax=YLIM[1], ls='--', lw=0.5, color='gray')
                axes[i, j].text(halfline[int(ibins-1)]+0.1, YLIM[1]*0.85, '{:0.3f}'.format(halfline[int(ibins-1)]), fontsize=8)


                # weight ratio
                if wgratio == 0:
                    wgratio = np.sum(wg)
                else:
                    wgratio = np.sum(wg)/wgratio

                # for legend
                if ibins == 1:
                    handles.append(p_tmp[0])

            axes[i, j].set_xlim(XLIM[0], XLIM[1])
            axes[i, j].set_ylim(YLIM[0], YLIM[1])
            axes[i, j].set_yticks([0.0, 0.3, 0.6, 0.9, 1.2])

            # bin label
            label = 'Bin' + str(ibins)
            x = XLIM[0] + 0.65*(XLIM[1]-XLIM[0])
            y = YLIM[0] + 0.7*(YLIM[1]-YLIM[0])
            axes[i,j].text(x, y, label)


            # wgratio label
            y2 = YLIM[0] + 0.6*(YLIM[1]-YLIM[0])
            axes[i,j].text(x, y2, '[{:0.2f}]'.format(wgratio), fontsize=10)

            if j == 0:
                axes[i, j].set_ylabel(YLABEL)
            if i == 1:
                axes[i, j].set_xlabel(XLABEL)

        else:
            axes[i, j].axis('off')
        ibins += 1

fig.legend(handles, LABELS, loc = 'upper right', bbox_to_anchor=(0.85, 0.35), frameon=False)
# fig.text(0.67, 0.18, r'[${\rm wg}_{\rm blue}/{\rm wg}_{\rm red}$]')
fig.text(0.67, 0.18, '[weight ratio]')


plt.savefig(outfile_png, dpi=300)
print("Plot saved in", outfile_png)

plt.savefig(outfile_pdf, dpi=300)
print("Plot saved in", outfile_pdf)



