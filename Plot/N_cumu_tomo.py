# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-06-09 18:54:52
# @Last Modified by:   lshuns
# @Last Modified time: 2020-06-10 23:35:46

"""
make plot of weight cumulative distribution.
    plot in 5 tomographic bins
"""

import numpy as np
import pandas as pd
import feather
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams["text.usetex"] =True


# +++++++++++++++ plot information
# input directory
inpathF = "/disks/shear15/ssli/KV450/tomo/all_tomo"
#

# half-half cut
halfline = [3.09, 3.18, 3.27, 3.27, 3.18]

# bin title
titles = [r'$0.1< z_{\rm B} \leq 0.3$', r'$0.3< z_{\rm B} \leq 0.5$', r'$0.5< z_{\rm B} \leq 0.7$', r'$0.7< z_{\rm B} \leq 0.9$', r'$0.9< z_{\rm B} \leq 1.2$']

# output directory
outfile_pdf = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/cumu_TB.pdf"


# plot related
LABELS = [r'$T_{\rm B} \leq 3$', r'$T_{\rm B} > 3$']
COLORS = ['red', 'blue']
XLABEL = r"$T_{\rm B}$"        
YLABEL = "Cumulative weights"
XLIM = [0.5, 5.5]
YLIM = [0, 1.0]

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

            inpath = inpathF + str(ibins) + '.feather'
            data = feather.read_dataframe(inpath)
            print("Data loaded from", inpath)

            data.sort_values(by=['T_B'], ascending=True, inplace=True)            
            TB = data['T_B'].values
            wg = data['recal_weight'].values
            # split 
            data1 = data[data['T_B']<=3]
            data2 = data[data['T_B']>3]
            TB1 = data1['T_B'].values
            wg1 = data1['recal_weight'].values
            TB2 = data2['T_B'].values
            wg2 = data2['recal_weight'].values

            # evaluate the histogram
            values1, base1 = np.histogram(TB1, bins=20, weights=wg1)
            values2, base2 = np.histogram(TB2, bins=20, weights=wg2)

            # evaluate the cumulative
            cumulative1 = np.cumsum(values1)/np.sum(wg)
            cumulative2 = np.cumsum(values2)/np.sum(wg) + cumulative1[-1]

            # for continuity
            cumulative1 = np.insert(cumulative1, 0, 0)
            base2[0] = base1[-1]
            cumulative2 = np.insert(cumulative2, 0, cumulative1[-1])
                    
            # plot the cumulative function
            l1 = axes[i, j].step(base1, cumulative1, c=COLORS[0], label=LABELS[0])
            l2 = axes[i, j].step(base2, cumulative2, c=COLORS[1], label=LABELS[1])
            if not handles:
                handles = [l1[0], l2[0]]

            axes[i, j].axvline(x=halfline[int(ibins-1)], ymin=0, ymax=0.5*(YLIM[1]-YLIM[0]), ls='--', lw=0.5, color='gray')
            axes[i, j].axhline(y=0.5, xmin=0, xmax=(halfline[int(ibins-1)]-XLIM[0])/(XLIM[1]-XLIM[0]), ls='--', lw=0.5, color='gray')

            # half value label
            axes[i, j].text(halfline[int(ibins-1)]+0.1, YLIM[1]*0.25, '{:0.2f}'.format(halfline[int(ibins-1)]), fontsize=8)


            # bin lable
            label = titles[int(ibins-1)]
            x = XLIM[0] + 0.1*(XLIM[1]-XLIM[0])
            y = YLIM[0] + 0.8*(YLIM[1]-YLIM[0])
            axes[i, j].text(x, y, label)

            # limits
            axes[i, j].set_xlim(XLIM[0], XLIM[1])
            axes[i, j].set_ylim(YLIM[0], YLIM[1])
            
            # ticks
            axes[i, j].set_yticks([0.1, 0.3, 0.5, 0.7, 0.9])
            axes[i, j].set_xticks([1, 3, 5])

            # lable
            # if j == 0:
                # axes[i, j].set_ylabel(YLABEL)
            if i == 1:
                axes[i, j].set_xlabel(XLABEL)
            if (i==0 and j==2):
                axes[i, j].set_xlabel(XLABEL)

        else:
            axes[i, j].axis('off')
        ibins += 1

# Y label
fig.text(0.04, 0.5, YLABEL, va='center', rotation='vertical')

fig.legend(handles, LABELS, loc = 'upper right', bbox_to_anchor=(0.88, 0.35), fancybox=True, shadow=True)


plt.savefig(outfile_pdf, dpi=300)
print("Plot saved in", outfile_pdf)



