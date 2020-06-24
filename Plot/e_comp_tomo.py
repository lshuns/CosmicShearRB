#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:12:58 2020

@author: ssli

make histogram plot of ellipticity from different samples.
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
# inpathF_whole = "/disks/shear15/ssli/KV450/tomo/all_tomo"
inpathF_split = "/disks/shear15/ssli/KV450/split/all_tomo"
# inpathP_whole = ".feather"
inpathP_less3 = "_T_B_less3.feather"
inpathP_greater3 = "_T_B_greater3.feather"
#
# inpathFs = [inpathF_whole, inpathF_split, inpathF_split]
# inpathPs = [inpathP_whole, inpathP_less3, inpathP_greater3]
inpathFs = [inpathF_split, inpathF_split]
inpathPs = [inpathP_less3, inpathP_greater3]

# output directory
outfile_pdf = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/hist_e.pdf"
outfile_png = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/hist_e.png"

# plot related
LABELS = [r'$T_{\rm B} \leq 3$', r'$T_{\rm B} > 3$']
COLORS = ['red', 'blue']
NBINS = [60, 60]
XLABEL = r'$\epsilon$'        
YLABEL = "Weighted density"
XLIM = [0, 1]
YLIM = [0, 2.4]
DENSITY = True
HISTTYPE = 'step'

# bin title
titles = [r'$0.1< z_{\rm B} \leq 0.3$', r'$0.3< z_{\rm B} \leq 0.5$', r'$0.5< z_{\rm B} \leq 0.7$', r'$0.7< z_{\rm B} \leq 0.9$', r'$0.9< z_{\rm B} \leq 1.2$']


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
            for k in range(len(inpathFs)):
                inpath = inpathFs[k] + str(ibins) + inpathPs[k]
                data = feather.read_dataframe(inpath)
                e = np.sqrt((data['bias_corrected_e1'].values)**2. + (data['bias_corrected_e2'].values)**2.)
                wg = data['recal_weight']
                print("Data loaded from", inpath)

                # plot
                CR = COLORS[k]
                NB = NBINS[k]
                LAB = LABELS[k]

                _, _, p_tmp = axes[i, j].hist(x=e, bins=NB, density=DENSITY, weights=wg, color=CR, label=LAB, histtype=HISTTYPE)
                
                # for legend
                if ibins == 1:
                    handles.append(p_tmp[0])

            axes[i, j].set_xlim(XLIM[0], XLIM[1])
            axes[i, j].set_ylim(YLIM[0], YLIM[1])
            axes[i, j].set_xticks([0.2, 0.4, 0.6, 0.8])

            # bin label
            label = titles[int(ibins-1)]
            x = XLIM[0] + 0.3*(XLIM[1]-XLIM[0])
            y = YLIM[0] + 0.8*(YLIM[1]-YLIM[0])
            axes[i,j].text(x, y, label)

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


plt.savefig(outfile_pdf)
print("Plot saved in", outfile_pdf)
plt.savefig(outfile_png, dpi=300)
print("Plot saved in", outfile_png)



