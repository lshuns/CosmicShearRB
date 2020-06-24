#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:12:58 2020

@author: ssli

plot of redshift distribution from different samples.
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
inpathF_less3 = "/disks/shear15/ssli/CosmicShear/redshift/less3/Nz_DIR/Nz_DIR_Mean/Nz_DIR_"
inpathF_greater3 = "/disks/shear15/ssli/CosmicShear/redshift/greater3/Nz_DIR/Nz_DIR_Mean/Nz_DIR_"
inpathFs = [inpathF_less3, inpathF_greater3]
inpathPs = ['z0.1t0.3.asc', 'z0.3t0.5.asc', 'z0.5t0.7.asc', 'z0.7t0.9.asc', 'z0.9t1.2.asc']

# redshift bins
z_mins = [0.1, 0.3, 0.5, 0.7, 0.9]
z_maxs = [0.3, 0.5, 0.7, 0.9, 1.2]

# redshift properties
z_mean_less3 = [0.351, 0.430, 0.546, 0.744, 0.909]
z_mean_greater3 = [0.437, 0.573, 0.791, 0.914, 1.081]

z_median_less3 = [0.282, 0.396, 0.531, 0.732, 0.894]
z_median_greater3 = [0.244, 0.431, 0.644, 0.842, 1.022]

dz_mean = []
dz_median = []
for i in range(len(z_mean_greater3)):
    dz_mean.append('{:+0.2f}'.format(z_mean_greater3[i]-z_mean_less3[i]))
    dz_median.append('{:+0.2f}'.format(z_median_greater3[i]-z_median_less3[i]))

# output directory
outfile_pdf = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/distri_z.pdf"
outfile_png = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/distri_z.png"

# plot related
LABELS = [r'$T_{\rm B} \leq 3$', r'$T_{\rm B} > 3$']
COLORS = ['red', 'blue']
XLABEL = r'$z$'        
YLABEL = r'$n(z)$'
XLIM = [0, 2.5]
YLIM = [0, 6]

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
            print("Plot for bin", ibins)
            for k in range(len(inpathFs)):
                inpath = inpathFs[k] + inpathPs[int(ibins-1)]
                data = np.loadtxt(inpath)
                x_val = data[:,0]
                y_val = data[:,1]
                print("Data loaded from", inpath)

                # plot
                CR = COLORS[k]
                LAB = LABELS[k]

                p_tmp = axes[i, j].plot(x_val, y_val, color=CR, label=LAB)
                
                # for legend
                if ibins == 1:
                    handles.append(p_tmp[0])

            # shadow
            axes[i, j].axvspan(z_mins[int(ibins-1)], z_maxs[int(ibins-1)], alpha=0.3, color='grey')
            #
            axes[i, j].set_xlim(XLIM[0], XLIM[1])
            axes[i, j].set_ylim(YLIM[0], YLIM[1])
            axes[i, j].set_xticks([0.5, 1.0, 1.5, 2.0])
            axes[i, j].set_yticks([0, 2, 4])

            # bin label
            label = titles[int(ibins-1)]
            x = XLIM[0] + 0.3*(XLIM[1]-XLIM[0])
            y = YLIM[0] + 0.8*(YLIM[1]-YLIM[0])
            x2 = XLIM[0] + 0.33*(XLIM[1]-XLIM[0])
            y2 = YLIM[0] + 0.65*(YLIM[1]-YLIM[0])

            axes[i,j].text(x, y, label)
            axes[i,j].text(x2, y2, r'$[{:}, {:}]$'.format(dz_mean[int(ibins-1)], dz_median[int(ibins-1)]))

            if j == 0:
                axes[i, j].set_ylabel(YLABEL)
            if i == 1:
                axes[i, j].set_xlabel(XLABEL)
            if (i==0 and j==2):
                axes[i, j].set_xlabel(XLABEL)

        else:
            axes[i, j].axis('off')
        ibins += 1

fig.legend(handles, LABELS, loc = 'upper right', bbox_to_anchor=(0.85, 0.4), frameon=False)
fig.text(0.67, 0.23, r'$[\delta z_{\rm mean}, \delta z_{\rm median}]$')

plt.savefig(outfile_pdf)
print("Plot saved in", outfile_pdf)
plt.savefig(outfile_png, dpi=300)
print("Plot saved in", outfile_png)



