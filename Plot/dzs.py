# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-06-04 11:32:20
# @Last Modified by:   lshuns
# @Last Modified time: 2020-07-12 13:45:43


# plot of dz,s results

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams["text.usetex"] =True


# +++ general settings for plot
mpl.use('Agg')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rcParams['xtick.top'] = True
#mpl.rcParams['ytick.right'] = True
plt.rc('font', size=16)


# +++ input

# KV450
D_z1_s = [0.0316, -0.0849, +0.1424]
D_z2_s = [0.0801, -0.0679, +0.0868]
D_z3_s = [0.0663, -0.0509, +0.0566] 
D_z4_s = [0.0024, -0.0498, +0.0482]
D_z5_s = [-0.0019, -0.0507, +0.0533]
D_z_s_K = [D_z1_s[0], D_z2_s[0], D_z3_s[0], D_z4_s[0], D_z5_s[0]]
D_z_s_K_error = [[-D_z1_s[1], -D_z2_s[1], -D_z3_s[1], -D_z4_s[1], -D_z5_s[1]], [D_z1_s[2], D_z2_s[2], D_z3_s[2], D_z4_s[2], D_z5_s[2]]]

# Planck
D_z1_s = [0.0036, -0.0981, +0.1222]
D_z2_s = [0.0389, -0.0529, +0.0591]
D_z3_s = [0.0403, -0.0393, +0.0421]
D_z4_s = [0.0087, -0.0389, +0.0394]
D_z5_s = [0.0076, -0.0464, +0.0477]
D_z_s_P = [D_z1_s[0], D_z2_s[0], D_z3_s[0], D_z4_s[0], D_z5_s[0]]
D_z_s_P_error = [[-D_z1_s[1], -D_z2_s[1], -D_z3_s[1], -D_z4_s[1], -D_z5_s[1]], [D_z1_s[2], D_z2_s[2], D_z3_s[2], D_z4_s[2], D_z5_s[2]]]

# mean difference
mean_diff = np.array([0.09,0.14,0.24,0.17,0.17])/2.

# Bin value
binvalue = np.array([1, 2, 3, 4, 5])

# output directory
# outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/Dzs.pdf"
outpath = "/Users/lshuns/surfdrive/Projects/6CosmicShear_RB/plot/publish/Dzs.pdf"


# plot related
LABELS = ['KV450', 'Planck']
COLORS = ['red', 'black']
SYMBOLS = ['o', 's']
# marker size
MSs = [4, 4]
# linewidth of the errorbar lines
ELWs = [1, 1]

YTICK = [1, 2, 3, 4, 5]
YTICKLABELS = [r'$0.1< z_{\rm B} \leq 0.3$', r'$0.3< z_{\rm B} \leq 0.5$', r'$0.5< z_{\rm B} \leq 0.7$', r'$0.7< z_{\rm B} \leq 0.9$', r'$0.9< z_{\rm B} \leq 1.2$']
XLABEL = r"$\delta_{z_i,{\rm s}}$"
YLIM = [0.5, 5.5]
XLIM = [-0.3, 0.3]

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.1, right=0.98, top=0.88, wspace=0, hspace=0)


# KV450                   
plt.errorbar(D_z_s_K, binvalue+0.05, xerr=D_z_s_K_error,
                color=COLORS[0], marker=SYMBOLS[0], markersize=MSs[0], elinewidth=ELWs[0], ls='none', label=LABELS[0])
# Planck                  
plt.errorbar(D_z_s_P, binvalue-0.05, xerr=D_z_s_P_error,
                color=COLORS[1], marker=SYMBOLS[1], markersize=MSs[1], elinewidth=ELWs[1], ls='none', label=LABELS[1])


ysc = [0,0.2,0.4,0.6,0.8,1.0]
for ibin in range(5):
    vx = mean_diff[ibin]
    print("mean diff in {:}: {:}".format(ibin,vx))
    plt.axvline(x=vx, ymin=ysc[ibin], ymax=ysc[ibin+1], color='blue', ls='--', lw=1)

    plt.axhline(y=1.5+ibin, color='black', ls='-', lw=1)
                        
plt.axvline(x=0, color='gray', ls='dotted', lw=1)
    
plt.yticks(ticks=YTICK, labels=YTICKLABELS)
plt.tick_params(axis='y', length=0, width=0)

plt.xlim(XLIM[0], XLIM[1])
plt.ylim(YLIM[0], YLIM[1])
plt.xlabel(XLABEL)

# invert y-axis
plt.gca().invert_yaxis()

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
             fancybox=True, shadow=True, ncol=2)

plt.savefig(outpath, dpi=300)
print('plot saved in', outpath)
plt.close()


