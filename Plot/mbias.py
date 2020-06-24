# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-05-24 12:25:48
# @Last Modified by:   lshuns
# @Last Modified time: 2020-06-11 14:55:11

# plot of m-bias results

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
plt.rc('font', size=12)


# +++ input
inDir = "/disks/shear15/ssli/CosmicShear/shear_bias/"
infiles = ["Summary_m_whole.csv", "Summary_m_less3.csv", "Summary_m_greater3.csv"]

m_error_budget = 0.02

# output directory
# outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/mbias/m_whole_rb.png"
outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/m_whole_rb.pdf"

# plot related
LABELS = ['whole', r'$T_{\rm B}\leq 3$', r'$T_{\rm B} > 3$']
COLORS = ['black', 'red', 'blue']
SYMBOLS = ['o', '<', '>']
# marker size
MSs = [4, 4, 4]
# linewidth of the errorbar lines
ELWs = [1, 1, 1]


YTICK = [1, 2, 3, 4, 5]
YTICKLABELS = [r'$0.1< z_{\rm B} \leq 0.3$', r'$0.3< z_{\rm B} \leq 0.5$', r'$0.5< z_{\rm B} \leq 0.7$', r'$0.7< z_{\rm B} \leq 0.9$', r'$0.9< z_{\rm B} \leq 1.2$']
XLABEL = r"$m$"
YLIM = [0.5, 5.5]
XLIM = [-0.05, 0.05]

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.2, bottom=0.1, right=0.98, top=0.9, wspace=0, hspace=0)


for i in range(len(infiles)):
    inpath = inDir + infiles[i]
    data = pd.read_csv(inpath)

    data.sort_values(by=['bin'], ascending=True, inplace=True)    

    # x position
    if i == 0:
        # whole in the middle
        binvalue = data['bin'].values
    elif i == 1: 
        # red in the left
        binvalue = data['bin'].values - 0.05
    elif i == 2: 
        # blue in the right
        binvalue = data['bin'].values + 0.05
    
    print('data sorted as bin', binvalue)

    mvalue = data['m'].values
    print('m-bias sorted as', mvalue)

    # error
    merror = data['m_err_BS'].values
    print('m error (BS)', merror)

    CR = COLORS[i]
    MK = SYMBOLS[i]
    MS = MSs[i]
    ELW = ELWs[i]
    LABEL = LABELS[i]
                   
    # error budget shadow
    if i == 0:
        for kk in range(len(binvalue)):
            x0 = mvalue[kk]
            y0 = binvalue[kk]

            x = np.array([x0-m_error_budget, x0+m_error_budget])
            y1 = y0 - 0.5
            y2 = y0 + 0.5
            plt.fill_between(x, y1, y2, edgecolor='gray', facecolor='white', alpha=0.3, hatch='/')

    plt.errorbar(mvalue, binvalue, xerr=merror,
                    color=CR, marker=MK, markersize=MS, elinewidth=ELW, ls='none', label=LABEL)


for ibin in range(5):
    plt.axhline(y=1.5+ibin, color='black', ls='-', lw=1)



plt.yticks(ticks=YTICK, labels=YTICKLABELS)
plt.tick_params(axis='y', length=0, width=0)

plt.xlim(XLIM[0], XLIM[1])
plt.ylim(YLIM[0], YLIM[1])
plt.xlabel(XLABEL)

# invert y-axis
plt.gca().invert_yaxis()

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.13),
             fancybox=True, shadow=True, ncol=3)

plt.savefig(outpath, dpi=300)
print('plot saved in', outpath)
plt.close()


