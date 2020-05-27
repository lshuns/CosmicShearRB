# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-05-24 12:25:48
# @Last Modified by:   lshuns
# @Last Modified time: 2020-05-27 22:09:29

# plot of m-bias results

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

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

m_error = 0.02

# output directory
# outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/mbias/m_whole_rb.png"
outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/m_whole_rb.pdf"

# plot related
LABELS = ['whole', r'$T_{\rm B}\leq 3$', r'$T_{\rm B} > 3$']
COLORS = ['black', 'red', 'blue']
SYMBOLS = ['o', 'v', '^']
# marker size
MSs = [4, 4, 4]
# linewidth of the errorbar lines
ELWs = [1, 1, 1]


XTICK = [1, 2, 3, 4, 5]
XTICKLABELS = ['B1', 'B2', 'B3', 'B4', 'B5']
YLABEL = r"$m$"
YLIM = [-0.05, 0.05]
XLIM = [0.5, 5.5]


for i in range(len(infiles)):
    inpath = inDir + infiles[i]
    data = pd.read_csv(inpath)


    data.sort_values(by=['bin'], ascending=True, inplace=True)    
    binvalue = data['bin'].values + 0.05*(i+1)
    print('data sorted as bin', binvalue)

    mvalue = data['m'].values
    print('m-bias sorted as', mvalue)

    CR = COLORS[i]
    MK = SYMBOLS[i]
    MS = MSs[i]
    ELW = ELWs[i]
    LABEL = LABELS[i]
                   
    plt.errorbar(binvalue, mvalue, yerr=m_error,
                    color=CR, marker=MK, markersize=MS, elinewidth=ELW, ls='none', label=LABEL)


for ibin in range(5):
    plt.axvline(x=1.5+ibin, color='gray', ls='--', lw=1)
                        

plt.xticks(ticks=XTICK, labels=XTICKLABELS)
plt.tick_params(axis='x', length=0, width=0)

plt.legend()
plt.xlim(XLIM[0], XLIM[1])
plt.ylim(YLIM[0], YLIM[1])
plt.ylabel(YLABEL)

plt.savefig(outpath, dpi=300)
print('plot saved in', outpath)
plt.close()


