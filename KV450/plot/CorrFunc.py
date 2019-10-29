#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 11:50:41 2019

@author: ssli

Plot 2-p correlation results
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# general settings for plot
# mpl.use('Agg')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rcParams['xtick.top'] = True
#mpl.rcParams['ytick.right'] = True
plt.rc('font', size=10)

# custom settings for plot
XLIM = [0.1, 300]
YLIM = [-3, 5]
# color
CR = ['black', 'blue', 'orange', 'red', 'cyan'] 
# marker
MK = [' ', 'o', 'o', 'o', 'o']
# marker size
MS = 4
# capthick
MW = 0.5
# linestyle
LS = ['-', '-', '--', '-', '--']
# linewidth
LW_H = 1.0
LW_G = 0.5

# plot xi_p or xi_m
# PORM = 'P'
PORM = 'M'

# input directory
inpath_HH = "/disks/shear15/ssli/KV450/CS/Hildebrandt2018/sheardata/KV450_COSMIC_SHEAR_DATA_RELEASE/DATA_VECTOR/KV450_xi_pm_tomographic_data_vector.dat"
# inpath_pre_slop = "/disks/shear15/ssli/KV450/CS/mine/TreeCorr/xi_pre_slop.csv"
inpath_pre = "/disks/shear15/ssli/KV450/CS/mine/TreeCorr/xi_pre.csv"
inpath_pre_withoutc = "/disks/shear15/ssli/KV450/CS/mine/TreeCorr/xi_pre_withoutc.csv"
inpath_new = "/disks/shear15/ssli/KV450/CS/mine/TreeCorr/xi_mine_new.csv"
inpath_old = "/disks/shear15/ssli/KV450/CS/mine/TreeCorr/xi_mine_old.csv"
# inpaths = [inpath_HH, inpath_pre_slop, inpath_pre, inpath_new, inpath_old]
inpaths = [inpath_HH, inpath_pre, inpath_pre_withoutc, inpath_new, inpath_old]
# names = ['HH', 'NoMagCut + MinusAveE + bin_slop', 'NoMagCut + MinusAveE', 'MagCut + MinusAveE', 'MagCut']
names = ['HH', '$pre$ with c-term', '$pre$ without c-term', '$new$ with c-term', '$new$ without c-term']

# output directory
outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/"

# data
paras = []
for path in inpaths:
    if path == inpath_HH:
        data = np.loadtxt(path)
        theta = data[:,1]
        xi = data[:,2]
        pm = data[:,3]
        ito = data[:,4]
        jto = data[:,5]
        para = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
    else:
        data = np.loadtxt(path, delimiter=',', skiprows=1)
        theta = data[:,1]
        xi = data[:,3]
        err = data[:,5]
        pm = data[:,6]
        ito = data[:,7]
        jto = data[:,8]
        para = pd.DataFrame({'theta': theta, 'xi': xi, 'err': err, 'pm': pm, 'ito': ito, 'jto': jto})
    paras.append(para)
    
# plot
fig, ax = plt.subplots(5, 5, sharey=True)
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
for i in range(5):
    for j in range(5):
        if j >= i:
            if PORM == 'P':
                l = 4 - j
                m = i
                PM = 1
            elif PORM == 'M':
                l = 4 - i
                m = j
                PM = 2
                XLIM = [5, 350]
            else:
                raise Exception("Unexpected PORM value.") 
                    
            for k in range(len(names)):
                para = paras[k]
                name = names[k]

                theta = para[(para.pm==PM) & (para.ito==i+1) & (para.jto==j+1)]['theta'].values
                xi = para[(para.pm==PM) & (para.ito==i+1) & (para.jto==j+1)]['xi'].values

                if name == 'HH':
                    if (i==0) and (j==0):
                        ax[l, m].errorbar(theta, theta*xi*1e4, color=CR[k], mfc='none', marker=MK[k], markersize=MS, mew=MW, linestyle=LS[k], linewidth=LW_H, label=name)
                    else:
                        ax[l, m].errorbar(theta, theta*xi*1e4, color=CR[k], mfc='none', marker=MK[k], markersize=MS, mew=MW, linestyle=LS[k], linewidth=LW_H)
                else:
                    err = para[(para.pm==PM) & (para.ito==i+1) & (para.jto==j+1)]['err'].values
                    if (i==0) and (j==0):
                        ax[l, m].errorbar(theta, theta*xi*1e4, yerr=theta*err*1e4, color=CR[k], mfc='none', marker=MK[k], markersize=MS, mew=MW, linestyle=LS[k], linewidth=LW_G, label=name)
                    else:
                        ax[l, m].errorbar(theta, theta*xi*1e4, yerr=theta*err*1e4, color=CR[k], mfc='none', marker=MK[k], markersize=MS, mew=MW, linestyle=LS[k], linewidth=LW_G)
                    
            # label
            label = str(i+1) + '-' + str(j+1)
            x = 10**(np.log10(XLIM[0]) + 0.1*(np.log10(XLIM[1])-np.log10(XLIM[0])))
            y = YLIM[0] + 0.7*(YLIM[1]-YLIM[0])
            ax[l,m].text(x, y, label)

            ax[l, m].set_xscale("log", nonposx='clip')
            ax[l, m].set_xlim(XLIM[0], XLIM[1])
            ax[l, m].set_ylim(YLIM[0], YLIM[1])
            ax[l, m].axhline(y=0., xmin=XLIM[0], xmax=XLIM[1], color = 'grey', linestyle = 'dotted')
            
            # tick
#            ax[l, m].tick_params(axis='both', which='major', width=1.00, length=5)
#            ax[l, m].tick_params(which='minor', width=0.75, length=2.5)
            if PORM == "P":
                ax[l, m].set_xticks([1, 10, 100])
            elif PORM == "M":
                ax[l, m].set_xticks([10, 100])

            ax[l, m].set_yticks([-2, 0, 2, 4])
            if (i==0) and (j==0):
                if PORM == 'P':
                    ax[l, m].set_xticklabels(['1', '10', '100'])
                else:
                    ax[l, m].set_xticklabels(['10', '100'])
                ax[l, m].set_yticklabels(['-2', '0', '2', '4'])
            else:
                ax[l, m].set_xticklabels(" ")
#                ax[l, m].set_yticklabels(" ")    
        else:
            if PORM == 'P':
                l = 4 - j
                m = i
                PM = 1
            elif PORM == 'M':
                l = 4 - i
                m = j
                PM = 2
            else:
                raise Exception("Unexpected PORM value.") 
            
            ax[l,m].axis('off')

if PORM == 'P':
    fig.legend(frameon=False, loc='lower right')
    fig.text(0.5, 0.05, r'$\theta$ [arcmin]', horizontalalignment='center', verticalalignment='center')
    fig.text(0.05, 0.5, r'$\theta \times \xi_+ [10^{-4}$ arcmin]', 
             horizontalalignment='center', verticalalignment='center',
             rotation='vertical')
    # plt.savefig(outpath+"compare_xip.png", dpi=300)
    plt.savefig(outpath+"compare_xip.pdf")
    plt.close()    
elif PORM == 'M':
    fig.legend(frameon=False, loc='upper left')
    fig.text(0.5, 0.05, r'$\theta$ [arcmin]', horizontalalignment='center', verticalalignment='center')
    fig.text(0.95, 0.5, r'$\theta \times \xi_- [10^{-4}$ arcmin]', 
             horizontalalignment='center', verticalalignment='center',
             rotation='vertical')
    # plt.savefig(outpath+"compare_xim.png", dpi=300)
    plt.savefig(outpath+"compare_xim.pdf")    
    plt.close()    
else:
    raise Exception("Unexpected PORM value.") 
             