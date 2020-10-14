#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 15:01:11 2020

@author: ssli

Plot of data vector (theory or signal)
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, LinearLocator, LogLocator


def XiPlotFunc(paras, names, nzbins,
                CRs, MKs, MSs, LSs, LWs, ELWs,
                YTYPE, 
                outpath,
                LABELs=None):
    """
    xi_pm plot
    """

    # general settings for plot
    plt.rcParams["text.usetex"] =True
    mpl.use('Agg')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    #mpl.rcParams['xtick.top'] = True
    #mpl.rcParams['ytick.right'] = True
    plt.rc('font', size=9)


    fig, ax = plt.subplots(nzbins+2, nzbins+2, sharey=True)
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    # general parameters for plot
    XLIM_P = [0.4, 97.]
    XLIM_M = [3.5, 360.]
    YLIM = [-6, 6]
    # limits for hline
    XLIM_H = [0., 360.] 

    # gap between xip and xim
    for j in range(nzbins+2):
        # middle diagonal
        l = nzbins + 1 - j
        m = j
        ax[l, m].axis('off')
        # left diagonal
        if j < nzbins+1:
            l = nzbins - j
            m = j
            ax[l, m].axis('off')
        # right diagonal
        if j > 0:
            l = nzbins + 2 - j
            m = j
            ax[l, m].axis('off')

    # actual plot
    # for legend
    handles_data = []
    handles_theory = []
    LABELs_data = []
    LABELs_theory = []
    for i in range(nzbins):
        for j in range(nzbins):
            if j >= i:
                # xip
                l_P = nzbins - 1 - j
                m_P = i
                id_P = 1
                # xim
                l_M = nzbins + 1 - i
                m_M = j + 2
                id_M = 2            
                #
                ls = [l_P, l_M] 
                ms = [m_P, m_M] 
                ids = [id_P, id_M] 

                for index_pm in range(2):
                    l = ls[index_pm]
                    m = ms[index_pm]
                    id_pm = ids[index_pm]
                    if id_pm == id_P:
                        XLIM = XLIM_P
                    else:
                        XLIM = XLIM_M
                    # plot set
                    label = str(i+1) + '-' + str(j+1)
                    x = 10**(np.log10(XLIM[0]) + 0.1*(np.log10(XLIM[1])-np.log10(XLIM[0])))
                    y = YLIM[0] + 0.7*(YLIM[1]-YLIM[0])
                    ax[l,m].text(x, y, label)
                    #
                    ax[l, m].set_xscale("log", nonposx='clip')
                    ax[l, m].set_xlim(XLIM[0], XLIM[1])
                    ax[l, m].set_ylim(YLIM[0], YLIM[1])
                    # hline
                    ax[l, m].axhline(y=0., xmin=XLIM_H[0], xmax=XLIM_H[1], color = 'grey', linestyle = 'dotted', linewidth=0.5)

                    # data plot        
                    for k in range(len(names)):
                        para = paras[k]
                        name = names[k]
                        #
                        CR = CRs[k]
                        MK = MKs[k]
                        MS = MSs[k]
                        #
                        LS = LSs[k]
                        LW = LWs[k]
                        ELW = ELWs[k]
                        # Labels
                        if LABELs != None:
                            LAB = LABELs[k]

                        # data
                        theta = para[(para.pm==id_pm) & (para.ito==i+1) & (para.jto==j+1)]['theta'].values
                        xi = para[(para.pm==id_pm) & (para.ito==i+1) & (para.jto==j+1)]['xi'].values

                        if name == 'data':
                            err = para[(para.pm==id_pm) & (para.ito==i+1) & (para.jto==j+1)]['error'].values
                            if LABELs != None:
                                l_tmp = ax[l, m].errorbar(theta, theta*xi*1e4, yerr=theta*err*1e4,
                                    color=CR, marker=MK, markersize=MS, elinewidth=ELW, 
                                    linestyle=LS, linewidth=LW, label=LAB)
                                if i==0 and j==0 and index_pm==0:
                                    handles_data.append(l_tmp[0])
                                    LABELs_data.append(LAB)
                            else:
                                ax[l, m].errorbar(theta, theta*xi*1e4, yerr=theta*err*1e4,
                                    color=CR, marker=MK, markersize=MS, elinewidth=ELW, 
                                    linestyle=LS, linewidth=LW)
                        elif name == 'theory':
                            if LABELs != None:
                                l_tmp = ax[l, m].plot(theta, theta*xi*1e4, 
                                    color=CR, linestyle=LS, linewidth=LW, label=LAB)
                                if i==0 and j==0 and index_pm==0:
                                    handles_theory.append(l_tmp[0])
                                    LABELs_theory.append(LAB)
                            else:
                                ax[l, m].plot(theta, theta*xi*1e4, 
                                    color=CR, linestyle=LS, linewidth=LW)

                    # plot set
                    # ticks
                    ax[l, m].set_yticks([-4, 0, 4])
                    if id_pm == id_P:
                        ax[l, m].set_xticks([1, 10, 100])
                    elif id_pm == id_M:
                        ax[l, m].set_xticks([10, 100])

                    ax[l, m].xaxis.set_major_formatter(FormatStrFormatter('%.0f')) #ScalarFormatter('%.2e'))
                    # ax[l, m].tick_params(bottom=False, top=False, left=False, right=False, which='both')
                    #
                    # ax[l, m].yaxis.set_minor_locator(LinearLocator(presets=[-1.5, -1, -0.5, 0.5, 1.0, 1.5, 2.5, 3.0, 3.5, 4.5, 5.0, 5.5]))
                    # ax[l, m].yaxis.set_minor_locator(LinearLocator())


                    # xticklabels
                    if id_pm == id_P:
                        if (i==j):
                            ax[l, m].set_xticklabels(['1', '10', '100'])
                        else:
                            ax[l,m].set_xticklabels(" ")
                    elif id_pm == id_M:
                        if (i==0):
                            ax[l, m].set_xticklabels(['10', '100'])
                        else:
                            ax[l, m].set_xticklabels(" ")

                    # yticklabels
                    if id_pm == id_P:
                        ax[l, m].yaxis.tick_left()
                        # if (i==0):
                        #     ax[l, m].set_yticklabels(['0', '2', '4'])
                        # else:
                        #     ax[l, m].set_yticklabels(" ")
                    elif id_pm == id_M:
                        ax[l, m].yaxis.tick_right()
                        ax[l, m].yaxis.set_label_position("right")
                        # if (j==4):
                        #     ax[l, m].set_yticklabels(['0', '2', '4'])
                        # else:
                        #     ax[l, m].set_yticklabels(" ")

                    # xlabel
                    if (i==0) and (j==0):
                        ax[l, m].set_xlabel(r'$\theta~[{\rm arcmin}]$')

    # legend
    if LABELs != None:
        fig.legend(handles_theory, LABELs_theory, loc = 'upper left', bbox_to_anchor=(0.68, 0.76), frameon=False)
        fig.legend(handles_data, LABELs_data, loc = 'upper left', bbox_to_anchor=(0.68, 0.84), frameon=False)

    # ylabel
    if YTYPE == 'diff':
        fig.text(0.07, 0.5, r'$\theta \times \Delta\xi_+ [10^{-4}~{\rm arcmin}]$', 
             horizontalalignment='center', verticalalignment='center',
             rotation='vertical')
        fig.text(0.93, 0.5, r'$\theta \times \Delta\xi_- [10^{-4}~{\rm arcmin}]$', 
                 horizontalalignment='center', verticalalignment='center',
                 rotation='vertical')
    elif YTYPE == 'orig':
        fig.text(0.07, 0.5, r'$\theta \times \xi_+ [10^{-4}~{\rm arcmin}]$', 
             horizontalalignment='center', verticalalignment='center',
             rotation='vertical')
        fig.text(0.93, 0.5, r'$\theta \times \xi_- [10^{-4}~{\rm arcmin}]$', 
             horizontalalignment='center', verticalalignment='center',
             rotation='vertical')
    else:
        raise Exception('Unexpected YTYPE!')
    
    plt.savefig(outpath, dpi=300)
    plt.close()    
