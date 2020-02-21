#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 10:45:46 2019

@author: ssli

histogram comparing for desired parameters
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

def HistPlotFunc(paras, wgs, outpath,
                XLABEL, YLABEL,
                COLORS, NBINS, LABELS=[], 
                DENSITY=True, HISTTYPE='step',
                TITLE=None,
                vline=None, vline_style=None):
    """
    Histogram plot for multiple parameters
    """

    mpl.use('Agg')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True

    fig = plt.figure()
    for i in range(len(paras)):
        para = paras[i]
        wg = wgs[i]
        CR = COLORS[i]
        NB = NBINS[i]
        if LABELS != []:
            LAB = LABELS[i]
            plt.hist(x=para, bins=NB, density=DENSITY, weights=wg, color=CR, label=LAB, histtype=HISTTYPE)
        else:
            plt.hist(x=para, bins=NB, density=DENSITY, weights=wg, color=CR, histtype=HISTTYPE)
        

    if LABELS != []:
        plt.legend(frameon=False)

    if vline != None:
        plt.axvline(x=vline, ls=vline_style)
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    if TITLE != None:
        plt.title(TITLE)

    plt.savefig(outpath, dpi=300)
    plt.close()
    print("Histogram plot saved in", outpath)
