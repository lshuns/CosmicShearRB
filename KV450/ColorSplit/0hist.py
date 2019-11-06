#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 12:11:21 2019

@author: ssli

Plot the hist of desired parameter

"""

import numpy as np
import pandas as pd
import feather
import matplotlib as mpl
import matplotlib.pyplot as plt

def histFunc(values, wgs, outdirs,
                COLORS, NBINS, LABELS,
                XLABEL, YLABEL):
    """
    Function for histogram plot
    """
    # mpl.use('Agg')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True

    
    fig = plt.figure()
    for i in range(len(values)):
        value = values[i]
        wg = wgs[i]
        CR = COLORS[i]
        NB = NBINS[i]
        LAB = LABELS[i]

        plt.hist(x=value, bins=NB, density=True, weights=wg, color=CR, label=LAB, histtype='step')

    if LABELS[0] != None:
        plt.legend(frameon=False)
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    for outdir in outdirs:
        plt.savefig(outdir, dpi=300)
        print("Plot saved in", outdir)

    plt.close()

if __name__=="__main__":
    
    bins = ["13", "35", "57", "79", "912"]
    patches = ["G9","G12","G15","G23","GS"]

    # input directory
    inpathF = "/disks/shear15/ssli/KV450/CorrFunc/data/feather/"
    inpathP = ".feather"

    para_n = "T_B"
    wg_n = "recal_weight"

    # output directory
    outpathF1 = "/disks/shear15/ssli/KV450/ColorSplit/hist_TB_"
    outpathF2 = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/ColorSplit/hist_TB_"
    outpathP = ".png"

    # plot related
    XLABEL = "$T_B$"
    YLABEL = "Distribution"

    for Bin in bins:
        data = []
        for patch in patches:
            inpath = inpathF + patch + "_" + Bin + inpathP
            tmp = feather.read_dataframe(inpath)
            tmp = tmp[["T_B", "recal_weight"]]
            data.append(tmp)
            print("Data loaded from", inpath)
        data = pd.concat(data)
        print("Data combined.")

        values = [data["T_B"]]
        wgs = [data["recal_weight"]]

        outpath1 = outpathF1 + Bin + outpathP
        outpath2 = outpathF2 + Bin + outpathP
        outdirs = [outpath1, outpath2]
        
        # plot related
        COLORS = ['blue']
        NBINS = [60]
        LABELS = [None]

        histFunc(values, wgs, outdirs,
                COLORS, NBINS, LABELS,
                XLABEL, YLABEL)