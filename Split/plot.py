#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:50:06 2019

@author: ssli

HistFunc:
Plot the hist of desired parameter for multiple samples in a single plot

CMDFunc:
Plot the color-magnitude diagram for multiple samples in a single plot
"""

import numpy as np
import pandas as pd
import feather
import matplotlib as mpl
import matplotlib.pyplot as plt


def HistFunc(values, wgs, outdirs,
                COLORS, NBINS, LABELS,
                XLABEL, YLABEL):
    """
    Function for histogram plot
    """
    print("Start histogram plot (histFunc)...")
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
        print("Histogram plot saved in", outdir)

    plt.close()
    print("Finished histogram plot (histFunc).")



def CMDFunc(samples, outpaths, bandX, bandY,
                    COLORS, LABELS, MARKSIZE, XLIM, YLIM, XLABEL, YLABEL):
    """
    Function for color-magnitude diagram
    """
    print("Start CMD plot (CMDFunc)...")

    # general settings for plot
    # mpl.use('Agg')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    plt.rc('font', size=10)


    cm = plt.cm.get_cmap('RdYlBu')
    for i in range(len(samples)):
        sample = samples[i]
        CR = COLORS[i]
        LA = LABELS[i]

        magx = sample[bandX].values
        magy = sample[bandY].values
        
        sc = plt.scatter(magx, magy-magx, c=CR, s=MARKSIZE, label=LA)
    
    plt.xlim(XLIM[0], XLIM[1])
    plt.ylim(YLIM[0], YLIM[1])
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.legend()
    
    for outpath in outpaths:   
        plt.savefig(outpath, dpi=300)
        print("CMD plot saved to", outpath)
    plt.close()

    print("Finished CMD plot (CMDFunc).")


if __name__=="__main__":

    # ++++++++++++++++++++++++++++++++++++++++++++++ HistFunc
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

        HistFunc(values, wgs, outdirs,
                COLORS, NBINS, LABELS,
                XLABEL, YLABEL)


    # ++++++++++++++++++++++++++++++++++++++++++++++ CMDFunc
    # All bands
    # ['MAG_GAAP_r_CALIB', 'MAG_GAAP_u_CALIB', 'MAG_GAAP_g_CALIB', 'MAG_GAAP_i_CALIB', 
    #     'MAG_GAAP_Z', 'MAG_GAAP_Y', 'MAG_GAAP_J', 'MAG_GAAP_H', 'MAG_GAAP_Ks']
    # band with good quality: g, r, i, Z

    # follow Johnston et al. 2019 # also backed by 'outlier' counts (99 or -99)
    bandX = 'MAG_GAAP_r_CALIB' # always use r-band for magnitude
    bandY = 'MAG_GAAP_g_CALIB' 
    TB = 'T_B'

    patches = ["G9","G12","G15","G23","GS"]
    bins = ["13", "35", "57", "79", "912"]

    # input directory
    inpathF = "/disks/shear15/ssli/KV450/ColorSplit/data/halfHalf/"
    inpathP = ".feather"

    # output directory
    outpathF1 = "/disks/shear15/ssli/KV450/ColorSplit/CMD/"
    outpathF2 = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/ColorSplit/CMD/"
    outpathP = ".png"

    # plot related
    MARKSIZE = 0.5
    XLIM = [18, 27]
    YLIM = [-6, 6]
    XLABEL = '$m_r$'
    YLABEL = '$m_g-m_r$'
    COLORS = ['red', 'blue']
    LABELS = ['small $T_B$', 'big $T_B$']


    for Bin in bins:
        for patch in patches:
            inpath1 = inpathF + patch + "_" + Bin + "_head" + inpathP
            sample1 = feather.read_dataframe(inpath1)
            print("Loaded data from", inpath1)

            inpath2 = inpathF + patch + "_" + Bin + "_tail" + inpathP
            sample2 = feather.read_dataframe(inpath2)
            print("Loaded data from", inpath2)

            samples = [sample1, sample2]

            # output path 
            outpath1 = outpathF1 + patch + "_" + Bin + outpathP
            outpath2 = outpathF2 + patch + "_" + Bin + outpathP
            outpaths = [outpath1, outpath2]

            CMDFunc(samples, outpaths, bandX, bandY,
                            COLORS, LABELS, MARKSIZE, XLIM, YLIM, XLABEL, YLABEL)


