#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:50:06 2019

@author: ssli

HistPlotFunc:
Plot the hist of desired parameter for multiple samples in a single plot

CMDPlotFunc:
Plot the color-magnitude diagram for multiple samples in a single plot
"""

import numpy as np
import pandas as pd
import feather
import matplotlib as mpl
import matplotlib.pyplot as plt


def HistPlotFunc(values, wgs, outdirs,
                vline, vline_style,
                COLORS, NBINS, LABELS,
                XLABEL, YLABEL, TITLE):
    """
    Function for histogram plot
    """
    print("Start histogram plot (HistPlotFunc)...")
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

    if vline != None:
        plt.axvline(x=vline, ls=vline_style)
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    if TITLE != None:
        plt.title(TITLE)

    for outdir in outdirs:
        plt.savefig(outdir, dpi=300)
        print("Histogram plot saved in", outdir)

    plt.close()
    print("Finished histogram plot (HistPlotFunc).")



def CMDPlotFunc(magxs, magys, outpaths,
                    COLORS, LABELS, MARKS, MARKSIZE, XLIM, YLIM, XLABEL, YLABEL, TITLE):
    """
    Function for color-magnitude diagram
    """
    print("Start CMD plot (CMDPlotFunc)...")

    # general settings for plot
    # mpl.use('Agg')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    plt.rc('font', size=10)


    cm = plt.cm.get_cmap('RdYlBu')
    for i in range(len(magxs)):

        MK = MARKS[i]
        LA = LABELS[i]
        CR = COLORS[i]

        magx = magxs[i]
        magy = magys[i]
        
        # sc = plt.scatter(magx, magy-magx, c=CR, s=MARKSIZE, label=LA)
        plt.plot(magx, magy-magx, MK, fillstyle='none', color=CR, markersize=MARKSIZE, label=LA)
    
    plt.xlim(XLIM[0], XLIM[1])
    plt.ylim(YLIM[0], YLIM[1])
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.legend()
    if TITLE != None:
        plt.title(TITLE)

    
    for outpath in outpaths:   
        plt.savefig(outpath, dpi=300)
        print("CMD plot saved to", outpath)
    plt.close()

    print("Finished CMD plot (CMDPlotFunc).")


if __name__=="__main__":
    import time
    Start = time.time()

    # ++++++++++++++++++++++++++++++++++++++++++++++ HistPlotFunc & CMDPlotFunc
    bins = ["13", "35", "57", "79", "912"]
    TITLES = [r"$0.1 < z_{\rm B} \leq 0.3$", r"$0.3 < z_{\rm B} \leq 0.5$", r"$0.5 < z_{\rm B} \leq 0.7$", r"$0.7 < z_{\rm B} \leq 0.9$", r"$0.9 < z_{\rm B} \leq 1.2$"]
    patches = ["G9","G12","G15","G23","GS"]

    # input directory
    inpathF = "/disks/shear15/ssli/KV450/Split/data/halfHalf/"
    inpathP = ".feather"

    # All bands
    # ['MAG_GAAP_r_CALIB', 'MAG_GAAP_u_CALIB', 'MAG_GAAP_g_CALIB', 'MAG_GAAP_i_CALIB', 
    #     'MAG_GAAP_Z', 'MAG_GAAP_Y', 'MAG_GAAP_J', 'MAG_GAAP_H', 'MAG_GAAP_Ks']
    # band with good quality: g, r, i, Z

    # follow Johnston et al. 2019 # also backed by 'outlier' counts (99 or -99)
    bandX = 'MAG_GAAP_r_CALIB' # always use r-band for magnitude
    bandY = 'MAG_GAAP_g_CALIB' 
    TB_n = 'T_B'
    wg_n = "recal_weight"


    # output directory
    outpathF_hist = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/Split/hist/TB_"
    outpathF_CMD = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/Split/CMD/TB_"
    outpathP = ".png"


    # plot related
    XLABEL_hist = "$T_B$"
    YLABEL_hist = "Distribution"
    LABELS_hist = [None, None]
    NBINS =[60, 60]
    vline = 3
    vline_style = '--'

    # plot related
    MARKSIZE = 2
    XLIM = [18, 27]
    YLIM = [-6, 6]
    XLABEL_CMD = '$m_r$'
    YLABEL_CMD = '$m_g-m_r$'
    COLORS = ['red', 'blue']
    LABELS_CMD = ['small $T_B$', 'big $T_B$']
    MARKS = ['x', 'o']


    for k in range(len(bins)):
        Bin = bins[k]
        data1 = []
        data2 = []
        for patch in patches:
            inpath1 = inpathF + patch + "_" + Bin + "_head_wg" + inpathP
            sample1 = feather.read_dataframe(inpath1)
            sample1 = sample1[[bandX, bandY, wg_n, TB_n]]
            data1.append(sample1)
            print("Loaded data from", inpath1)

            inpath2 = inpathF + patch + "_" + Bin + "_tail_wg" + inpathP
            sample2 = feather.read_dataframe(inpath2)
            sample2 = sample2[[bandX, bandY, wg_n, TB_n]]
            data2.append(sample2)
            print("Loaded data from", inpath2)

        data1 = pd.concat(data1)
        data2 = pd.concat(data2)
        print("Data combined.")

        TITLE = TITLES[k]

        # CMD
        magxs = [data1[bandX], data2[bandX]]
        magys = [data1[bandY], data2[bandY]]
        # output path 
        outpath = outpathF_CMD + Bin + outpathP
        outpaths = [outpath]
        CMDPlotFunc(magxs, magys, outpaths,
                        COLORS, LABELS_CMD, MARKS, MARKSIZE, XLIM, YLIM, XLABEL_CMD, YLABEL_CMD, TITLE)

        # hist
        values = [data1[TB_n], data2[TB_n]]
        wgs = [data1[wg_n], data2[wg_n]]
        #
        outpath = outpathF_hist + Bin + outpathP
        outdirs = [outpath]
        HistPlotFunc(values, wgs, outdirs,
                vline, vline_style,
                COLORS, NBINS, LABELS_hist,
                XLABEL_hist, YLABEL_hist, TITLE)

    print("Finished in", time.time()-Start)
    # 414.486989736557