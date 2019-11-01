#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 13:41:27 2019

@author: ssli

Plot the color-magnitude diagram
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


def CMdiagramFunc(patches, inpaths, outpath, Z, bins, bandX, bandY, TB, 
                    XLABEL, YLABEL):
    """
    Function for color-magnitude diagram
    """

    # general settings for plot
    # mpl.use('Agg')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    #mpl.rcParams['xtick.top'] = True
    #mpl.rcParams['ytick.right'] = True
    plt.rc('font', size=10)


    # data
    data = []
    for i in range(len(patches)):
        patch = patches[i]
        inpath = inpaths[i]

        hdf = pd.HDFStore(inpath, mode='r')
        print("hdf built from", inpath)

        for i in range(len(bins)-1):
        
            key = Z + '__' + str(bins[i]) + str(bins[i+1])

            tmp = hdf.select(key=key, columns=[bandX, bandY, TB])
            print("Selected data from", key)
            data.append(tmp)
        hdf.close()
        print("Done with data selection for", patch)
    data = pd.concat(data)

    magx = data[bandX].values
    magy = data[bandY].values
    vTB = data[TB].values

    cm = plt.cm.get_cmap('RdYlBu')
    sc = plt.scatter(magx, magy-magx, c=vTB, vmin=0, vmax=np.amax(TB), cmap=cm)
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.colorbar(sc)
    for outpath in outpaths:   
        plt.savefig(outpath)
    plt.close()


if __name__ == "__main__":

    # All bands
    # ['MAG_GAAP_r_CALIB', 'MAG_GAAP_u_CALIB', 'MAG_GAAP_g_CALIB', 'MAG_GAAP_i_CALIB', 
    #     'MAG_GAAP_Z', 'MAG_GAAP_Y', 'MAG_GAAP_J', 'MAG_GAAP_H', 'MAG_GAAP_Ks']
    # band with good quality: g, r, i, Z

    # follow Johnston et al. 2019 # also backed by 'outlier' counts (99 or -99)
    bandX = 'MAG_GAAP_r_CALIB' # always use r-band for magnitude
    bandY = 'MAG_GAAP_g_CALIB' 
    TB = 'T_B'

    bins = [1, 3, 5, 7, 9, 12]
    patches = ["G9","G12","G15","G23","GS"]
    key = "_whole"
    # column used as redshift
    Z = 'Z_B'

    # input path
    pathF = "/disks/shear15/ssli/KV450/CorrFunc/data/"
    pathP = ".h5"
    inpaths = []
    for patch in patches:
        inpath = pathF + patch + key + pathP
        inpaths.append(inpath)

    # output path 
    outpath1 = "/disks/shear15/ssli/KV450/ColorSplit/CMdiagram.png"
    outpath2 = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/ColorSplit/CMdiagram.png"
    outpaths = [outpath1, outpath2]

    # plot related
    XLABEL = '$m_r$'
    YLABEL = '$m_g-m_r$'

    CMdiagramFunc(patches, inpaths, outpaths, Z, bins, bandX, bandY, TB, 
                    XLABEL, YLABEL)


