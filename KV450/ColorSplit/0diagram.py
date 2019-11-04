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


def CMdiagramFunc(patches, inpaths, outpaths, bandX, bandY, TB, 
                    MARKSIZE, XLIM, YLIM, XLABEL, YLABEL):
    """
    Function for color-magnitude diagram
    """

    # general settings for plot
    # mpl.use('Agg')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True
    plt.rc('font', size=10)


    # data
    data = []
    for i in range(len(patches)):
        patch = patches[i]
        inpath = inpaths[i]

        tmp = pd.read_feather(inpath)

        data.append(tmp)
        print("Data loaded from", inpath)

    data = pd.concat(data)

    magx = data[bandX].values
    magy = data[bandY].values
    vTB = data[TB].values

    # cm = plt.cm.get_cmap('RdYlBu')
    # sc = plt.scatter(magx, magy-magx, c=vTB, s=MARKSIZE, vmin=0, vmax=np.amax(vTB), cmap=cm)
    # plt.xlim(XLIM[0], XLIM[1])
    # plt.ylim(YLIM[0], YLIM[1])
    # plt.xlabel(XLABEL)
    # plt.ylabel(YLABEL)
    # cbar = plt.colorbar(sc)
    # cbar.ax.set_ylabel(r'$T_B$')


    magx_r = magx[vTB<=1.9]
    magy_r = magy[vTB<=1.9]
    vTB_r = vTB[vTB<=1.9]

    magx_b = magx[vTB>1.9]
    magy_b = magy[vTB>1.9]
    vTB_b = vTB[vTB>1.9]

    sc_r = plt.scatter(magx_r, magy_r-magx_r, c='red', s=MARKSIZE, label=r'$T_B\leq 1.9$')
    sc_b = plt.scatter(magx_b, magy_b-magx_b, facecolors='none', edgecolors='blue', linewidth=0.1, s=MARKSIZE, label=r'$T_B > 1.9$')
    plt.xlim(XLIM[0], XLIM[1])
    plt.ylim(YLIM[0], YLIM[1])
    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.legend()

    
    for outpath in outpaths:   
        plt.savefig(outpath, dpi=300)
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

    patches = ["G9","G12","G15","G23","GS"]

    # input path
    pathF = "/disks/shear15/ssli/KV450/CorrFunc/data/feather/"
    pathP = "_allBins.feather"
    inpaths = []
    for patch in patches:
        inpath = pathF + patch + pathP
        inpaths.append(inpath)

    # output path 
    outpath1 = "/disks/shear15/ssli/KV450/ColorSplit/CMdiagram_RB.png"
    outpath2 = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/ColorSplit/CMdiagram_RB.png"
    outpaths = [outpath1, outpath2]

    # plot related
    MARKSIZE = 0.5
    XLIM = [18, 27]
    YLIM = [-6, 6]
    XLABEL = '$m_r$'
    YLABEL = '$m_g-m_r$'


    CMdiagramFunc(patches, inpaths, outpaths, bandX, bandY, TB, 
                    MARKSIZE, XLIM, YLIM, XLABEL, YLABEL)


