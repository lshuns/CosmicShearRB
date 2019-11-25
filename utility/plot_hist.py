#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 13:50:48 2019

@author: ssli

make the histogram plot
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

def HistPlotFunc(values, wgs, outpaths, 
    DENSITY, HISTTYPE, NBs, COLORs, LABLEs,
    XLABEL, YLABEL):
    """
    Function for histogram plot
    """
    print("Make the histogram plot (HistPlotFunc)...")

    # mpl.use('Agg')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True

    fig = plt.figure()
    for i in range(len(values)):
        value = values[i]
        wg = wgs[i]

        NB = NBs[i]
        COLOR = COLORs[i]
        LABLE = LABLEs[i]
    
        plt.hist(x=value, bins=NB, density=DENSITY, weights=wg, color=COLOR, label=LABLE, histtype=HISTTYPE)

    if LABELS[0] != None:
        plt.legend(frameon=False)

    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)

    for path in outpaths:
        plt.savefig(path, dpi=300)
        print("Hist plot saved to", path)
    plt.close()

    print("Finished the histogram plot (HistPlotFunc).")
