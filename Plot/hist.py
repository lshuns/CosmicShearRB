#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 10:45:46 2019

@author: ssli

histogram comparing simulated images with KV450 
"""

import numpy as np
import pandas as pd
import feather
import matplotlib as mpl
import matplotlib.pyplot as plt

import os
import sys
sys.path.insert(0,os.path.realpath('..')) 
from utility import HistPlotFunc


bins = ["13", "35", "57", "79", "912"]
patches = ["G9","G12","G15","G23","GS"]

# input directory
inpathF_KV450 = "/disks/shear15/ssli/KV450/CorrFunc/data/feather/"
inpathF_Sim = "/disks/shear15/ssli/SimCat/"
inpathP = ".feather"

# parameters to be plotted
para_KV450 = "T_B"
wg_KV450 = "recal_weight"
#
para_Sim = "TB9_in"
wg_Sim = "LFweight"

# output directory
outpathF = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/compare_KV450_Sim/hist_TB_"
outpathP = ".png"

# plot related
XLABEL = "$T_B$"
YLABEL = "Distribution"
COLORs = ['black', 'red']
NBs = [60, 60]
LABELs = ['KV450', 'Simulated data']
DENSITY = True
HISTTYPE = 'step'

for Bin in bins:
    data = []
    # KV450 data
    for patch in patches:
        inpath = inpathF_KV450 + patch + "_" + Bin + inpathP
        tmp = feather.read_dataframe(inpath)
        tmp = tmp[[para_KV450, wg_KV450]]
        data.append(tmp)
        print("Data loaded from", inpath)
    data = pd.concat(data)
    print("Patches combined.")
    # Simulated data
    inpath = inpathF_Sim + 'Selec__ZB9_in__' + Bin + inpathP
    tmp = feather.read_dataframe(inpath)
    sim = tmp[[para_Sim,wg_Sim]]
    print("Data loaded from", inpath)

    values = [data[para_KV450], sim[para_Sim]]
    wgs = [data[wg_KV450], sim[wg_Sim]]

    outpath = outpathF + Bin + outpathP
    outdirs = [outpath]
        

    HistPlotFunc(values, wgs, outdirs, 
        DENSITY, HISTTYPE, NBs, COLORs, LABELs,
        XLABEL, YLABEL)

