#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 19:21:53 2020

@author: ssli

Script to run the data process for simulated catalogue
"""

import numpy as np
import pandas as pd
import feather
from astropy.io import fits

import time

import os
import sys
# Self-defined package
sys.path.insert(0,os.path.realpath('..')) 
from SimCat import SelecFunc, TomoBinFunc




# ++++++++++++++++++++++++++ Selection
start_SelecFunc = time.time() 

# input
inpath = "/disks/shear15/ssli/SimCat/MasterCat_TSTnewinputglobalRecalbluered_all_13_PSF.fits"

# output
outpath = "/disks/shear15/ssli/SimCat/SimCatSelec.feather"

# data information
logpath = "/disks/shear15/ssli/SimCat/log/number.csv"

SelecFunc(inpath, outpath, logpath)

print("SelecFunc finished in", time.time()-start_SelecFunc)
# SelecFunc finished in 106.43442273139954


# ++++++++++++++++++++++++++ Binning
start_BinFunc = time.time()

zbins_min = [0.1, 0.3, 0.5, 0.7, 0.9]
zbins_max = [0.3, 0.5, 0.7, 0.9, 1.2]
patches = ["G9","G12","G15","G23","GS"]
# column used as redshift
z_col = 'ZB9_in'

# input
inpath = "/disks/shear15/ssli/SimCat/SimCatSelec.feather"
indata = feather.read_dataframe(inpath)

# outpath
outdir = "/disks/shear15/ssli/SimCat/"
save_name_prefix = "SimCatSelec"

# running information
logpath = "/disks/shear15/ssli/SimCat/log/number_tomo.csv"

TomoBinFunc(indata, z_col, zbins_min, zbins_max, 
        outdir, save_name_prefix, logpath)

print("TomoBinFunc finished in", time.time()-start_BinFunc)
# TomoBinFunc finished in 11.625519037246704
