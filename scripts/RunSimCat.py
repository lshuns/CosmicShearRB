#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 19:21:53 2020

@author: ssli

Script to run the data process for simulated catalogue

Package:
    SimCat:
        SelecFunc: selection
        TomoBinFunc: redshift binning
Data location: SimCat
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

# ++++++++++++++++++++++++++ General information
# Parent path
data_path = "/disks/shear15/ssli/SimCat/"
zbins_min = [0.1, 0.3, 0.5, 0.7, 0.9]
zbins_max = [0.3, 0.5, 0.7, 0.9, 1.2]
# column used as redshift
z_col = 'ZB9_in'
# prefix for all the saved data
save_name_prefix = "SimCatSelec"


# ++++++++++++++++++++++++++ Selection
start_SelecFunc = time.time() 

# input
inpath = data_path + "MasterCat_TSTnewinputglobalRecalbluered_all_13_PSF.fits"

# output
outpath = data_path + save_name_prefix + ".feather"

# data information
logpath = data_path + "log/number.csv"

SelecFunc(inpath, outpath, logpath)
end_SelecFunc = time.time()

# ++++++++++++++++++++++++++ Redshift binning
start_BinFunc = time.time()

# input
inpath = data_path + save_name_prefix + ".feather"
indata = feather.read_dataframe(inpath)

# outpath
outdir = data_path

# running information
logpath = data_path + "log/number_tomo.csv"

TomoBinFunc(indata, z_col, zbins_min, zbins_max, 
        outdir, save_name_prefix, logpath)

end_BinFunc = time.time()


print("SelecFunc finished in", end_SelecFunc-start_SelecFunc)
print("TomoBinFunc finished in", end_BinFunc-start_BinFunc)
# eemmeer (March 15, 2020)
# SelecFunc finished in 92.3946418762207
# TomoBinFunc finished in 12.4992196559906