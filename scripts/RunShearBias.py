#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 17:49:17 2020

@author: ssli

script to run the shear bias calibration
"""

import multiprocessing as mp
import numpy as np
from scipy import optimize
import pandas as pd
from astropy.io import fits
import feather

import time

import os
import sys
# Self-defined package
sys.path.insert(0,os.path.realpath('..')) 
from ShearBias import mcCalFunc


Start = time.time()

# directory to SimCat
inDirSim = "/disks/shear15/ssli/SimCat/"
# directory to KV450 data
inDirReal = "/disks/shear15/ssli/KV450/"

# Number of tomographic bins
Nbins = 5

# Number of re-weighting bins
Nbin1 = 20
Nbin2 = 20

# output
outpath_whole = "/disks/shear15/ssli/CosmicShear/shear_bias/Summary_multiplicative_whole.csv"
outpath_red = "/disks/shear15/ssli/CosmicShear/shear_bias/Summary_multiplicative_red.csv"
outpath_blue = "/disks/shear15/ssli/CosmicShear/shear_bias/Summary_multiplicative_blue.csv"
#
outfile_whole = open(outpath_whole, "w")
print("bin,m,m_err_BS,m1,m2,m1_err,m2_err,m1_err_BS,m2_err_BS", file=outfile_whole)
outfile_red = open(outpath_red, "w")
print("bin,m,m_err_BS,m1,m2,m1_err,m2_err,m1_err_BS,m2_err_BS", file=outfile_red)
outfile_blue = open(outpath_blue, "w")
print("bin,m,m_err_BS,m1,m2,m1_err,m2_err,m1_err_BS,m2_err_BS", file=outfile_blue)


# for mp
jobs_whole = []
pq_whole = mp.Queue()

jobs_red = []
pq_red = mp.Queue()

jobs_blue = []
pq_blue = mp.Queue()


# MC calculation
for i in range(Nbins):
# for i in range(1):
    # simulated data
    inpathSim_whole = inDirSim + "SimCatSelec_tomo" + str(i+1) +'.feather'
    inpathSim_red = inDirSim + "SimCatSelec_tomo" + str(i+1) +'_TB9_in_less3.feather'
    inpathSim_blue = inDirSim + "SimCatSelec_tomo" + str(i+1) +'_TB9_in_greater3.feather'
    #
    dataSim_whole = feather.read_dataframe(inpathSim_whole)
    dataSim_red = feather.read_dataframe(inpathSim_red)
    dataSim_blue = feather.read_dataframe(inpathSim_blue)
    
    # real data
    inpathReal_whole = inDirReal + 'tomo/all_tomo' + str(i+1) +'.feather'
    inpathReal_red = inDirReal + 'split/all_tomo' + str(i+1) +'_T_B_less3.feather'
    inpathReal_blue = inDirReal + 'split/all_tomo' + str(i+1) +'_T_B_greater3.feather'
    #        
    dataReal_whole = feather.read_dataframe(inpathReal_whole)
    dataReal_red = feather.read_dataframe(inpathReal_red)
    dataReal_blue = feather.read_dataframe(inpathReal_blue)

    p_whole = mp.Process(target=mcCalFunc, args=(dataSim_whole, dataReal_whole, Nbin1, Nbin2, pq_whole))
    p_red = mp.Process(target=mcCalFunc, args=(dataSim_red, dataReal_red, Nbin1, Nbin2, pq_red))
    p_blue = mp.Process(target=mcCalFunc, args=(dataSim_blue, dataReal_blue, Nbin1, Nbin2, pq_blue))

    jobs_whole.append(p_whole)
    p_whole.start()
    print("Start running for bin", str(i+1), "in whole data.")

    jobs_red.append(p_red)
    p_red.start()
    print("Start running for bin", str(i+1), "in red data.")

    jobs_blue.append(p_blue)
    p_blue.start()
    print("Start running for bin", str(i+1), "in blue data.")


for p_whole in jobs_whole:
    p_whole.join()
for p_red in jobs_red:
    p_red.join()
for p_blue in jobs_blue:
    p_blue.join()

print("Finished running all.")
print("Start saving data information...")

i = 1
while not pq_whole.empty():
    tmp = pq_whole.get()
    print(i, tmp["m_final"], tmp['m_err_BS_final'], \
        tmp['m1_final'], tmp['m2_final'], \
        tmp['m1_err_final'], tmp['m2_err_final'], \
        tmp['m1_err_BS_final'], tmp['m2_err_BS_final'], \
        sep=',', file=outfile_whole)
    i += 1
outfile_whole.close()

i = 1
while not pq_red.empty():
    tmp = pq_red.get()
    print(i, tmp["m_final"], tmp['m_err_BS_final'], \
        tmp['m1_final'], tmp['m2_final'], \
        tmp['m1_err_final'], tmp['m2_err_final'], \
        tmp['m1_err_BS_final'], tmp['m2_err_BS_final'], \
        sep=',', file=outfile_red)
    i += 1
outfile_red.close()

i = 1
while not pq_blue.empty():
    tmp = pq_blue.get()
    print(i, tmp["m_final"], tmp['m_err_BS_final'], \
        tmp['m1_final'], tmp['m2_final'], \
        tmp['m1_err_final'], tmp['m2_err_final'], \
        tmp['m1_err_BS_final'], tmp['m2_err_BS_final'], \
        sep=',', file=outfile_blue)
    i += 1
outfile_blue.close()

print("All Finished in", time.time()-Start)
# All Finished in 147.10100865364075
