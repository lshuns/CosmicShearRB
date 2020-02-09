#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 15:27:14 2020

@author: ssli

Script to run the data process for downloaded KV450 catalogues
"""

import multiprocessing as mp
import numpy as np
import pandas as pd
import feather
from astropy.io import fits

import time

import os
import sys
# Self-defined package
sys.path.insert(0,os.path.realpath('..')) 
from KV450 import SelecFunc, TomoBinFunc, CombPatchFunc


# ++++++++++++++++++++++++++ selection
start_SelecFunc = time.time()

# input
inpathF = "/disks/shear15/ssli/KV450/original/KV450_"
inpathP = "_reweight_3x4x4_v2_good.cat"
patches = ["G9","G12","G15","G23","GS"]

# output
outdir = "/disks/shear15/ssli/KV450/selected/"
    
# data information
log_data = open("/disks/shear15/ssli/KV450/log/number_patch.csv", "w")

# for mp
jobs = []
pq = mp.Queue()

for patch in patches:
    inpath = inpathF + patch + inpathP

    p = mp.Process(target=SelecFunc, args=(inpath, outdir, patch, pq))
    jobs.append(p)
    p.start()
    print("Start running selection for patch", patch)

for p in jobs:
    p.join()

print("Finished running all selections.")
print("Start saving data information...")

# data information
print("patch,total_number,selected_number", file=log_data)
while not pq.empty():
    tmp = pq.get()
    print(tmp["file_name"], tmp['total_number'], tmp['selected_number'], sep=',', file=log_data)

log_data.close()
print("All done for selections.")
end_SelecFunc = time.time()


# ++++++++++++++++++++++++++ redshift binning
start_TomoBinFunc = time.time()

zbins_min = [0.1, 0.3, 0.5, 0.7, 0.9]
zbins_max = [0.3, 0.5, 0.7, 0.9, 1.2]
patches = ["G9","G12","G15","G23","GS"]
# column used as redshift
z_col = 'Z_B'

# inpath directory
inpathF = "/disks/shear15/ssli/KV450/selected/"

# outpath
outdir = "/disks/shear15/ssli/KV450/tomo/"

# running information
log_bins = open("/disks/shear15/ssli/KV450/log/number_tomo_patch.csv", "w")

# for mp
jobs = []
pq = mp.Queue()

for patch in patches:

    inpath = inpathF + patch + ".feather"
    indata = feather.read_dataframe(inpath)

    p = mp.Process(target=TomoBinFunc, args=(indata, z_col, zbins_min, zbins_max, outdir, patch, pq))
    jobs.append(p)
    p.start()
    print("Start binning in patch", patch)

for p in jobs:
    p.join()

print("Finished binning for all patches.")
print("Start saving data information...")

# data information
print("patch,id_tomo,number", file=log_bins)
while not pq.empty():
    tmp = pq.get()
    print(tmp["save_name_prefix"], tmp['id_tomo'], tmp['number'], sep=',', file=log_bins)

log_bins.close()
print("All done for binning.")

end_TomoBinFunc = time.time()


# ++++++++++++++++++++++++++ Combine all the patches
start_CombPatchFunc = time.time()

indir = "/disks/shear15/ssli/KV450/tomo/"
outdir = "/disks/shear15/ssli/KV450/tomo/"
patches =  ["G9","G12","G15","G23","GS"]
nbins = 5

outpath_inf = "/disks/shear15/ssli/KV450/log/number_tomo.csv"
aORw = 'w'


# running information
log_patches = open("/disks/shear15/ssli/KV450/log/number_tomo.csv", "w")


# for mp
jobs = []
pq = mp.Queue()

for i in range(nbins):

    p = mp.Process(target=CombPatchFunc, args=(indir, i+1, patches, outdir, pq))
    jobs.append(p)
    p.start()
    print("Start combining patches for bin", str(i+1))
    
for p in jobs:
    p.join()

print("Finished combining for all bins.")
print("Start saving data information...")

# data information
print("id_tomo,number", file=log_patches)
while not pq.empty():
    tmp = pq.get()
    print(tmp["id_tomo"], tmp['number'], sep=',', file=log_patches)

log_patches.close()

print("All done for combining patches.")

end_CombPatchFunc = time.time()

# +++++++++++++++++++++++++ Running time
print("Running time for SelecFunc", end_SelecFunc-start_SelecFunc, "seconds.")
print("Running time for TomoBinFunc", end_TomoBinFunc-start_TomoBinFunc, "seconds.")
print("Running time for CombPatchFunc", end_CombPatchFunc-start_CombPatchFunc, "seconds.")
# Running on markermeer
# Running time for SelecFunc 114.11349105834961 seconds.
# Running time for TomoBinFunc 43.105056285858154 seconds.
# Running time for CombPatchFunc 57.36755561828613 seconds.
