#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 13:43:07 2020

@author: ssli

Script to run the split
"""

import multiprocessing as mp
import pandas as pd 
import numpy as np
import feather

import time

import os
import sys
# Self-defined package
sys.path.insert(0,os.path.realpath('..')) 
from Split import SplitByValFunc


# +++++++++++++++++++++++ KV450
Start = time.time()

indir = "/disks/shear15/ssli/KV450/tomo/"

nbins = 5
patches = ["all","G9","G12","G15","G23","GS"]

para_col = 'T_B'
cvalue = 3
wg_col = 'recal_weight'

# out
outdir = "/disks/shear15/ssli/KV450/split/"
logpath = '/disks/shear15/ssli/KV450/log/inf_split.csv'
log_data = open(logpath, 'w')

# for mp
jobs = []
pq = mp.Queue()

for patch in patches:
    for i in range(nbins):
        inpath = indir + patch + "_tomo" + str(i+1) + '.feather'
        data = feather.read_dataframe(inpath)
        print("Data loaded from", inpath)

        save_name_prefix = patch + '_tomo' + str(i+1)

        p = mp.Process(target=SplitByValFunc, args=(data, para_col, cvalue, outdir, save_name_prefix, wg_col, pq))
        jobs.append(p)
        p.start()
        print(f"Start running split for patch {patch} of bin {i+1}")

for p in jobs:
    p.join()


print("Finished running all split.")
print("Start saving data information...")

# data information
print("save_name_prefix,wg_tot,wg_less,wg_greater,number_tot,number_less,number_greater", file=log_data)
while not pq.empty():
    tmp = pq.get()
    print(tmp["save_name_prefix"], tmp['wg_tot'], tmp['wg_less'], tmp["wg_greater"], tmp['number_tot'], tmp['number_less'], tmp['number_greater'], sep=',', file=log_data)

log_data.close()
print("All done for selections.")

print("Finished in", time.time()-Start)
# Finished in 47.23604655265808


# +++++++++++++++++++++++ SimCat
Start = time.time()

indir = "/disks/shear15/ssli/SimCat/"

nbins = 5

para_col = 'TB9_in'
cvalue = 3
wg_col = 'LFweight'

# out
outdir = "/disks/shear15/ssli/SimCat/"
logpath = '/disks/shear15/ssli/SimCat/log/inf_split.csv'
log_data = open(logpath, 'w')

# for mp
jobs = []
pq = mp.Queue()

for i in range(nbins):
    inpath = indir + "SimCatSelec_tomo" + str(i+1) + '.feather'
    data = feather.read_dataframe(inpath)
    print("Data loaded from", inpath)

    save_name_prefix = 'SimCatSelec_tomo' + str(i+1)

    p = mp.Process(target=SplitByValFunc, args=(data, para_col, cvalue, outdir, save_name_prefix, wg_col, pq))
    jobs.append(p)
    p.start()
    print(f"Start running split for bin {i+1}")

for p in jobs:
    p.join()


print("Finished running all split.")
print("Start saving data information...")

# data information
print("save_name_prefix,wg_tot,wg_less,wg_greater,number_tot,number_less,number_greater", file=log_data)
while not pq.empty():
    tmp = pq.get()
    print(tmp["save_name_prefix"], tmp['wg_tot'], tmp['wg_less'], tmp["wg_greater"], tmp['number_tot'], tmp['number_less'], tmp['number_greater'], sep=',', file=log_data)

log_data.close()
print("All done for selections.")

print("Finished in", time.time()-Start)
# Finished in 4.039673566818237
