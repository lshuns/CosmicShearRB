#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 13:43:07 2020

@author: ssli

Script to run the split for KV450 data & SimCat

Package:
    Split: SplitByValFunc

Data location:
    KV450:
        Input: tomo
        Output: split, log
    SimCat
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

# +++++++++++++++++++++++ General information
nbins = 5
cvalue = 3

data_path_KV450 = "/disks/shear15/ssli/KV450/"
patches = ["all","G9","G12","G15","G23","GS"]
para_col_KV450 = 'T_B'
wg_col_KV450 = 'recal_weight'

data_path_SimCat = "/disks/shear15/ssli/SimCat/"
para_col_SimCat = 'TB9_in'
wg_col_SimCat = 'LFweight'


# +++++++++++++++++++++++ KV450
start_KV450 = time.time()

# input
indir = data_path_KV450 + "tomo/"

# out
outdir = data_path_KV450 + "split/"
logpath = data_path_KV450 + 'log/inf_split.csv'
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

        p = mp.Process(target=SplitByValFunc, args=(data, para_col_KV450, cvalue, outdir, save_name_prefix, wg_col_KV450, pq))
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
print("All done for selection in KV450 data.")
end_KV450 = time.time()

# +++++++++++++++++++++++ SimCat
start_SimCat = time.time()

indir = data_path_SimCat

# out
outdir = data_path_SimCat
logpath = data_path_SimCat + 'log/inf_split.csv'
log_data = open(logpath, 'w')

# for mp
jobs = []
pq = mp.Queue()

for i in range(nbins):
    inpath = indir + "SimCatSelec_tomo" + str(i+1) + '.feather'
    data = feather.read_dataframe(inpath)
    print("Data loaded from", inpath)

    save_name_prefix = 'SimCatSelec_tomo' + str(i+1)

    p = mp.Process(target=SplitByValFunc, args=(data, para_col_SimCat, cvalue, outdir, save_name_prefix, wg_col_SimCat, pq))
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
print("All done for selection in SimCat.")
end_SimCat = time.time()


print("KV450 finished in", end_KV450-start_KV450)
print("SimCat finished in", end_SimCat-start_SimCat)
# eemmeer (March 15, 2020)
# KV450 finished in 41.25113773345947
# SimCat finished in 4.159558057785034