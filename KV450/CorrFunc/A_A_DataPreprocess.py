#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 10:20:07 2019

@author: ssli

SelecFunc:
Read the original data (fits)
Set selection criteria
Save for further use (feather)

BinFunc:
build the five tomographic redshift bins:
    "13": (0.1, 0.3]
    "35": (0.3, 0.5]
    "57": (0.5, 0.7]
    "79": (0.7, 0.9]
    "912": (0.9, 1.2]

BinCountFunc:
Function for number counting for binning results
"""

import numpy as np
import pandas as pd
import feather
from astropy.io import fits

# def setup():
#     """
#     CallFunction
#     """

def SelecFunc(patch, inpath, outdir, pq):
    """
    Function for selection of targets
    """
    print(patch, "Start selection of targets (SelecFunc)...")
        
    # with scope making the file be closed automatically
    with fits.open(inpath, memmap=True) as hdul:        
        df = pd.DataFrame(hdul[1].data)

    print(patch, "Loaded data from", inpath)
            
    # Drop undesired columns 
    df.drop(['THELI_NAME','2D_measurement_variance', '2D_measurement_variance_corr'], 
            axis=1, inplace=True)
    
    print(patch, "Undesired columns dropped.")

    # Total objects
    Ntot = len(df)

    # Selection
    df = df[df.GAAP_Flag_ugriZYJHKs == 0]
    
    print(patch, "Succeed in selection.")

    # objects after selection
    Ns = len(df)

    # save data
    outpath = outdir + patch + ".feather"
    df = df.reset_index(drop=True)
    df.to_feather(outpath)
    print(patch, "Data saved to", outpath)

    # data information
    logdata = {"patch": patch, 'totalNumber': Ntot, 'selectedNumber': Ns}
    pq.put(logdata)

    print(patch, "Finished selection of targets (SelecFunc).")


def BinFunc(patch, inpath, outdir, Z, bins, pq):
    """
    Function for redshift binning
    """
    print(patch, "Start redshift binning (BinFunc)...")

    data = feather.read_dataframe(inpath)
    print(patch, "Loaded data from", inpath)

    for i in range(len(bins)-1):
    
        zmin = bins[i] * 1e-1
        zmax = bins[i+1] * 1e-1
        dat = data[(data[Z] > zmin) & (data[Z] <= zmax)]
        print(patch, "Succeed in binning", str(bins[i]) + str(bins[i+1]))

        outpath = outdir + patch + '__' + Z + '__' + str(bins[i]) + str(bins[i+1]) + '.feather'
        dat = dat.reset_index(drop=True)
        dat.to_feather(outpath)
        print(patch, "Data saved to", outpath)

        # selected numbers
        N_bin = len(dat)

        # data information
        save_key = Z + '__' + str(bins[i]) + str(bins[i+1])
        logdata = {"patch": patch, 'bin': save_key, 'number': N_bin}
        pq.put(logdata)

    print(patch, "Finished redshift binning (BinFunc).")


def BinCountFunc(infpath, outpath, Z, bins):
    """
    Function for number counting for binning results
    """
    print("Start number counting (BinCountFunc)...")

    data = pd.read_csv(infpath)

    # data information
    log = open(outpath, "w")
    print("bin,number", file=log)

    res = []
    for i in range(len(bins)-1):
        key = Z + '__' + str(bins[i]) + str(bins[i+1])
        
        N = np.sum(data[data.bin==key]['number'].values)
        print("Number in bin (", str(bins[i] * 1e-1), ",", str(bins[i+1] * 1e-1), "]:", N)
        print(key, N, sep=',', file=log)

        res.append(N)

    N_tot = np.sum(res)
    print("Total Number in all bins:", N_tot)
    print('total', N_tot, sep=',', file=log)
    log.close()

    print("Finished number counting (BinCountFunc)...")



if __name__ == "__main__":
    import multiprocessing as mp

    # ++++++++++++++++++++++++++ SelecFunc
    import time
    start_SelecFunc = time.time()

    # input
    inpathF = "/disks/shear15/ssli/KV450/KV450_"
    inpathP = "_reweight_3x4x4_v2_good.cat"
    pathGL = ["G9","G12","G15","G23","GS"]

    # output
    outdir = "/disks/shear15/ssli/KV450/selected/"
    
    # data information
    log_data = open("/disks/shear15/ssli/KV450/selected/log/log_data.csv", "w")

    # for mp
    jobs = []
    pq = mp.Queue()

    for s in pathGL:
        inpath = inpathF + s + inpathP

        p = mp.Process(target=SelecFunc, args=(s, inpath, outdir, pq))
        jobs.append(p)
        p.start()

    for p in jobs:
        p.join()

    print("Start saving data information.")

    # data information
    print("patch,totalNumber,selectedNumber", file=log_data)
    while not pq.empty():
        tmp = pq.get()
        print(tmp["patch"], tmp['totalNumber'], tmp['selectedNumber'], sep=',', file=log_data)

    log_data.close()
    print("Done for SelecFunc.")
    end_SelecFunc = time.time()


    # ++++++++++++++++++++++++++ BinFunc
    start_BinFunc = time.time()

    bins = [1, 3, 5, 7, 9, 12]
    patches = ["G9","G12","G15","G23","GS"]
    # column used as redshift
    Z = 'Z_B'

    # inpath directory
    inpathF = "/disks/shear15/ssli/KV450/selected/"

    # outpath
    outdir = "/disks/shear15/ssli/KV450/selected/"

    ### running information
    log_bins = open("/disks/shear15/ssli/KV450/selected/log/log_bins.csv", "w")

    # for mp
    jobs = []
    pq = mp.Queue()

    for patch in patches:
        inpath = inpathF + patch + ".feather"

        p = mp.Process(target=BinFunc, args=(patch, inpath, outdir, Z, bins, pq))
        jobs.append(p)
        p.start()

    for p in jobs:
        p.join()

    print("Start saving data information.")

    # data information
    print("patch,bin,number", file=log_bins)
    while not pq.empty():
        tmp = pq.get()
        print(tmp["patch"], tmp['bin'], tmp['number'], sep=',', file=log_bins)

    log_bins.close()
    print("Done for BinFunc.")

    end_BinFunc = time.time()

    # ++++++++++++++++++++++++++ BinCountFunc
    start_BinCountFunc = time.time()

    bins = [1, 3, 5, 7, 9, 12]
    # column used as redshift
    Z = 'Z_B'

    # path directory
    infpath = "/disks/shear15/ssli/KV450/selected/log/log_bins.csv"

    # outpath
    outpath = "/disks/shear15/ssli/KV450/selected/log/log_BinCount.csv"

    BinCountFunc(infpath, outpath, Z, bins)

    end_BinCountFunc = time.time()

    # +++++++++++++++++++++++++ Running time
    print("Running time for SelecFunc", end_SelecFunc-start_SelecFunc, "seconds.")
    print("Running time for BinFunc", end_BinFunc-start_BinFunc, "seconds.")
    print("Running time for BinCountFunc", end_BinCountFunc-start_BinCountFunc, "seconds.")
