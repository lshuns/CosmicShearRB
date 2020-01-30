#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:25:35 2020

@author: ssli

Functions dealing with the KV450 data

SelecFunc (mp supported):
    set selection and re-save the original data

BinFunc (mp supported):
    binning data and save
BinCountFunc
    Bin number counts
        associated with BinFunc with pq!=None

PatchFunc:
    combining all patches

"""

import numpy as np
import pandas as pd
import feather
from astropy.io import fits

def SelecFunc(inpath, outdir, save_key, pq=None):
    """
    Function for selection of targets
    """

    # with scope making the file be closed automatically
    with fits.open(inpath, memmap=True) as hdul:        
        df = pd.DataFrame(hdul[1].data)
            
    # Drop undesired columns 
    df.drop(['THELI_NAME','2D_measurement_variance', '2D_measurement_variance_corr'], 
            axis=1, inplace=True)
    
    # Selection
    df_selection = df[df.GAAP_Flag_ugriZYJHKs == 0]
    
    # save data
    outpath = outdir + save_key + ".feather"
    df_selection = df_selection.reset_index(drop=True)
    df_selection.to_feather(outpath)
    
    # data information
    if pq != None:
        # Total objects
        Ntot = len(df)
        # objects after selection
        Ns = len(df_selection)
        #
        logdata = {"save_key": save_key, 'totalNumber': Ntot, 'selectedNumber': Ns}
        pq.put(logdata)


def BinFunc(indata, Z, bins, outdir, save_key, pq=None):
    """
    Function for redshift binning
    """

    for i in range(len(bins)-1):
        zmin = bins[i]
        zmax = bins[i+1]
        indata = indata[(indata[Z] > zmin) & (indata[Z] <= zmax)]
        #
        outpath = outdir + save_key + '_bin' + str(i) + '.feather'
        indata = indata.reset_index(drop=True)
        indata.to_feather(outpath)


        if pq != None:
            # selected numbers
            N_bin = len(indata)
            # data information
            logdata = {"save_key": save_key, 'bin': str(i), 'number': N_bin}
            pq.put(logdata)

def BinCountFunc(infpath, Nbins, 
                   outprint=None, outpath=None):
    """
    Function for number counts
        associated with BinFunc with pq!=None
    """

    data = pd.read_csv(infpath)

    # data information
    if outpath != None:
        log = open(outpath, "w")
        print("bin,number", file=log)

    res = []
    for i in range(Nbins):
        key = i       
        N = np.sum(data[data.bin==key]['number'].values)
        if outprint != None:
            print(f"Number in bin {key}: {N}")
        if outpath != None:
            print(key, N, sep=',', file=log)

        res.append(N)

    N_tot = np.sum(res)
    if outprint != None:
        print(f"Number in total bins: {N_tot}")
    if outpath != None:
        print('total', N_tot, sep=',', file=log)
        log.close()



def PatchFunc(inDir, patches, read_key, outpath):
    """
    Function for combining all patches
    """

    res = []
    for patch in patches:
        inpath = inDir + patch + read_key
        data = feather.read_dataframe(inpath)
        res.append(data)
    res = pd.concat(res, ignore_index=True)

    # save data
    res = res.reset_index(drop=True)
    res.to_feather(outpath)
    

if __name__ == "__main__":

    import multiprocessing as mp

    # ++++++++++++++++++++++++++ SelecFunc
    import time
    start_SelecFunc = time.time()

    # input
    inpathF = "/disks/shear15/ssli/KV450/data/KV450_"
    inpathP = "_reweight_3x4x4_v2_good.cat"
    pathGL = ["G9","G12","G15","G23","GS"]

    # output
    outdir = "/disks/shear15/ssli/KV450/data/feather/"
    
    # data information
    log_data = open("/disks/shear15/ssli/KV450/data/log/numberPatches.csv", "w")

    # for mp
    jobs = []
    pq = mp.Queue()

    for s in pathGL:
        inpath = inpathF + s + inpathP
        save_key = s + '_selected'

        p = mp.Process(target=SelecFunc, args=(inpath, outdir, save_key, pq))
        jobs.append(p)
        p.start()

    for p in jobs:
        p.join()

    print("Start saving data information.")

    # data information
    print("save_key,totalNumber,selectedNumber", file=log_data)
    while not pq.empty():
        tmp = pq.get()
        print(tmp["save_key"], tmp['totalNumber'], tmp['selectedNumber'], sep=',', file=log_data)

    log_data.close()
    print("Done for SelecFunc.")
    end_SelecFunc = time.time()


    # ++++++++++++++++++++++++++ BinFunc
    start_BinFunc = time.time()

    bins = [0.1, 0.3, 0.5, 0.7, 0.9, 1.2]
    patches = ["G9","G12","G15","G23","GS"]
    # column used as redshift
    Z = 'Z_B'

    # inpath directory
    inpathF = "/disks/shear15/ssli/KV450/data/feather/"

    # outpath
    outdir = "/disks/shear15/ssli/KV450/data/feather/"

    ### running information
    log_bins = open("/disks/shear15/ssli/KV450/data/log/numberBinsPatch.csv", "w")

    # for mp
    jobs = []
    pq = mp.Queue()


    for patch in patches:
        save_key = patch + '_selected'

        inpath = inpathF + save_key + ".feather"
        indata = feather.read_dataframe(inpath)

        p = mp.Process(target=BinFunc, args=(indata, Z, bins, outdir, save_key, pq))
        jobs.append(p)
        p.start()

    for p in jobs:
        p.join()

    print("Start saving data information.")

    # data information
    print("save_key,bin,number", file=log_bins)
    while not pq.empty():
        tmp = pq.get()
        print(tmp["save_key"], tmp['bin'], tmp['number'], sep=',', file=log_bins)

    log_bins.close()
    print("Done for BinFunc.")

    end_BinFunc = time.time()


    # ++++++++++++++++++++++++++ PatchFunc
    start_PatchFunc = time.time()

    inDir = "/disks/shear15/ssli/KV450/data/feather/"
    outDir = "/disks/shear15/ssli/KV450/data/feather/"
    patches =  ["G9","G12","G15","G23","GS"]
    Nbins = 5

    for i in range(Nbins):
        read_key = '_selected_bin' + str(i) + '.feather'
        outpath = outDir + 'AllPatch_selected_bin' +str(i) +'.feather'

        PatchFunc(inDir, patches, read_key, outpath)

    end_PatchFunc = time.time()


    # ++++++++++++++++++++++++++ BinCountFunc
    start_BinCountFunc = time.time()

    Nbins = 5

    # path directory
    infpath = "/disks/shear15/ssli/KV450/data/log/numberBinsPatch.csv"

    # outpath
    outpath = "/disks/shear15/ssli/KV450/data/log/numberBins.csv"

    BinCountFunc(infpath, Nbins, 
                   outprint=None, outpath=outpath)
    
    end_BinCountFunc = time.time()

    # +++++++++++++++++++++++++ Running time
    print("Running time for SelecFunc", end_SelecFunc-start_SelecFunc, "seconds.")
    print("Running time for BinFunc", end_BinFunc-start_BinFunc, "seconds.")
    print("Running time for PatchFunc", end_PatchFunc-start_PatchFunc, "seconds.")
    print("Running time for BinCountFunc", end_BinCountFunc-start_BinCountFunc, "seconds.")

    # Running on markermeer
# Running time for SelecFunc 263.9306733608246 seconds.
# Running time for BinFunc 20.77444362640381 seconds.
# Running time for PatchFunc 26.071428060531616 seconds.
# Running time for BinCountFunc 0.01830768585205078 seconds.
