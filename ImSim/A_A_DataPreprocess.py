#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:18:01 2019

@author: ssli

SelecFunc:
Read the original data
If it is fits, save original data as feather
Select & save for further use (feather)

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


def SelecFunc(inpath, DATAFORM, outdir, logpath):
    """
    Function for selection of targets
    """
    print("Start selection of targets (SelecFunc)...")

    if DATAFORM == 'FITS':
        # with scope making the file be closed automatically
        with fits.open(inpath, memmap=True) as hdul:        
            df = pd.DataFrame(hdul[1].data)
        outpath = outdir + "SimCatOrigi.feather"
        df = df.astype('float32')
        df = df.reset_index(drop=True)
        df.to_feather(outpath)
        print("Original data saved to", outpath)

    elif DATAFORM == 'FEATHER':
        df = feather.read_dataframe(inpath)

    else:
        raise Exception("Unsupported data form.")

    print("Loaded data from", inpath)

    # Total objects
    Ntot = len(df)
    print("Total number originally", Ntot)

    # Selection
    df = df[(df.fitclass==0) | (df.fitclass==-6) | (df.fitclass==-9)]
    df = df[(df.size_out>0) 
            & (df.contamination_radius>4.25) 
            & (df.prior_matched==1) & (df.size_in>0) 
            & (df.LFweight>0)
            & (df.TB9_in>0)]

    # objects after selection
    Ns = len(df)
    print("Total number after selection", Ns)

    # save data
    outpath = outdir + "SimCatSelec.feather"
    df = df.reset_index(drop=True)
    df.to_feather(outpath)
    print("Selected data saved to", outpath)

    # data information
    log = open(logpath, "w")
    print("TotalNumber,SelectedNumber", file=log)
    print(f"{Ntot},{Ns}", file=log)
    log.close()
    print("Log information saved to", logpath)

    print("Finished selection of targets (SelecFunc).")


def BinFunc(inpath, Z, bins, outdir, save_key, logpath):
    """
    Function for redshift binning
    """
    print("Start redshift binning (BinFunc)...")

    data = feather.read_dataframe(inpath)
    print("Loaded data from", inpath)

    # data info
    log = open(logpath, 'w')
    print("Bin,Number", file=log)

    Ntot = 0

    for i in range(len(bins)-1):
    
        zmin = bins[i] * 1e-1
        zmax = bins[i+1] * 1e-1
        dat = data[(data[Z] > zmin) & (data[Z] <= zmax)]
        print("Succeed in binning", str(bins[i]) + str(bins[i+1]))

        outpath = outdir + save_key + '__' + Z + '__' + str(bins[i]) + str(bins[i+1]) + '.feather'
        dat = dat.reset_index(drop=True)
        dat.to_feather(outpath)
        print("Data saved to", outpath)

        # selected numbers
        N_bin = len(dat)
        # Counting total number
        Ntot += N_bin

        # data information
        Bin = Z + '__' + str(bins[i]) + str(bins[i+1])
        print(f"{Bin},{N_bin}", file=log)
    print(f"AllBins,{Ntot}", file=log)
    log.close()
    print("Log information saved to", logpath)

    print("Finished redshift binning (BinFunc).")



if __name__ == "__main__":

    # ++++++++++++++++++++++++++ SelecFunc
    import time
    start_SelecFunc = time.time()

    # input
    inpath = "/disks/shear15/KiDS/ImSim/pipeline/archive/TSTnewinputglobalRecalbluered/MasterCat_TSTnewinputglobalRecalbluered_all_5_PSF.fits"
    DATAFORM = 'FITS'
    # output
    outdir = "/disks/shear15/ssli/SimCat/"
    
    # data information
    logpath = "/disks/shear15/ssli/SimCat/log/Nwhole.csv"

    SelecFunc(inpath, DATAFORM, outdir, logpath)

    print("SelecFunc finished in", time.time()-start_SelecFunc)
    # SelecFunc finished in 69.32869291305542


    # ++++++++++++++++++++++++++ BinFunc
    start_BinFunc = time.time()

    # input
    inpath = "/disks/shear15/ssli/SimCat/SimCatSelec.feather"

    # column used as redshift
    Z = 'ZB9_in'

    bins = [1, 3, 5, 7, 9, 12]

    # outpath
    outdir = "/disks/shear15/ssli/SimCat/"
    save_key = "Selec"

    ### running information
    logpath = "/disks/shear15/ssli/SimCat/log/Nbin_selec.csv"

    BinFunc(inpath, Z, bins, outdir, save_key, logpath)

    print("BinFunc finished in", time.time()-start_BinFunc)
    # BinFunc finished in 12.345204830169678

