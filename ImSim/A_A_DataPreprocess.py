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


def SelecFunc(inpath, DATAFORM, SELECTION, outdir, outkey, logpath):
    """
    Function for selection of targets
    """
    print("Start selection of targets (SelecFunc)...")
    print("Selection:", SELECTION)

    if DATAFORM == 'FITS':
        # with scope making the file be closed automatically
        with fits.open(inpath, memmap=True) as hdul:        
            df = pd.DataFrame(hdul[1].data)
        df = df.astype('float32')

    elif DATAFORM == 'FEATHER':
        df = feather.read_dataframe(inpath)

    else:
        raise Exception("Unsupported data form!")

    print("Loaded data from", inpath)

    # Total objects
    Ntot = len(df)
    print("Total number originally", Ntot)

    # Selection
    if SELECTION == 'NONE':
        outpath = outdir + "SimCatOrigi_" + outkey + ".feather"
        df = df.reset_index(drop=True)
        df.to_feather(outpath)
        print("Original data saved to", outpath)
        print("Finished selection of targets (SelecFunc).")
        return 0

    elif SELECTION == 'ARUN':
        fitcuts = ((df.fitclass==0) | (df.fitclass==-6) | (df.fitclass==-9))
        star_cuts = (df.size_in>0)
        weight_cuts = (df.LFweight>0) ## because measurements must be valid
        valid_measurement = (df.size_out>0)
        blend_cuts = (df.contamination_radius>4.25)

        binary_star_cuts = ((np.hypot(df.e1, df.e2)<=0.8) | (df.size_out>=0.5*np.exp(0.65788*(24.2-df.mag_out))))
        tiny_cuts = (df.size_out-df.size_corr>0.5)
        snr_cuts = (df.snr_model<210)

        df = df[fitcuts & star_cuts & weight_cuts & valid_measurement & blend_cuts & binary_star_cuts & tiny_cuts & snr_cuts]
        outpath = outdir + "SimCatSelec_Arun_" + outkey + ".feather"

    elif SELECTION == 'NEW':
        fitcuts = ((df.fitclass==0) | (df.fitclass==-6) | (df.fitclass==-9))
        star_cuts = (df.size_in>0)
        weight_cuts = (df.LFweight>0) ## because measurements must be valid
        valid_measurement = (df.size_out>0)
        blend_cuts = (df.contamination_radius>4.25)
        TB_cuts = ((df.prior_matched==1) & (df.TB9_in>0))

        binary_star_cuts = ((np.hypot(df.e1, df.e2)<=0.8) | (df.size_out>=0.5*np.exp(0.65788*(24.2-df.mag_out))))
        tiny_cuts = (df.size_out-df.size_corr>0.5)
        snr_cuts = (df.snr_model<210)

        df = df[fitcuts & star_cuts & weight_cuts & valid_measurement & blend_cuts & TB_cuts & binary_star_cuts & tiny_cuts & snr_cuts]
        outpath = outdir + "SimCatSelec_" + outkey + ".feather"

    else:
        raise Exception("Unsupported selection criterion!")

    # objects after selection
    Ns = len(df)
    print("Total number after selection", Ns)

    # save data
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
    return 1

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

    # # ++++++++++++++++++++++++++ SelecFunc
    import time

    # ++++ from fits to feather
    start_SelecFunc = time.time()
    # input
    inpath = "/disks/shear15/KiDS/ImSim/pipeline/archive/TSTnewinputglobalRecalbluered/MasterCat_TSTnewinputglobalRecalbluered_all_13_PSF.fits"
    # Total number originally 37470125

    outkey = "all_13_PSF"
    DATAFORM = 'FITS'
    SELECTION = 'NONE'
    # output
    outdir = "/disks/shear15/ssli/SimCat/"

    # data information
    logpath = "None"

    SelecFunc(inpath, DATAFORM, SELECTION, outdir, outkey, logpath)
    print("SelecFunc (NONE) finished in", time.time()-start_SelecFunc)
    # finished in 98.01019883155823

    # ++++ Arun 
    # # 5 PSF
    # start_SelecFunc = time.time()
    # # input
    # outkey = "all_5_PSF"
    # # inpath = "/disks/shear15/ssli/SimCat/SimCatOrigi_"+outkey+".feather"
    # inpath = "/disks/shear15/KiDS/ImSim/pipeline/archive/TSTnewinputglobalRecal/MasterCat_TSTnewinputglobalRecal_all_5_PSF.fits"
    # # Total number originally 14345700
    # # Total number after selection 5848519
    # DATAFORM = 'FITS'
    # SELECTION = 'ARUN'
    # # output
    # outdir = "/disks/shear15/ssli/SimCat/"
    
    # # data information
    # logpath = "/disks/shear15/ssli/SimCat/log/Nwhole_Arun_"+outkey+".csv"

    # SelecFunc(inpath, DATAFORM, SELECTION, outdir, outkey, logpath)

    # print("SelecFunc (Arun) finished in", time.time()-start_SelecFunc)
    # # finished in 20.329578638076782

    # 13 PSF
    start_SelecFunc = time.time()
    # input
    outkey = "all_13_PSF"
    inpath = "/disks/shear15/ssli/SimCat/SimCatOrigi_"+outkey+".feather"
    # inpath = "/disks/shear15/KiDS/ImSim/pipeline/archive/TSTnewinputglobalRecal/MasterCat_TSTnewinputglobalRecal_all_13_PSF.fits"
    # Total number originally 37470125
    # Total number after selection 15140564 
    DATAFORM = 'FEATHER'
    SELECTION = 'ARUN'
    # output
    outdir = "/disks/shear15/ssli/SimCat/"
    
    # data information
    logpath = "/disks/shear15/ssli/SimCat/log/Nwhole_Arun_"+outkey+".csv"

    SelecFunc(inpath, DATAFORM, SELECTION, outdir, outkey, logpath)
    print("SelecFunc (Arun) finished in", time.time()-start_SelecFunc)
    # finished in 20.448835849761963


    # # ++++ new
    # 13 PSF
    start_SelecFunc = time.time()
    # input
    outkey = "all_13_PSF"
    inpath = "/disks/shear15/ssli/SimCat/SimCatOrigi_"+outkey+".feather"
    # inpath = "/disks/shear15/KiDS/ImSim/pipeline/archive/TSTnewinputglobalRecal/MasterCat_TSTnewinputglobalRecal_all_13_PSF.fits"
    # Total number originally 37470125
    # Total number after selection 13128786
    DATAFORM = 'FEATHER'
    SELECTION = 'NEW'
    # output
    outdir = "/disks/shear15/ssli/SimCat/"
    
    # data information
    logpath = "/disks/shear15/ssli/SimCat/log/Nwhole_new_"+outkey+".csv"

    SelecFunc(inpath, DATAFORM, SELECTION, outdir, outkey, logpath)
    print("SelecFunc (new) finished in", time.time()-start_SelecFunc)
    # finished in 16.999256372451782


    # ++++++++++++++++++++++++++ BinFunc
    # Arun
    start_BinFunc = time.time()
    # input
    outkey = "Arun_all_13_PSF"
    inpath = "/disks/shear15/ssli/SimCat/SimCatSelec_"+outkey+".feather"

    # column used as redshift
    Z = 'ZB9_in'

    bins = [1, 3, 5, 7, 9, 12]

    # outpath
    outdir = "/disks/shear15/ssli/SimCat/"
    save_key = "Bin_"+outkey

    ### running information
    logpath = "/disks/shear15/ssli/SimCat/log/Nbin_"+outkey+".csv"

    BinFunc(inpath, Z, bins, outdir, save_key, logpath)

    print("BinFunc finished in", time.time()-start_BinFunc)
    # BinFunc finished in 14.208211183547974


    # new
    start_BinFunc = time.time()
    # input
    outkey = "all_13_PSF"
    inpath = "/disks/shear15/ssli/SimCat/SimCatSelec_"+outkey+".feather"

    # column used as redshift
    Z = 'ZB9_in'

    bins = [1, 3, 5, 7, 9, 12]

    # outpath
    outdir = "/disks/shear15/ssli/SimCat/"
    save_key = "Bin_"+outkey

    ### running information
    logpath = "/disks/shear15/ssli/SimCat/log/Nbin_"+outkey+".csv"

    BinFunc(inpath, Z, bins, outdir, save_key, logpath)

    print("BinFunc finished in", time.time()-start_BinFunc)
    # BinFunc finished in 11.965852499008179

