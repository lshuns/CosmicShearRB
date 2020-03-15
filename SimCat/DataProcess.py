#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:18:01 2019

@author: ssli

Functions dealing with the simulated catalogues

SelecFunc:
    Apply mask and re-save the original catalogues.

TomoBinFunc:
    Build and save the tomographic catalogues.

    The fiducial KV450 cosmic shear bins:
    Bin1: (0.1, 0.3]
    Bin2: (0.3, 0.5]
    Bin3: (0.5, 0.7]
    Bin4: (0.7, 0.9]
    Bin5: (0.9, 1.2]
"""

import numpy as np
import pandas as pd
import feather
from astropy.io import fits


def SelecFunc(inpath, outpath, logpath=None):
    """
    Apply mask and re-save the original catalogues.

    Parameters
    ----------
    inpath : str
        Path to the original catalogue.
    
    outpath: str
        Path to the re-saved catalogue.

    logpath: str, optional
        Path to save the log information
    """

    # with scope making the file be closed automatically
    with fits.open(inpath, memmap=True) as hdul:        
        df = pd.DataFrame(hdul[1].data)
    df = df.astype('float32')

    # selection
    # lensfit cuts with 0 means no issue and -9 means large galaxies
    fitcuts = ((df.fitclass==0) | (df.fitclass==-9))
    # # remove stars
    # star_cuts = (df.size_in>0)
    # remove invalid measurements
    weight_cuts = (df.LFweight>0)
    # valid_measurement = (df.size_out>0)
    # remove potentially blended sources
    blend_cuts = (df.contamination_radius>4.25)
    # required for meaningful T_B value
    TB_cuts = ((df.prior_matched==1) & (df.TB9_in>0))
    # remove unresolved binary stars
    binary_star_cuts = ((np.hypot(df.e1, df.e2)<=0.8) | (df.size_out>=0.5*np.exp(0.65788*(24.2-df.mag_out))))
    # remove tiny sources
    tiny_cuts = (df.size_out-df.size_corr>0.5)
    # # some random cut to the maximum snr
    # snr_cuts = (df.snr_model<210)
    #
    # # selection for full sample (more conservative)
    # df_select = df[fitcuts & star_cuts & weight_cuts & blend_cuts & TB_cuts & binary_star_cuts & tiny_cuts & snr_cuts]
    # selection used for calibration
    df_select = df[fitcuts & weight_cuts & blend_cuts & TB_cuts & binary_star_cuts & tiny_cuts]
    
    # save data
    df_select = df_select.reset_index(drop=True)
    df_select.to_feather(outpath)
    print("Selected data saved to", outpath)

    # data information    
    if logpath != None:
        # Total objects
        Ntot = len(df)
        
        # objects after selection
        Ns = len(df_select)
    
        log = open(logpath, "w")
        print("total_number,selected_number", file=log)
        print(f"{Ntot},{Ns}", file=log)
        log.close()
        print("Log information saved to", logpath)

    
def TomoBinFunc(indata, z_col, zbins_min, zbins_max, 
        outdir, save_name_prefix, logpath=None):

    """
    Build and save the tomographic catalogues.

    Parameters
    ----------
    indata : pandas.DataFrame() or numpy.recarray
        Data used for redshift binning.
    
    z_col: str
        Column name of the redshift.
    
    zbins_min: list
        Lower bounds of the redshit bins.

    zbins_max: list
        Upper bounds of the redshift bins.

    outdir: str
        Directory to the saved catalogues.

    save_name_prefix: str
        File name prefix of the saved catalogues.

    logpath: str, optional
        Path to save the log information
    """

    # data info
    if logpath != None:
        log = open(logpath, 'w')
        print("id_tomo,number", file=log)

    for i in range(len(zbins_min)):
        zmin = zbins_min[i]
        zmax = zbins_max[i]
        indataBin = indata[(indata[z_col] > zmin) & (indata[z_col] <= zmax)]
        #
        outpath = outdir + save_name_prefix + '_tomo' + str(i+1) + '.feather'
        indataBin = indataBin.reset_index(drop=True)
        indataBin.to_feather(outpath)
        print("Binned data saved in ", outpath)

        # data information
        if logpath != None:
            # selected numbers
            N_bin = len(indataBin)
            print(f"{i+1},{N_bin}", file=log)

    log.close()
    print("Log information saved to", logpath)