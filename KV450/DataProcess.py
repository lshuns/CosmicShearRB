#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:25:35 2020

@author: ssli

Functions dealing with the downloaded KV450 lensing catalogues

SelecFunc (mp supported):
    Apply mask and re-save the downloaded catalogues.

TomoBinFunc (mp supported):
    Build and save the tomographic catalogues.

    The fiducial KV450 cosmic shear bins:
    Bin1: (0.1, 0.3]
    Bin2: (0.3, 0.5]
    Bin3: (0.5, 0.7]
    Bin4: (0.7, 0.9]
    Bin5: (0.9, 1.2]

CombPatchFunc (mp supported):
    Combine all patches.

"""

import numpy as np
import pandas as pd
import feather
from astropy.io import fits

def SelecFunc(inpath, 
        outdir, save_name, pq=None):
    """
    Apply mask and re-save the downloaded catalogues.

    Parameters
    ----------
    inpath : str
        Path to the original catalogue.
    
    outdir: str
        Directory to the re-saved catalogue.

    save_name: str
        File name of the re-saved catalogue.

    pq: multiprocessing.Queue(), optional
        Object used for communication 
        between parallel processes.
    """


    # with scope making the file be closed automatically
    with fits.open(inpath, memmap=True) as hdul:        
        df = pd.DataFrame(hdul[1].data)
            
    # Drop undesired columns 
    df.drop(['THELI_NAME','2D_measurement_variance', '2D_measurement_variance_corr'], 
            axis=1, inplace=True)
    
    # 9-band mask
    # the only required mask if 9-band photo-z is used
    df_selection = df[df.GAAP_Flag_ugriZYJHKs == 0]
    
    # save data
    outpath = outdir + save_name + ".feather"
    df_selection = df_selection.reset_index(drop=True)
    df_selection.to_feather(outpath)
    
    # data information
    if pq != None:
        # Total objects
        Ntot = len(df)
        # objects after selection
        Ns = len(df_selection)
        #
        logdata = {"file_name": save_name, 'total_number': Ntot, 'selected_number': Ns}
        pq.put(logdata)


def TomoBinFunc(indata, z_col, zbins_min, zbins_max, 
        outdir, save_name_prefix, pq=None):
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

    pq: multiprocessing.Queue() object, optional
        Object used for communication 
        between parallel processes.
    """

    for i in range(len(zbins_min)):
        zmin = zbins_min[i]
        zmax = zbins_max[i]
        indataBin = indata[(indata[z_col] > zmin) & (indata[z_col] <= zmax)]
        #
        outpath = outdir + save_name_prefix + '_tomo' + str(i+1) + '.feather'
        indataBin = indataBin.reset_index(drop=True)
        indataBin.to_feather(outpath)

        # data information
        if pq != None:
            # selected numbers
            N_bin = len(indataBin)
            # data information
            logdata = {"save_name_prefix": save_name_prefix, 'id_tomo': str(i+1), 'number': N_bin}
            pq.put(logdata)


def CombPatchFunc(indir, id_tomo, patches, outdir, pq=None):
    """
    Combine all patches.

    Parameters
    ----------
    indir : str
        Directory to the catalogues to be combined.

    id_tomo: int
        Tomographic bin ID subjected to the patch combination.    

    patches: list
        A list of patch names.

    outdir: str
        Directory to the saved catalogue.

    pq: multiprocessing.Queue() object, optional
        Object used for communication 
        between parallel processes.
    """

    res = []
    for patch in patches:
        inpath = indir + patch + '_tomo' + str(id_tomo) + '.feather'
        data = feather.read_dataframe(inpath)
        res.append(data)
    res = pd.concat(res, ignore_index=True)

    # save data
    outpath = outdir + 'all_tomo' + str(id_tomo) +'.feather'
    res = res.reset_index(drop=True)
    res.to_feather(outpath)

    # data information
    if pq != None:
        N = len(res)
        #
        logdata = {"id_tomo": id_tomo, 'number': N}
        pq.put(logdata)
