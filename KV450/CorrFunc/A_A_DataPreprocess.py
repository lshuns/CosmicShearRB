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

# Parameters
# ['SeqNr', 'FLUX_AUTO_THELI', 'FLUXERR_AUTO_THELI', 'MAG_AUTO', 'MAGERR_AUTO_THELI', 'KRON_RADIUS_THELI', 
#	 'BackGr_THELI', 'Level_THELI', 'MU_THRESHOLD_THELI', 'MaxVal_THELI', 'MU_MAX_THELI', 
#	 'ISOAREA_WORLD_THELI', 'Xpos_THELI', 'Ypos_THELI', 'ALPHA_J2000', 'DELTA_J2000', 
#	 'A_WORLD', 'B_WORLD', 'THETA_J200memmap=True0', 'ERRA_WORLD_THELI', 'ERRA_IMAGE_THELI', 'ERRB_WORLD_THELI', 
#	 'ERRB_IMAGE_THELI', 'THETA_SKY_THELI', 'ERRTHETA_SKY_THELI', 'THETA_WORLD_THELI', 'ERRTHETA_WORLD_THELI', 
#	 'FWHM_IMAGE_THELI', 'FWHM_WORLD_THELI', 'Fla g_THELI', 'FLUX_RADIUS_THELI', 'NIMAFLAGS_ISO_THELI', 
#	 'CLASS_STAR_THELI', 'IMAFLAGS_ISO_THELI', 'X_WORLD_THELI', 'Y_WORLD_THELI', 'XY_WORLD_THELI', 
#	 'X2_WORLD_THELI', 'Y2_WORLD_THELI', 'ERRX2_IMAGE_THELI', 'ERRX2_WORLD_THELI', 'ERRXY_WORLD_THELI', 
#	 'ERRY2_WORLD_THELI', 'ERRXY_IMAGE_THELI', 'ERRY2_IMAGE_THELI', 'XM2_THELI', 'YM2_THELI', 'Corr_THELI', 
#	 'CXX_IMAGE_THELI', 'CXY_IMAGE_THELI', 'CYY_IMAGE_THELI', 'CXX_WORLD_THELI', 'CXY_WORLD_THELI', 'CYY_WORLD_THELI', 
#	 'ERRCXX_IMAGE_THELI', 'ERRCXY_IMAGE_THELI', 'ERRCYY_IMAGE_THELI', 'ERRCXX_WORLD_THELI', 'ERRCXY_WORLD_THELI', 
#	 'ERRCYY_WORLD_THELI', 'NPIX_THELI', 'XMIN_IMAGE_THELI', 'XMAX_IMAGE_THELI', 'YMIN_IMAGE_THELI', 'YMAX_IMAGE_THELI', 
#	 'A_THELI', 'B_THELI', 'POSANG_THELI', 'ERRTHETA_IMAGE_THELI', 'ELLIPTICITY_THELI', 'ELONGATION_THELI', 
#	 'ISO0_THELI', 'ISO1_THELI', 'ISO2_THELI', 'ISO3_THELI', 'ISO4_THELI', 'ISO5_THELI', 'ISO6_THELI', 'ISO7_THELI', 
#	 'EXTINCTION_u', 'EXTINCTION_g', 'EXTINCTION_r', 'EXTINCTION_i', 
#	 'MAG_LIM_u', 'MAG_LIM_g', 'MAG_LIM_r', 'MAG_LIM_i', 
#	 'MAG_ISO_THELI', 'MAGERR_ISO_THELI', 'MAG_ISOCOR_THELI', 'MAGERR_ISOCOR_THELI', 
#	 'MAG_BEST_THELI', 'MAGERR_BEST_THELI', 'MAG_APER_THELI', 'MAGERR_APER_THELI', 
#	 'FLUX_BEST_THELI', 'FLUXERR_BEST_THELI', 'FLUX_ISO_THELI', 'FLUXERR_ISO_THELI', 
#	 'FLUX_ISOCOR_THELI', 'FLUXERR_ISOCOR_THELI', 
#	 'FLUX_APER_4_THELI', 'FLUX_APER_5_THELI', 'FLUX_APER_6_THELI', 'FLUX_APER_7_THELI', 'FLUX_APER_8_THELI', 'FLUX_APER_9_THELI', 'FLUX_APER_10_THELI', 'FLUX_APER_11_THELI', 'FLUX_APER_12_THELI', 'FLUX_APER_13_THELI', 'FLUX_APER_14_THELI', 'FLUX_APER_15_THELI', 'FLUX_APER_16_THELI', 'FLUX_APER_17_THELI', 'FLUX_APER_18_THELI', 'FLUX_APER_19_THELI', 'FLUX_APER_20_THELI', 'FLUX_APER_25_THELI', 'FLUX_APER_30_THELI', 'FLUX_APER_35_THELI', 'FLUX_APER_40_THELI', 'FLUX_APER_45_THELI', 'FLUX_APER_50_THELI', 'FLUX_APER_55_THELI', 
#	 'FLUXERR_APER_4_THELI', 'FLUXERR_APER_5_THELI', 'FLUXERR_APER_6_THELI', 'FLUXERR_APER_7_THELI', 'FLUXERR_APER_8_THELI', 'FLUXERR_APER_9_THELI', 'FLUXERR_APER_10_THELI', 'FLUXERR_APER_11_THELI', 'FLUXERR_APER_12_THELI', 'FLUXERR_APER_13_THELI', 'FLUXERR_APER_14_THELI', 'FLUXERR_APER_15_THELI', 'FLUXERR_APER_16_THELI', 'FLUXERR_APER_17_THELI', 'FLUXERR_APER_18_THELI', 'FLUXERR_APER_19_THELI', 'FLUXERR_APER_20_THELI', 'FLUXERR_APER_25_THELI', 'FLUXERR_APER_30_THELI', 'FLUXERR_APER_35_THELI', 'FLUXERR_APER_40_THELI', 'FLUXERR_APER_45_THELI', 'FLUXERR_APER_50_THELI', 'FLUXERR_APER_55_THELI', 
#	 'Agaper_u', 'Bgaper_u', 'PAgaap_u', 'MAG_GAAP_u', 'MAGERR_GAAP_u', 'FLUX_GAAP_u', 'FLUXERR_GAAP_u', 
#	 'Agaper_g', 'Bgaper_g', 'PAgaap_g', 'MAG_GAAP_g', 'MAGERR_GAAP_g', 'FLUX_GAAP_g', 'FLUXERR_GAAP_g', 
#	 'Agaper_r', 'Bgaper_r', 'PAgaap_r', 'MAG_GAAP_r', 'MAGERR_GAAP_r', 'FLUX_GAAP_r', 'FLUXERR_GAAP_r', 
#	 'Agaper_i', 'Bgaper_i', 'PAgaap_i', 'MAG_GAAP_i', 'MAGERR_GAAP_i', 'FLUX_GAAP_i', 'FLUXERR_GAAP_i', 
#	 'Z_B_ugri', 'Z_B_MIN_ugri', 'Z_B_MAX_ugri', 'T_B_ugri', 
#	 'ODDS_ugri', 'Z_ML_ugri', 'T_ML_ugri', 'CHI_SQUARED_BPZ_ugri', 'M_0_ugri', 'FIELD_POS', 
#	 'MAG_GAAP_Z', 'MAGERR_GAAP_Z', 'FLUX_GAAP_Z', 'FLUXERR_GAAP_Z', 'GAAP_Flag_Z', 'GAAP_nexp_Z', 
#	 'MAG_GAAP_Y', 'MAGERR_GAAP_Y', 'FLUX_GAAP_Y', 'FLUXERR_GAAP_Y', 'GAAP_Flag_Y', 'GAAP_nexp_Y', 
#	 'MAG_GAAP_J', 'MAGERR_GAAP_J', 'FLUX_GAAP_J', 'FLUXERR_GAAP_J', 'GAAP_Flag_J', 'GAAP_nexp_J', 
#	 'MAG_GAAP_H', 'MAGERR_GAAP_H', 'FLUX_GAAP_H', 'FLUXERR_GAAP_H', 'GAAP_Flag_H', 'GAAP_nexp_H', 
#	 'MAG_GAAP_Ks', 'MAGERR_GAAP_Ks', 'FLUX_GAAP_Ks', 'FLUXERR_GAAP_Ks', 'GAAP_Flag_Ks', 'GAAP_nexp_Ks', 
#	 'EXTINCTION_Z', 'EXTINCTION_Y', 'EXTINCTION_J', 'EXTINCTION_H', 'EXTINCTION_Ks', 
#	 'MAG_LIM_Z', 'MAG_LIM_Y', 'MAG_LIM_J', 'MAG_LIM_H', 'MAG_LIM_Ks', 
#	 'GAAP_Flag_u', 'GAAP_Flag_g', 'GAAP_Flag_r', 'GAAP_Flag_i', 
#	 'Z_B', 'Z_B_MIN', 'ZKV450_G12_reweight_3x4x4_v2_good.cat_B_MAX', 
#	 'T_B', 'ODDS', 'Z_ML', 'T_ML', 'CHI_SQUARED_BPZ', 'M_0', 
#	 'BPZ_FILT', 'NBPZ_FILT', 'BPZ_NONDETFILT', 'NBPZ_NONDETFILT', 'BPZ_FLAGFILT', 'NBPZ_FLAGFILT', 
#	 'GAAP_Flag_ugriZYJHKs', 'wcsx', 'wcsy', 'bias_corrected_e1', 'bias_corrected_e2', 'weight', 'fitclass', 
#	 'bias_corrected_scalelength_pixels', 'bulge_fraction', 'model_flux', 'pixel_SNratio', 'model_SNratio', 
#	 'contamination_radius', 'PSF_e1', 'PSF_e2', 'PSF_Strehl_ratio', 'PSF_Q11', 'PSF_Q22', 'PSF_Q12', 
#	 'star_galaxy_f_probability', 'r_correction', '2D_measurement_variance', 'mean_likelihood_e', 
#	 'e1_correction', 'e2_correction', 'neighbour_mag', 'neighbour_distance', 
#	 'catmag', 'n_exposures_used', 'cat_ID', 
#	 'PSF_e1_exp1', 'PSF_e2_exp1', 'PSF_e1_exp2', 'PSF_e2_exp2', 'PSF_e1_exp3', 'PSF_e2_exp3', 'PSF_e1_exp4', 'PSF_e2_exp4', 'PSF_e1_exp5', 'PSF_e2_exp5', 
#	 'MASK', 'SG_FLAG', 'MAG_GAAP_u_CALIB', 'MAG_GAAP_g_CALIB', 'MAG_GAAP_r_CALIB', 'MAG_GAAP_i_CALIB', 
#	 'THELI_NAME', 'THELI_INT', 'SeqNr_field', 'weight_A', 'gmr', 'imr', 'eabs', 
#	 'binary_cut', 'recal_weight', '2D_measurement_variance_corr']
#


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
