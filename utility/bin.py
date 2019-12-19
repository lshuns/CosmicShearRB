#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 15:21:31 2019

@author: ssli

Functions for binning

SimpleBinFunc:
(Not tested)
simple binning given bounds

WgQuantile1DFunc:
compute the weighted quantile of a 1D numpy array

WgBin2DFunc:
binning using weighted quantiles
"""

import feather
import numpy as np
import pandas as pd

def SimpleBinFunc(data, para, bins, Bounds='oc',
    saveRes=False, dirRes=None, saveLog=False, pathLog=None):
    """
    Function for binning (simple case)
    """
    print("Start binning (SimpleBinFunc)...")

    if saveLog:
        log = open(pathLog, "w")
        print("bin,number", file=log)

    data_final = []
    for i in range(len(bins)-1):
        vmin = bins[i]
        vmax = bins[i+1]

        if Bounds == 'oc':
            data = data[(data[para] > vmin) & (data[para] <= vmax)]
            N = len(data)
            print(f"Number in bin ({vmin}, {vmax}]: {N}")
        elif Bounds == 'co':
            data = data[(data[para] >= vmin) & (data[para] < vmax)]
            N = len(data)
            print(f"Number in bin [{vmin}, {vmax}): {N}")

        if saveRes:
            outpath = dirRes + str(vmin) + '_' + str(vmax) + '.feather'
        
            data = data.reset_index(drop=True)
            data.to_feather(outpath)
            print("Results saved to", outpath)
        else:
            data_final.append(data)

        if saveLog:
            key = str(vmin) + '_' + str(vmax)
            print(key, N, sep=',', file=log)

    print("Finished binning (SimpleBinFunc).")

    if ~saveRes:
        return data_final

def WgQuantile1DFunc(values, weights, pq):
    """
    Function for the weighted quantile of a 1D numpy array.

    """
    print("Start computing the weighted quantile (WgQuantile1DFunc)...")

    # Sort the data
    ind_sorted = np.argsort(values)
    v_sorted = values[ind_sorted]
    wg_sorted = weights[ind_sorted]
    
    # Compute the auxiliary arrays
    Sn = np.cumsum(wg_sorted)
    Pn = (Sn-0.5*wg_sorted)/np.sum(wg_sorted)

    # Get the quantiles
    res = np.interp(pq, Pn, v_sorted)
    
    print("Finished computing the weighted quantile (WgQuantile1DFunc).")
    return res


def WgBin2DFunc(v1, v2, wgs, Nbin1, Nbin2):
    """
    Function for binning using weighted quantiles
    """
    print("Start binning (WgBin2DFunc)...")

    # Define the probabilities for the quantiles based on the number of bins
    pq1 = np.linspace(0,1.0,Nbin1+1)
    pq2 = np.linspace(0,1.0,Nbin2+1)

    # Calculate quantiles for v1
    q1 = WgQuantile1DFunc(v1, wgs, pq1)
        
    #Compute quantiles for v2 in each v1 bin 
    q1s=[]
    q2s=[]
    for i in range(len(q1)-1):
        mask = (v1>=q1[i])&(v1<q1[i+1])
        q2 = WgQuantile1DFunc(v2[mask], wgs[mask], pq2)

        # q1s.append(q1)
        q2s.append(q2)

    print("Finished binning (WgQuantile1DFunc).")
    # return np.array(q1s), np.array(q2s)
    return q1, np.array(q2s)


if __name__=='__main__':

    import time

    # +++++++++++++ WgBin2DFunc
    Start = time.time()

    # data
    inpath = "/disks/shear15/ssli/SimCat/old/Selec__ZB9_in__13.feather"
    data = feather.read_dataframe(inpath)

    # Number of bins
    Nbin1 = 20
    Nbin2 = 20

    # binning values
    # SNR
    v1 = data['snr_model'].values

    # measured ellipticity
    e1_meas = data['e1'].values
    e2_meas = data['e2'].values
    mod_e_meas = np.sqrt(e1_meas**2.+e2_meas**2.)
    q_meas = (1.-mod_e_meas)/(1.+mod_e_meas)
    # r_ab
    size_ab = (data['size_out'].values)*np.sqrt(q_meas)
    # r_psf
    size_psf = data['psf_size_in'].values
    # R
    v2 = size_psf**2./(size_ab**2.+size_psf**2.)

    # weights
    wg = data['LFweight'].values

    # calculation
    v1_bounds, v2_bounds = WgBin2DFunc(v1, v2, wg, Nbin1, Nbin2)
    # print(len(v1_bounds), len(v2_bounds))
    # print(f"v1_bounds: {v1_bounds}")
    # print(f"v2_bounds: {v2_bounds}")
    # print(f"Running time: {time.time()-Start}")
    # Running time: 0.440138578414917

    # calculate the sum of lensfit weights in each bin
    # KV450
    inpath = "/disks/shear15/ssli/KV450/selected/G12__Z_B__13.feather"
    dat = feather.read_dataframe(inpath)
    # SNR
    v1_KV450 = dat['model_SNratio'].values + 0.000001

    # measured ellipticity ????????????????
    e1_meas_KV450 = dat['bias_corrected_e1'].values
    e2_meas_KV450 = dat['bias_corrected_e2'].values
    mod_e_meas_KV450 = np.sqrt(e1_meas_KV450**2.+e2_meas_KV450**2.)
    q_meas_KV450 = (1.-mod_e_meas_KV450)/(1.+mod_e_meas_KV450)
    # r_ab
    size_ab_KV450 = (dat['bias_corrected_scalelength_pixels'].values)*np.sqrt(q_meas_KV450)
    # r_psf
    psf_Q11 = dat['PSF_Q11'].values
    psf_Q12 = dat['PSF_Q12'].values
    psf_Q22 = dat['PSF_Q22'].values
    size_psf_KV450 = (psf_Q11*psf_Q22-psf_Q12**2)**0.5
    # R
    v2_KV450 = size_psf_KV450**2./(size_ab_KV450**2.+size_psf_KV450**2.)

    # weights
    wg_KV450 = dat['recal_weight'].values



    for i in range(Nbin1):
        mask1 = (v1>=v1_bounds[i])&(v1<v1_bounds[i+1])
        mask1_KV450 = (v1_KV450>=v1_bounds[i])&(v1_KV450<v1_bounds[i+1])
        v2_bound = v2_bounds[i]
        for j in range(Nbin2):
            mask2 = (v2>=v2_bound[j])&(v2<v2_bound[j+1])
            mask = mask1 & mask2

            wg_bin = np.sum(wg[mask])

            mask2_KV450 = (v2_KV450>=v2_bound[j])&(v2_KV450<v2_bound[j+1])
            mask_KV450 = mask1_KV450 & mask2_KV450

            wg_bin_KV450 = np.sum(wg_KV450[mask_KV450])

            print(wg_bin/wg_bin_KV450)






