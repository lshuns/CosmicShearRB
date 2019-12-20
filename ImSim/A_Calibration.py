#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 13:28:44 2019

@author: ssli

MCFunc:
    fitting m and c from data
MCalFunc:
    calculate the residual shear bias (m only)
    Method C of Fenech Conti et al. 2017 
"""

import feather
import pandas as pd
import numpy as np
from scipy import stats, optimize

import os
import sys
sys.path.insert(0,os.path.realpath('..')) 
from utility import WgBin2DFunc


def mcFitFunc(x, m, c):
    """
    Function required for fitting m and c when METHOD=='curve_fit'
    """
    return (1.+m)*x+c


def MCFunc(g1, g2, e1, e2, wg,
            METHOD='linregress'):
    """
    Function for fitting m and c from data
    """
    # print("Start fitting m and c (MCFunc)...")

    # All possible g1,g2 combinations
    g1Range = np.array([-0.04,0.00,0.04,0.00,-0.0283,+0.0283,+0.0283,-0.0283])
    g2Range = np.array([0.00,0.04,0.00,-0.04,+0.0283,+0.0283,-0.0283,-0.0283])

    g1_meas, g2_meas, wt_g = [], [], []
    for i in range(len(g1Range)):
        cuts = ((g1==g1Range[i])&(g2==g2Range[i]))

        # print(f"Object Number in {i} shear combination:", cuts.sum())

        g1hat = np.average(e1[cuts], weights=wg[cuts])
        g2hat = np.average(e2[cuts], weights=wg[cuts])

        g1_meas.append(g1hat)
        g2_meas.append(g2hat)

    g1_meas = np.array(g1_meas)
    g2_meas = np.array(g2_meas)

    if METHOD == 'linregress':
        # Using linregress
        # print("Using linregress from scipy.stats")
        slope, intercept, r_value, p_value, stderr = stats.linregress(g1Range, g1_meas)
        slope2, intercept2, r_value2, p_value2, stderr2 = stats.linregress(g2Range, g2_meas)
        #
        m1 = slope - 1.
        sm1 = stderr
        m2 = slope2 - 1. 
        sm2 = stderr2
        #
        c1 = intercept
        sc1 = None
        c2 = intercept2
        sc2 = None

    elif METHOD == 'curve_fit':
        # Using curve_fit
        # print("Using curve_fit from scipy.optimize")
        m1c1, err1 = optimize.curve_fit(mcFitFunc, xdata=g1Range, ydata=g1_meas)
        m2c2, err2 = optimize.curve_fit(mcFitFunc, xdata=g2Range, ydata=g2_meas)

        m1 = m1c1[0]
        sm1 = (err1[0,0])**0.5
        m2 = m2c2[0]
        sm2 = (err2[0,0])**0.5
        #
        c1 = m1c1[1]
        sc1 = (err1[1,1])**0.5
        c2 = m2c2[1]
        sc2 = (err2[1,1])**0.5

    else:
        raise Exception("Unsupported fitting METHOD!")

    # print(f"m1 = {m1} +- {sm1}")
    # print(f"m2 = {m2} +- {sm2}")
    # print(f"c1 = {c1} +- {sc1}")
    # print(f"c2 = {c2} +- {sc2}")

    # print("Finished fitting m and c (MCFunc).")

    return m1, sm1, m2, sm2, c1, sc1, c2, sc2


def MCalFunc(data_sim, data_real):
    """
    Function for calculating the residual shear bias (m only)
    """

    # Number of bins
    Nbin_SNR = 20
    Nbin_R = 20

    # SNR
    SNR_sim = data_sim['snr_model'].values
    #
    SNR_real = data_real['model_SNratio'].values

    # measured ellipticity
    e1_sim = data_sim['e1'].values
    e2_sim = data_sim['e2'].values
    mod_e_sim = np.sqrt(e1_sim**2.+e2_sim**2.)
    q_sim = (1.-mod_e_sim)/(1.+mod_e_sim)
    #
    e1_real = data_real['bias_corrected_e1'].values
    e2_real = data_real['bias_corrected_e2'].values
    mod_e_real = np.sqrt(e1_real**2.+e2_real**2.)
    q_real = (1.-mod_e_real)/(1.+mod_e_real)

    # r_ab
    r_ab_sim = (data_sim['size_out'].values)*np.sqrt(q_sim)
    #
    r_ab_real = (data_real['bias_corrected_scalelength_pixels'].values)*np.sqrt(q_real)

    # r_psf
    size_psf_sim = data_sim['psf_size_in'].values
    #
    psf_Q11 = data_real['PSF_Q11'].values
    psf_Q12 = data_real['PSF_Q12'].values
    psf_Q22 = data_real['PSF_Q22'].values
    r2_psf_real = (psf_Q11*psf_Q22-psf_Q12**2.)**0.5

    # R
    R_sim = size_psf_sim/(r_ab_sim**2.+size_psf_sim)
    #
    R_real = r2_psf_real/(r_ab_real**2.+r2_psf_real)

    # weights
    wg_sim = data_sim['LFweight'].values
    #
    wg_real = data_real['recal_weight'].values

    # calculation bounds for binning
    SNR_bounds, R_bounds = WgBin2DFunc(SNR_sim, R_sim, wg_sim, Nbin_SNR, Nbin_R)

    # parameters for bias estimation
    g1 = data_sim['g1'].values
    g2 = data_sim['g2'].values
    e1 = data_sim['e1'].values
    e2 = data_sim['e2'].values

    # bias estimation for each bin
    wgs_real = []
    m1s = []
    m2s = []
    # sm1s = []
    # sm2s = []
    for i in range(Nbin_SNR):
        mask1_sim = (SNR_sim>=SNR_bounds[i])&(SNR_sim<SNR_bounds[i+1])
        mask1_real = (SNR_real>=SNR_bounds[i])&(SNR_real<SNR_bounds[i+1])

        R_bound = R_bounds[i]
        for j in range(Nbin_R):
            mask2_sim = (R_sim>=R_bound[j])&(R_sim<R_bound[j+1])
            mask2_real = (R_real>=R_bound[j])&(R_real<R_bound[j+1])

            mask_sim = mask1_sim & mask2_sim
            mask_real = mask1_real & mask2_real

            # weight difference
            tmp = np.sum(wg_real[mask_real])
            wgs_real.append(tmp)

            # bias
            m1, sm1, m2, sm2, c1, sc1, c2, sc2 = \
                MCFunc(g1[mask_sim], g2[mask_sim], e1[mask_sim], e2[mask_sim], wg_sim[mask_sim])
            m1s.append(m1)
            m2s.append(m2)
            # sm1s.append(sm1)
            # sm2s.append(sm2)

    m1f = np.average(m1s, weights=wgs_real)
    m2f = np.average(m2s, weights=wgs_real)

    mf = (m1f + m2f)/2.

    return m1f, m2f, mf


if __name__ == "__main__":

    import time

#     # ++++++++++++++++++++++++++++++ MCFunc
#     # data
#     # inpath = "/disks/shear15/ssli/SimCat/SimCatOrigi.feather"
#     inpath_Arun = "/disks/shear15/ssli/SimCat/SimCatSelec_Arun_all_13_PSF.feather"
#     inpath_new = "/disks/shear15/ssli/SimCat/SimCatSelec_all_13_PSF.feather"
#     inpath_low = "/disks/shear15/ssli/SimCat/split/SimCatSelec_all_13_PSF_TB9_in_3_less.feather"
#     inpath_high = "/disks/shear15/ssli/SimCat/split/SimCatSelec_all_13_PSF_TB9_in_3_greater.feather"
#     inpaths = [inpath_Arun, inpath_new, inpath_high, inpath_low]

#     # linregress
#     METHOD = 'linregress'
#     Start = time.time()

#     for path in inpaths:
#         data = feather.read_dataframe(path)
#         g1 = data['g1'].values
#         g2 = data['g2'].values

#         wg = data['LFweight'].values
#         e1 = data['e1'].values
#         e2 = data['e2'].values

#         print(path)
    
#         MCFunc(g1, g2, e1, e2, wg, METHOD)

#     print("Running time", time.time()-Start)
    

# # /disks/shear15/ssli/SimCat/SimCatSelec_Arun_all_13_PSF.feather
# # m1 = -0.008226157898461817 +- 0.0011035142956440916
# # m2 = -0.005615895528964132 +- 0.0011894159333708442

# # /disks/shear15/ssli/SimCat/SimCatSelec_all_13_PSF.feather
# # m1 = -0.006390492795793401 +- 0.0009419661670895168
# # m2 = -0.004149017271736111 +- 0.0012030771371226198

# # /disks/shear15/ssli/SimCat/split/SimCatSelec_all_13_PSF_TB9_in_3_greater.feather
# # m1 = 0.003443817090559742 +- 0.0011026532287912445
# # m2 = 0.009581121655187985 +- 0.0013579449791004464

# # /disks/shear15/ssli/SimCat/split/SimCatSelec_all_13_PSF_TB9_in_3_less.feather
# # m1 = -0.01973750966868082 +- 0.001940890385471888
# # m2 = -0.02278890408800005 +- 0.0020136098135490133


    # ++++++++++++++++++++++++++++++++ MCalFunc

    indir_sim = "/disks/shear15/ssli/SimCat/"
    inge_sim = "Bin_Arun_all_13_PSF__ZB9_in__"

    indir_real = "/disks/shear15/ssli/KV450/selected/"
    inge_real = "__Z_B__"
    patches = ["G9","G12","G15","G23","GS"]

    bins = ['13', '35', '57', '79', '912']


    for Bin in bins:
        inpath = indir_sim + inge_sim + Bin + '.feather'
        data_sim = feather.read_dataframe(inpath)

        data_tmp = []
        for patch in patches:
            inpath = indir_real + patch + inge_real + Bin + '.feather'
            data_tmp.append(feather.read_dataframe(inpath))
        data_real = pd.concat(data_tmp, ignore_index=True)

        m1f, m2f, mf = MCalFunc(data_sim, data_real)

        print(m1f, m2f, mf)
