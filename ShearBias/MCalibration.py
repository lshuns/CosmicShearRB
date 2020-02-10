#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 13:24:43 2020

@author: ssli

Module to calculate the m bias 

mcFitFunc:
    Shear bias function.

WgQuantile1DFunc:
    Calculate the weighted quantile by given probabilities 
        designed for 1D numpy array.

WgBin2DFunc:
    Calculate the weighted quantile by given bin numbers
        designed for 2D numpy array

mcCalFunc:    
    Calculating the residual shear bias in 2-d bins
"""

import numpy as np
from scipy import optimize
import pandas as pd
from astropy.io import fits

## All possible g1,g2 combinations
g1Range = np.array([-0.04,0.00,0.04,0.00,-0.0283,+0.0283,+0.0283,-0.0283])
g2Range = np.array([0.00,0.04,0.00,-0.04,+0.0283,+0.0283,-0.0283,-0.0283])


def mcFitFunc(x, m, c):
    """
    Shear bias function.
    """
    return (1.+m)*x+c


def WgQuantile1DFunc(values, weights, pq):
    """
    Calculate the weighted quantile by given probabilities 
        designed for 1D numpy array.


    Parameters
    ----------


    """

    # Sort the data
    ind_sorted = np.argsort(values)
    v_sorted = values[ind_sorted]
    wg_sorted = weights[ind_sorted]
    
    # Compute the auxiliary arrays
    Sn = np.cumsum(wg_sorted)
    Pn = (Sn-0.5*wg_sorted)/np.sum(wg_sorted)

    # Get the quantiles
    res = np.interp(pq, Pn, v_sorted)
    
    return res


def WgBin2DFunc(v1, v2, wgs, Nbin1, Nbin2):
    """
    Calculate the weighted quantile by given bin numbers
        designed for 2D numpy array
    """

    # Define the probabilities for the quantiles based on the number of bins
    pq1 = np.linspace(0,1.0,Nbin1+1)
    pq2 = np.linspace(0,1.0,Nbin2+1)

    # Calculate quantiles for v1
    q1 = WgQuantile1DFunc(v1, wgs, pq1)
        
    #Compute quantiles for v2 in each v1 bin
    q2s=[]
    for i in range(len(q1)-1):
        mask = (v1>=q1[i])&(v1<q1[i+1])
        q2 = WgQuantile1DFunc(v2[mask], wgs[mask], pq2)
        q2s.append(q2)

    return q1, np.array(q2s)


def mcCalFunc(dataSim, dataReal, 
                Nbin1, Nbin2, pq):    
    """
    Calculating the residual shear bias in 2-d bins
    """


    # helper quantities
    # Simulation
    snrSim = dataSim['snr_model'].values
    #
    g1_inSim = dataSim['g1'].values
    g2_inSim = dataSim['g2'].values
    #
    e1Sim = dataSim['e1'].values
    e2Sim = dataSim['e2'].values
    eSim = np.sqrt(e1Sim**2 + e2Sim**2)
    #
    size_out_altSim = dataSim['size_out'].values*np.sqrt((1.-eSim)/(1.+eSim))
    #
    RSim = dataSim['psf_size_in'].values/(size_out_altSim**2+dataSim['psf_size_in'].values)
    #
    wgSim= dataSim['LFweight'].values
    #
    # Data
    snrReal = dataReal['model_SNratio'].values
    #Define PSF size
    size_psfReal = np.sqrt(dataReal['PSF_Q11'].values*dataReal['PSF_Q22'].values - dataReal['PSF_Q12'].values**2)
    #Define |e| for the 3 blindings
    eReal = np.sqrt(dataReal['bias_corrected_e1'].values**2 + dataReal['bias_corrected_e2'].values**2)
    #Define circularised galaxy size
    size_abReal = dataReal['bias_corrected_scalelength_pixels'].values*np.sqrt((1.-eReal)/(1.+eReal))
    #Define galaxy 'resolution'
    RReal = size_psfReal/(size_abReal**2+size_psfReal)
    # weight
    wgReal = dataReal['recal_weight'].values

    # 2D binning
    #Calculate the bins such that each bin contains the same number of points.
    bin1_bounds, bin2_bounds = WgBin2DFunc(snrSim, RSim, wgSim, Nbin1, Nbin2)

    wgRealSums = []
    wgReal2Sums = []
    m1s = []
    m2s = []
    m1_errs = []
    m2_errs = []
    m1_err_BSs = []
    m2_err_BSs = []
    m_err_BSs = []
    for i in range(Nbin1):
        lower1 = bin1_bounds[i]
        upper1 = bin1_bounds[i+1]
        #
        mask1Sim = (snrSim>=lower1)&(snrSim<upper1)
        mask1Real = (snrReal>=lower1)&(snrReal<upper1)

        for j in range(Nbin2):
            lower2 = bin2_bounds[i][j]
            upper2 = bin2_bounds[i][j+1]
            #
            mask2Sim = (RSim>=lower2)&(RSim<upper2)
            mask2Real = (RReal>=lower2)&(RReal<upper2)
            #
            maskSim = mask1Sim & mask2Sim
            maskReal = mask1Real & mask2Real

            # mask input parameters
            # Simulation
            wgSim_mask = wgSim[maskSim]
            #
            e1Sim_mask = e1Sim[maskSim]
            e2Sim_mask = e2Sim[maskSim]
            #
            g1_inSim_mask = g1_inSim[maskSim]
            g2_inSim_mask = g2_inSim[maskSim]
            # data
            wgReal_mask = wgReal[maskReal]
            wgRealSums.append(np.sum(wgReal_mask))

            # prepare shear parameters for mc fitting
            g1_out=[]
            g2_out=[]
            g_out_w=[]
            g1_in_used=[]
            g2_in_used=[]
            for kk in range(len(g1Range)):
                maskShear=(g1_inSim_mask==g1Range[kk])&(g2_inSim_mask==g2Range[kk])
                numMasked=len(e1Sim_mask[maskShear])
                if (numMasked >0):
                    #Calculating bin average for calibration quantities
                    g1_out.append(np.average(e1Sim_mask[maskShear], weights=wgSim_mask[maskShear]))
                    g2_out.append(np.average(e2Sim_mask[maskShear], weights=wgSim_mask[maskShear]))
                    g_out_w.append(1./(np.sum(wgSim_mask[maskShear]))**0.5)
                    #
                    g1_in_used.append(g1Range[kk])
                    g2_in_used.append(g2Range[kk])

            # Start mc fitting
            numShear=len(g1_out)
            if(numShear<3):
                print('Cannot do the regression in bin ', \
                    bin1_bounds[i], bin1_bounds[i+1], bin2_bounds[i][j], bin2_bounds[i][j+1], \
                    ' less than 3 shear values! (', numShear, ')')
                exit()
            else:
                g1_in_used = np.array(g1_in_used)
                g2_in_used = np.array(g2_in_used)
                g1_out = np.array(g1_out)
                g2_out = np.array(g2_out)
                g_out_w = np.array(g_out_w)

                m1c1, err1 = optimize.curve_fit(mcFitFunc, xdata=g1_in_used, ydata=g1_out, sigma=g_out_w)
                m2c2, err2 = optimize.curve_fit(mcFitFunc, xdata=g2_in_used, ydata=g2_out, sigma=g_out_w)

                m1 = m1c1[0]
                m1_err = (err1[0,0])**0.5
                # # #
                # c1 = m1c1[1]
                # c1_err = (err1[1,1])**0.5
                #
                m2 = m2c2[0]
                m2_err = (err2[0,0])**0.5
                # # #
                # c2 = m2c2[1]
                # c2_err =(err2[1,1])**0.5
                #
                #
                # m = (m1 + m2)/2.

                # Performing Bootstrap
                nboot = 30
                m1_sample = np.zeros(nboot)
                m2_sample = np.zeros(nboot)
                m_sample = np.zeros(nboot)
                # c1_sample = np.zeros(nboot)
                # c2_sample = np.zeros(nboot)
                for BS_index in range(nboot):
                    # Retrieving random shears
                    index = np.random.randint(0,numShear,numShear)

                    BS_g1_in = g1_in_used[index]
                    BS_g2_in = g2_in_used[index] 
                    BS_g1_out = g1_out[index]
                    BS_g2_out = g2_out[index] 
                    BS_g_out_w = g_out_w[index]

                    m1c1, err1 = optimize.curve_fit(mcFitFunc, xdata=BS_g1_in, ydata=BS_g1_out, sigma=BS_g_out_w)
                    m2c2, err2 = optimize.curve_fit(mcFitFunc, xdata=BS_g2_in, ydata=BS_g2_out, sigma=BS_g_out_w)

                    m1_sample[BS_index] = m1c1[0]
                    m2_sample[BS_index] = m2c2[0]
                    m_sample[BS_index] = (m1c1[0]+m2c2[0])/2.
                    # c1_sample[BS_index] = m1c1[1]
                    # c2_sample[BS_index] = m2c2[1]
                        
                m1_err_BS = np.std(m1_sample)
                m2_err_BS = np.std(m2_sample)
                m_err_BS = np.std(m_sample)
                # c1_err_BS = np.std(c1_sample)
                # c2_err_BS = np.std(c2_sample)
                #
                m1s.append(m1)
                m2s.append(m2)
                m1_errs.append(m1_err)
                m2_errs.append(m2_err)
                m1_err_BSs.append(m1_err_BS)
                m2_err_BSs.append(m2_err_BS)
                m_err_BSs.append(m_err_BS)


    wgRealSums = np.array(wgRealSums)
    m1s = np.array(m1s)
    m2s = np.array(m2s)
    m1_errs = np.array(m1_errs)
    m2_errs = np.array(m2_errs)
    m1_err_BSs = np.array(m1_err_BSs)
    m2_err_BSs = np.array(m2_err_BSs)
    m_err_BSs = np.array(m_err_BSs)


    m1_final = np.average(m1s, weights=wgRealSums)
    m2_final = np.average(m2s, weights=wgRealSums)
    m_final = (m1_final+m2_final)/2.

    m1_err_final = np.sqrt(np.sum((wgRealSums*m1_errs)**2.))/np.sum(wgRealSums)
    m2_err_final = np.sqrt(np.sum((wgRealSums*m2_errs)**2.))/np.sum(wgRealSums)

    m1_err_BS_final = np.sqrt(np.sum((wgRealSums*m1_err_BSs)**2.))/np.sum(wgRealSums)
    m2_err_BS_final = np.sqrt(np.sum((wgRealSums*m2_err_BSs)**2.))/np.sum(wgRealSums)
    m_err_BS_final = np.sqrt(np.sum((wgRealSums*m_err_BSs)**2.))/np.sum(wgRealSums)

    #
    res = {"m_final": m_final, 'm_err_BS_final': m_err_BS_final, \
        'm1_final': m1_final, 'm2_final': m2_final, \
        'm1_err_final': m1_err_final, 'm2_err_final': m2_err_final, \
        'm1_err_BS_final': m1_err_BS_final, 'm2_err_BS_final': m2_err_BS_final}
    pq.put(res)
