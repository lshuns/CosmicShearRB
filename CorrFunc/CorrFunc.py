#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 16:09:28 2019

@author: ssli

Module to calculate the data vector (correlation function)

MeanFunc:
    weighted average of e for a given data set (for c-term determination)

CorrFunc:
    Using TreeCorr to calculate the 2-point shear correlation function
    e is modified using weighted e 
    bin_slop is set to be 0.05

CorrPlotFunc:
    Rearrange the results to be used in plot.

CorrCosmoFunc:
    Rearrange the results to be used in cosmo analysis.

"""


import feather
import treecorr
import pandas as pd
import numpy as np

def MeanFunc(e1, e2, wt):
    """
    Function for weighted average calculation
    """

    # boots
    nboot = 30
    e1w = np.zeros(nboot)
    e2w = np.zeros(nboot)

    ngals = len(e1)
    # weighted mean in each boot
    for i in range(nboot):
        idx = np.random.randint(0, ngals, ngals)
        sow = np.sum(wt[idx])
        e1w[i] = np.dot(e1[idx], wt[idx]) / sow
        e2w[i] = np.dot(e2[idx], wt[idx]) / sow
    e1_ave_b = np.mean(e1w)
    e1_err_b = np.std(e1w)
    e2_ave_b = np.mean(e2w)
    e2_err_b = np.std(e2w)

    # weight mean from direct calculation
    e1_ave = np.dot(e1, wt) / np.sum(wt)
    e2_ave = np.dot(e2, wt) / np.sum(wt)

    # output
    out = {'e1_ave': e1_ave, 'e2_ave': e2_ave, 
        'e1_ave_b': e1_ave_b, 'e1_err_b': e1_err_b, 
        'e2_ave_b': e2_ave_b, 'e2_err_b': e2_err_b}

    return out


def CorrFunc(cat1, cat2, 
                outpath, nbins, mins, maxs, units, bin_slop, nthr):
    """
    Function for correlation calculation (auto & cross)
    """

    gg = treecorr.GGCorrelation(nbins=nbins, min_sep=mins, max_sep=maxs, 
                                sep_units=units, bin_slop=bin_slop)

    if cat2 == None:
        # auto-correlation
        gg.process(cat1, num_threads=nthr)
    else:
        # cross-correlation
        gg.process(cat1, cat2, num_threads=nthr)
    gg.write(outpath)
    print("TreeCorr results saved in", outpath)


def CorrPlotFunc(Nbins=None, 
                        inpathF=None, inpathP=None, 
                        m_list=[],
                        outpath=None):
    """
    Rearrange the results to be used in plot
    """

    # Initialise all output colunms
    thetas = []
    xi_pm = []
    idx_pm = []
    idx_tomo_z1 = []
    idx_tomo_z2 = []

    for idx_z1 in range(Nbins):
        for idx_z2 in range(idx_z1, Nbins):

            inpath = inpathF + 'tomo_' + str(idx_z1+1) + '_' + str(idx_z2+1) + inpathP
            dat = np.loadtxt(inpath)

            # theta
            theta_bins = np.concatenate((dat[0:7,1],dat[-6:,1]))
            thetas.append(theta_bins)
            # xi_pm
            xi_pm_bins = np.concatenate((dat[0:7,3], dat[-6:,4]))
            if m_list != []:
                Kplus1 = (1. + m_list[idx_z1]) * (1 + m_list[idx_z2]) 
        
                xi_pm_bins = xi_pm_bins / Kplus1

            xi_pm.append(xi_pm_bins)
            # idx_pm
            idx_pm_bins = np.concatenate((np.ones(7), np.ones(6) + 1))
            idx_pm.append(idx_pm_bins)
            # idx_tomo
            idx_tomo_z1.append(np.ones(7+6)*(idx_z1+1))
            idx_tomo_z2.append(np.ones(7+6)*(idx_z2+1))

    thetas = np.concatenate(thetas)
    xi_pm = np.concatenate(xi_pm)
    idx_pm = np.concatenate(idx_pm)
    idx_tomo_z1 = np.concatenate(idx_tomo_z1)
    idx_tomo_z2 = np.concatenate(idx_tomo_z2)

    idx_run = np.arange(len(idx_pm)) + 1
    savedata = np.column_stack((idx_run, thetas, xi_pm, idx_pm, idx_tomo_z1, idx_tomo_z2))
    header = ' i    theta(i)\'        xi_p/m(i)  (p=1, m=2)  itomo   jtomo'
    np.savetxt(outpath, savedata, header=header, delimiter='\t', fmt=['%4i', '%.5e', '%12.5e', '%i', '%i', '%i'])
    

    print("TreeCorr Results rearranged for plot saved to", outpath)


def CorrCosmoFunc(Nbins=None, 
                        inpathF=None, inpathP=None, 
                        m_list=[],
                        outpathF=None, outpathP=None):
    """
    Rearrange the results to be used in cosmo analysis
    """

    header = ' theta(arcmin)   xi_p    xi_m'

    for idx_z1 in range(Nbins):
        for idx_z2 in range(idx_z1, Nbins):

            inpath = inpathF + 'tomo_' + str(idx_z1+1) + '_' + str(idx_z2+1) + inpathP
            dat = np.loadtxt(inpath)

            # theta
            theta = dat[:,1]
            # xi_p
            xi_p = dat[:,3]
            xi_m = dat[:,4]

            if m_list != []:                
                Kplus1 = (1. + m_list[idx_z1]) * (1 + m_list[idx_z2]) 
        
                xi_p = xi_p / Kplus1
                xi_m = xi_m / Kplus1

                outpath = outpathF + 'tomo_' + str(idx_z1+1) + '_' + str(idx_z2+1) + '_withK' + outpathP
            else:
                outpath = outpathF + 'tomo_' + str(idx_z1+1) + '_' + str(idx_z2+1) + '_withoutK' + outpathP

            savedata = np.column_stack((theta, xi_p, xi_m))
            np.savetxt(outpath, savedata, header=header, delimiter='\t', fmt=['%.5e', '%.5e', '%.5e'])
    
            print("TreeCorr Results rearranged for cosmo analysis saved to", outpath)