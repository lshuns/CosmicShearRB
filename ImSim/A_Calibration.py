#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 13:28:44 2019

@author: ssli

MCFunc:
    Function for fitting m and c from data


Bug:
    Error from METHOD = 'self' (MCFunc) is wrong
"""

import feather
import pandas as pd
import numpy as np
from scipy import stats, optimize

import os
import sys
sys.path.insert(0,os.path.realpath('..')) 
from utility import LinearFitFunc

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
    print("Start fitting m and c (MCFunc)...")

    # All possible g1,g2 combinations
    g1Range = np.array([-0.04,0.00,0.04,0.00,-0.0283,+0.0283,+0.0283,-0.0283])
    g2Range = np.array([0.00,0.04,0.00,-0.04,+0.0283,+0.0283,-0.0283,-0.0283])

    g1_meas, g2_meas, wt_g = [], [], []
    for i in range(len(g1Range)):
        cuts = ((g1==g1Range[i])&(g2==g2Range[i]))

        print(f"Object Number in {i} shear combination:", cuts.sum())

        g1hat = np.average(e1[cuts], weights=wg[cuts])
        g2hat = np.average(e2[cuts], weights=wg[cuts])

        g1_meas.append(g1hat)
        g2_meas.append(g2hat)

    g1_meas = np.array(g1_meas)
    g2_meas = np.array(g2_meas)

    if METHOD == 'linregress':
        # Using linregress\
        print("Using linregress from scipy.stats")
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

    elif METHOD == 'self':
        # Using self-defined LinearFitFunc
        print("Using self-defined LinearFitFunc from utility/fitting.py")
        a, b, sa, sb, rchi2, dof = LinearFitFunc(g1Range, g1_meas)
        a2, b2, sa2, sb2, rchi22, dof2 = LinearFitFunc(g2Range, g2_meas)

        m1 = a - 1.
        sm1 = sa
        m2 = a2 - 1.
        sm2 = sa2
        #
        c1 = b
        sc1 = sb
        c2 = b2
        sc2 = sb2

    elif METHOD == 'curve_fit':
        # Using curve_fit
        print("Using curve_fit from scipy.optimize")
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

    print(f"m1 = {m1} +- {sm1}")
    print(f"m2 = {m2} +- {sm2}")
    print(f"c1 = {c1} +- {sc1}")
    print(f"c2 = {c2} +- {sc2}")

    print("Finished fitting m and c (MCFunc).")

    return m1, sm1, m2, sm2, c1, sc1, c2, sc2


if __name__ == "__main__":

    import time

    # ++++++++++++++++++++++++++++++ MCFunc
    # data
    # inpath = "/disks/shear15/ssli/SimCat/SimCatOrigi.feather"
    # inpath = "/disks/shear15/ssli/SimCat/SimCatSelec_Arun_all_13_PSF.feather"
    inpath = "/disks/shear15/ssli/SimCat/SimCatSelec_all_13_PSF.feather"


    data = feather.read_dataframe(inpath)
    g1 = data['g1'].values
    g2 = data['g2'].values

    wg = data['LFweight'].values
    e1 = data['e1'].values
    e2 = data['e2'].values

    # linregress
    METHOD = 'linregress'
    Start = time.time()
    
    MCFunc(g1, g2, e1, e2, wg, METHOD)
    print("Running time for linregress:", time.time()-Start)
    # Running time for linregress: 0.8835294246673584

# m1 = -0.008226157898461817 +- 0.0011035142956440916
# m2 = -0.005615895528964132 +- 0.0011894159333708442
# c1 = 0.00010176980867981911 +- None
# c2 = 0.0007877680473029613 +- None


    # self
    METHOD = 'self'
    Start = time.time()
    
    MCFunc(g1, g2, e1, e2, wg, METHOD)
    print("Running time for self:", time.time()-Start)
    # Running time for self: 0.886962890625

# m1 = -0.008226157898461817 +- 12.496524887206434
# m2 = -0.005615895528964132 +- 12.496524887206434
# c1 = 0.00010176980867981912 +- 0.3535533905932738
# c2 = 0.0007877680473029613 +- 0.3535533905932738


    # curve_fit
    METHOD = 'curve_fit'
    Start = time.time()
    
    MCFunc(g1, g2, e1, e2, wg, METHOD)
    print("Running time for curve_fit:", time.time()-Start)
    # Running time for curve_fit: 0.8620116710662842

# m1 = -0.008226157956150568 +- 0.001103514332243542
# m2 = -0.005615895240730806 +- 0.001189417034808952
# c1 = 0.00010176959425471485 +- 3.1220791879571026e-05
# c2 = 0.0007877682913796833 +- 3.365112068114336e-05
