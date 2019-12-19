#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 13:28:44 2019

@author: ssli

MCFunc:
    Function for fitting m and c from data


"""

import feather
import pandas as pd
import numpy as np
from scipy import stats, optimize

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
    inpath_Arun = "/disks/shear15/ssli/SimCat/SimCatSelec_Arun_all_13_PSF.feather"
    inpath_new = "/disks/shear15/ssli/SimCat/SimCatSelec_all_13_PSF.feather"
    inpath_low = "/disks/shear15/ssli/SimCat/split/SimCatSelec_all_13_PSF_TB9_in_3_less.feather"
    inpath_high = "/disks/shear15/ssli/SimCat/split/SimCatSelec_all_13_PSF_TB9_in_3_greater.feather"
    inpaths = [inpath_Arun, inpath_new, inpath_high, inpath_low]

    # linregress
    METHOD = 'linregress'
    Start = time.time()

    for path in inpaths:
        data = feather.read_dataframe(path)
        g1 = data['g1'].values
        g2 = data['g2'].values

        wg = data['LFweight'].values
        e1 = data['e1'].values
        e2 = data['e2'].values

        print(path)
    
        MCFunc(g1, g2, e1, e2, wg, METHOD)

    print("Running time", time.time()-Start)
    

# /disks/shear15/ssli/SimCat/SimCatSelec_Arun_all_13_PSF.feather
# m1 = -0.008226157898461817 +- 0.0011035142956440916
# m2 = -0.005615895528964132 +- 0.0011894159333708442

# /disks/shear15/ssli/SimCat/SimCatSelec_all_13_PSF.feather
# m1 = -0.006390492795793401 +- 0.0009419661670895168
# m2 = -0.004149017271736111 +- 0.0012030771371226198

# /disks/shear15/ssli/SimCat/split/SimCatSelec_all_13_PSF_TB9_in_3_greater.feather
# m1 = 0.003443817090559742 +- 0.0011026532287912445
# m2 = 0.009581121655187985 +- 0.0013579449791004464

# /disks/shear15/ssli/SimCat/split/SimCatSelec_all_13_PSF_TB9_in_3_less.feather
# m1 = -0.01973750966868082 +- 0.001940890385471888
# m2 = -0.02278890408800005 +- 0.0020136098135490133
