#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:00:57 2020

@author: ssli

Function to calculate n_density, n_eff & sigma_e, 
    which are required by covariance matrix calculation
"""

import numpy as np
import math


def NeffSigmaeFunc(id_zbin, e1, e2, wg, area, outpath, WorA='w'):
    """
    Calculate n_eff & sigma_e for covariance matrix calculation purpose
    """
    log = open(outpath, WorA)
    if WorA == 'w':
        print("# ========================================================================", file=log)
        print("# 'n_density' is the 'raw' galaxy number density without weights", file=log)
        print("# 'n_eff' is the effective number density defined by Heymans et al. (2012)", file=log)
        print("# two types of 'sigma_e*'' are calculated, ", file=log)
        print("#        those with '_wsq' are calculated with the square weights, while the other two are using the single weights", file=log)
        print("# ========================================================================", file=log)
        print("# zbin n_obj n_density[1/arcmin^2] n_eff[1/arcmin^2] sigma_e1 sigma_e2 sigma_e1_wsq sigma_e2_wsq", file=log)
    
    number = len(wg)

    n_density = number/area
    
    n_eff = np.sum(wg)**2/np.sum(wg**2)/area
    sigma_e1 = weighted_avg_and_std(e1, wg)[1]
    sigma_e2 = weighted_avg_and_std(e2, wg)[1]

    sigma_e1_wsq = weighted_avg_and_std(e1, wg**2)[1]
    sigma_e2_wsq = weighted_avg_and_std(e2, wg**2)[1]

    print(id_zbin, number, n_density, n_eff, sigma_e1, sigma_e2, sigma_e1_wsq, sigma_e2_wsq, sep=' ', file=log)
    log.close()


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))

if __name__ == '__main__':
    import feather
    import time

    # # +++++++++++++++++++++++++++ test
    # start = time.time()
    # # input directory
    # indir = "/disks/shear15/ssli/KV450/tomo/all_tomo"

    # area = 341.3 * 3600. # 1/arcmin^2
    # outpath = './test/whole.txt'

    # WorA = 'w'
    # for i in range(5):
    #     inpath = indir + str(i+1) + '.feather'

    #     indata = feather.read_dataframe(inpath)

    #     id_zbin = i + 1
    #     e1 = indata['bias_corrected_e1']
    #     e2 = indata['bias_corrected_e2']
    #     wg = indata['recal_weight']

    #     NeffSigmaeFunc(id_zbin, e1, e2, wg, area, outpath, WorA)
    #     WorA = 'a'
    #     print("Finished in", id_zbin)

    # print("all done in", time.time()-start)
    # # all done in 37.25660729408264

    # # +++++++++++++++++++++++++++ main running
    # input directory
    indir = "/disks/shear15/ssli/KV450/split/all_tomo"
    # input postfix
    inP_r = "_T_B_less3"
    inP_b = "_T_B_greater3"
    inPs = [inP_r, inP_b]


    area = 341.3 * 3600. # 1/arcmin^2

    outdir = "/disks/shear15/ssli/CosmicShear/covariance/prepare/Ndensity_sigmae"
    # output postfix
    outPs = ["_red", "_blue"]

    for k in range(len(inPs)):
        WorA = 'w' 
        for i in range(5):
            inpath = indir + str(i+1) + inPs[k] + '.feather'
            indata = feather.read_dataframe(inpath)
            outpath = outdir + outPs[k] + '.txt'

            id_zbin = i + 1
            e1 = indata['bias_corrected_e1']
            e2 = indata['bias_corrected_e2']
            wg = indata['recal_weight']

            NeffSigmaeFunc(id_zbin, e1, e2, wg, area, outpath, WorA)
            WorA = 'a'
            print("Finished in", id_zbin, inPs[k])

