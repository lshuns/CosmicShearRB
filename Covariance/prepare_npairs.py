#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 19:46:22 2020

@author: ssli

Get the number of pairs in correlation function calculation
    required by covariance calculation (for shot noise)
"""

import numpy as np

inpathF = "/disks/shear15/ssli/CosmicShear/data_vector/cross_subsamples/xi_blue_red_tomo_"
outpath = "/disks/shear15/ssli/CosmicShear/data_vector/cross_subsamples/npairs_blue_red_all.dat"

nzbin = 10
ntheta = 9

theta_list  = []
npairs_list = []
for i in range(10):
    for j in range(i, 10):
        inpath = inpathF + str(i+1) + '_' + str(j+1) + '.dat'
        print(str(i+1), str(j+1))
        data = np.loadtxt(inpath)

        theta_list.append(data[:,1])
        # twice the unique pairs in auto-correlations
        if i==j:
            npairs_list.append(data[:,10]*2)
        else:
            npairs_list.append(data[:,10])

# average theta to smear the slight difference
theta = np.average(theta_list, axis=0)

# combine npairs with theta
data_final = [theta]
for npairs in npairs_list:
    data_final.append(npairs)

data_final = np.array(data_final).transpose()

header = "theta_mean    1-1  1-2 ... 1-10  2-2  2-3 ... 2-10  3-3 ... 10-10"
np.savetxt(outpath, data_final, fmt='%.4e', header=header)