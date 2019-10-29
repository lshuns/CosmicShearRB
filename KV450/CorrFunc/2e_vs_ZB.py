#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 16:21:10 2019

@author: ssli

weighted average of e for each tomographic bin
"""

import numpy as np
import pandas as pd

# input path
inpathF = "/disks/shear15/ssli/KV450/selected/mine/bins/zb_"
inpathL = ["13","35","57","79","912"]
inpathP = ".h5"

# redshift median
z_median = [0.2, 0.4, 0.6, 0.8, 1.0]

nboot = 30
e1w = np.zeros(nboot)
e2w = np.zeros(nboot)

# output path 
out = open("/disks/shear15/ssli/KV450/CS/mine/TreeCorr/full/mine/e_vs_ZB.dat", mode='w')
print("# z_median e1_ave e1_err e2_ave e2_err e1_ave_n e2_ave_n", file=out)

for j in range(len(inpathL)):
    print("In z_median =", z_median[j])

    inpath = inpathF + inpathL[j] + inpathP
    data = pd.read_hdf(inpath, key='whole', 
                        columns=["bias_corrected_e1", "bias_corrected_e2", "recal_weight"])

    e1 = data['bias_corrected_e1'].values
    e2 = data['bias_corrected_e2'].values
    wt = data['recal_weight'].values

    ngals = len(e1)
    # weighted mean in each boot
    for i in range(nboot):
    	idx = np.random.randint(0, ngals, ngals)
    	sow = np.sum(wt[idx])
    	e1w[i] = np.dot(e1[idx], wt[idx]) / sow
    	e2w[i] = np.dot(e2[idx], wt[idx]) / sow

    e1_ave = np.mean(e1w)
    e1_err = np.std(e1w)
    e2_ave = np.mean(e2w)
    e2_err = np.std(e2w)

    # weight mean from direct calculation
    e1_ave_n = np.dot(e1, wt) / np.sum(wt)
    e2_ave_n = np.dot(e2, wt) / np.sum(wt)

    print(z_median[j], e1_ave, e1_err, e2_ave, e2_err, e1_ave_n, e2_ave_n, file=out)
    print("Done with z_median =", z_median[j])

print("All done.")