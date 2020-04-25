#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 16:41:11 2020
@author: ssli
Post process of covariance calculation
    re-arrange the results to desired forms
"""
import numpy as np
import pandas as pd

# cut matrix into three different parts:
#       bb, rr, rb (=br)

inpath = "/disks/shear15/ssli/CosmicShear/covariance/apr8_new/original/thps_cov_apr8_list.dat"
outDir = "/disks/shear15/ssli/CosmicShear/covariance/apr8_new/"

tmp_raw = np.loadtxt(inpath)

# discard undesired columns
tmp_raw = np.delete(tmp_raw, [2, 3, 8, 9], axis=1)

# build dataframe
df = pd.DataFrame(tmp_raw, \
    columns=['s1_bin1','s1_bin2','s1_xip0_xim1', 's1_theta', \
            's2_bin1','s2_bin2','s2_xip0_xim1', 's2_theta', \
            'Gaussian', 'non_Gaussian'])

df = df.astype(dtype= {"s1_bin1":"int32",\
                        "s1_bin2":"int32",\
                        "s1_xip0_xim1":"int32",\
                        "s1_theta":"float64",\
                        "s2_bin1":"int32",\
                        "s2_bin2":"int32",\
                        "s2_xip0_xim1":"int32",\
                        "s2_theta":"float64",\
                        "Gaussian":"float64",\
                        "non_Gaussian":"float64",\
                        })

# get blue-blue part
mask_bb = (df.s1_bin1<=5) & (df.s1_bin2<=5) & (df.s2_bin1<=5) & (df.s2_bin2<=5)
df_bb = df[mask_bb]

# get red-red part
mask_rr = (df.s1_bin1>5) & (df.s1_bin2>5) & (df.s2_bin1>5) & (df.s2_bin2>5)
df_rr = df[mask_rr]

# get blue-red part
mask_br = (df.s1_bin1<=5) & (df.s1_bin2<=5) & (df.s2_bin1>5) & (df.s2_bin2>5)
df_br = df[mask_br]


# output
outpath = outDir + 'thps_cov_apr8_bb_list.dat'
df_bb.to_csv(outpath, sep=' ', index=False, header=False)
#
outpath = outDir + 'thps_cov_apr8_rr_list.dat'
df_rr.to_csv(outpath, sep=' ', index=False, header=False)
#
outpath = outDir + 'thps_cov_apr8_br_list.dat'
df_br.to_csv(outpath, sep=' ', index=False, header=False)