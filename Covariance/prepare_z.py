#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 14:24:29 2020

@author: ssli

Rearrange the redshift distribution file 
    required by covariance calculation
"""

import numpy as np 

z_min = ['0.1', '0.3', '0.5', '0.7', '0.9']
z_max = ['0.3', '0.5', '0.7', '0.9', '1.2']

inpathF = '/disks/shear15/ssli/CosmicShear/redshift/red_blue_renamed/Nz_DIR_z'
outDir = '/disks/shear15/ssli/CosmicShear/covariance/prepare/'
data = []
for i in range(len(z_min)):
    fname_b = inpathF + z_min[i] + 't' + z_max[i] + '_greater3.asc'
    data_b = np.loadtxt(fname_b)
    if i==0:
        # shift to middle
        shift_to_midpoint = (data_b[:,0][1]-data_b[:,0][0]) / 2.
        data.append(data_b[:,0]+shift_to_midpoint)

    data.append(data_b[:,1])

for i in range(len(z_min)):
    fname_r = inpathF + z_min[i] + 't' + z_max[i] + '_less3.asc'
    data_r = np.loadtxt(fname_r)

    data.append(data_r[:,1])

data = np.stack(data, axis=-1)
header = 'z_mid nz_bin1_blue nz_bin2_blue nz_bin3_blue nz_bin4_blue nz_bin5_blue nz_bin1_red nz_bin2_red nz_bin3_red nz_bin4_red nz_bin5_red'
np.savetxt(outDir+'Nz_DIR_shift_to_middle.asc', data, fmt=['%.3f', '%.6f', '%.6f', '%.6f', \
    '%.6f', '%.6f', '%.6f', '%.6f', '%.6f', '%.6f', '%.6f'], header=header)