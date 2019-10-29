#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 14:14:45 2019

@author: ssli

Plot the magnitude distribution of the selected data
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# path directory
inpathL = ["G9","G12","G15","G23","GS"]
inpathP = ".h5"
# with magnitude cut
inpathF_w = "/disks/shear15/ssli/KV450/selected/mine/"
# without magnitude cut
inpathF_wo = "/disks/shear15/ssli/KV450/selected/pre/"

def ReadData(s, inpathF, inpathP):
    inpath = inpathF + s + inpathP
    dat = pd.read_hdf(inpath,key='whole')
    return dat['MAG_GAAP_r_CALIB'].values
    
mag_w = np.concatenate([ReadData(s, inpathF_w, inpathP) for s in inpathL], axis=None)
mag_wo = np.concatenate([ReadData(s, inpathF_wo, inpathP) for s in inpathL], axis=None)

print(np.min(mag_wo))
print(np.max(mag_wo))

plt.figure()
plt.hist(x=mag_wo, bins=100, range=(0,50), color='blue', label='total', histtype='step')
plt.hist(x=mag_w, bins=100, range=(0,50), color='red', label='after magnitude cut', histtype='step')
plt.legend(frameon=False)
plt.xlabel('r magnitude')
plt.ylabel('Counts')
plt.savefig("/net/raam/data1/surfdrive_ssli/Projects/6wl_SimIm/plot/CorrFunc/hist_magr.pdf")
plt.close()
