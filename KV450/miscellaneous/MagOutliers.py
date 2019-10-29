#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 13:01:39 2019

@author: ssli

Count number of outliers in terms of magnitude

Things need to be improved: better output
"""

import pandas as pd

def Counts(s_mag, v_mag, s_path, inpathF, inpathP):
    inpath = inpathF + s_path + inpathP
    dat = pd.read_hdf(inpath,key='whole')
    return (dat[s_mag].values == v_mag).sum()


# path directory
inpathL = ["G9","G12","G15","G23","GS"]
inpathP = ".h5"
# without magnitude cut
inpathF = "/disks/shear15/ssli/KV450/selected/pre/"

# magnitude bands
smagL = ['MAG_GAAP_r_CALIB', 'MAG_GAAP_u_CALIB', 'MAG_GAAP_g_CALIB', 'MAG_GAAP_i_CALIB', 
        'MAG_GAAP_Z', 'MAG_GAAP_Y', 'MAG_GAAP_J', 'MAG_GAAP_H', 'MAG_GAAP_Ks']

# Number counts
c_m = 0
c_p = 0

for s_mag in smagL:
    tmp_m = 0
    tmp_p = 0
    for s_path in inpathL:
        tmp_m += Counts(s_mag, -99., s_path, inpathF, inpathP)
        tmp_p += Counts(s_mag,  99., s_path, inpathF, inpathP)
    
    print("Outliers in", s_mag, "(-99):", tmp_m)
    print("Outliers in", s_mag, "(99):", tmp_p)

    c_m += tmp_m
    c_p += tmp_p

print("Total outliers (-99):", c_m)
print("Total outliers (99):", c_p)
