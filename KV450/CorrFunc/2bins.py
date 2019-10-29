#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 14:53:56 2019

@author: ssli

build the five tomographic redshift bins:
    "13": (0.1, 0.3]
    "35": (0.3, 0.5]
    "57": (0.5, 0.7]
    "79": (0.7, 0.9]
    "912": (0.9, 1.2]

Things need to be improved: the form of the out log
"""


import pandas as pd
    
# path directory
inpathF = "/disks/shear15/ssli/KV450/selected/pre/"
inpathL = {"G9","G12","G15","G23","GS"}
inpathP = ".h5"

outpathF = "/disks/shear15/ssli/KV450/selected/pre/bins/zb_"
outpathL = {"13","35","57","79","912"} 
outpathP = ".h5"

### running information
log = open("/disks/shear15/ssli/KV450/selected/pre/bins/log.txt", "w")

# out bins
hdfL = []
for s in outpathL:
    outpath = outpathF + s + outpathP
    hdfL.append(pd.HDFStore(outpath, mode='w'))


for s in inpathL:
    print("In ", s)
    print(s, file=log)
    
    inpath = inpathF + s + inpathP
    dat = pd.read_hdf(inpath,key='whole')
    
    print("Total objects: ", len(dat), file=log)
    print("Succeed in loading data.")
    
    # Selection
    zmin = 0.1
    dz = 0.2
    for i in range(4):
        zmax = zmin + dz
        tmp = dat[(dat.Z_B > zmin) & (dat.Z_B <= zmax)]
        hdfL[i].put(key="whole", value=tmp, format='table',append=True, data_columns=True)
        
        print("Succeed in selection and save for bin (", zmin, ",", zmax, "]")
        print("Objects in bin (", zmin, ",", zmax, "]", len(tmp), file=log)

        zmin = zmax
        
    tmp = dat[(dat.Z_B > zmin) & (dat.Z_B <= 1.2)]
    hdfL[4].put(key="whole", value=tmp, format='table',append=True, data_columns=True)

    print("Succeed in selection and save for bin (", zmin, ", 1.2]")
    print("Objects in bin (", zmin, ", 1.2]", len(tmp), file=log)
    
# close all the hdf
for hdf in hdfL:
    hdf.close()