#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 15:09:03 2019

@author: ssli

saving the binning results to feather forms

need to be added into 2bins.py in the future

"""

import pandas as pd


bins = [1, 3, 5, 7, 9, 12]
patches = ["G9","G12","G15","G23","GS"]
# column used as redshift
Z = 'Z_B'

# input path
inpathF = "/disks/shear15/ssli/KV450/CorrFunc/data/"
inpathP = ".h5"

# output path
outpathF = "/disks/shear15/ssli/KV450/CorrFunc/data/feather/"
outpathP = ".feather"


for patch in patches:

    inpath = inpathF + patch + inpathP
    hdf = pd.HDFStore(inpath, mode='r')
    print("hdf built from", inpath)

    for i in range(len(bins)-1):
    
        key = Z + '__' + str(bins[i]) + str(bins[i+1])
        savekey = '_' + str(bins[i]) + str(bins[i+1])

        data = hdf.get(key=key)
        print("Selected data from", key)

        outpath = outpathF + patch + savekey + outpathP

        data = data.reset_index(drop=True)
        data.to_feather(outpath)
        print("Saved data to", outpath)

    data = hdf.get(key='allBins')
    print("Selected data from allBins")

    outpath = outpathF + patch + '_allBins' + outpathP
    
    data = data.reset_index(drop=True)
    data.to_feather(outpath)
    print("Saved data to", outpath)
    
    hdf.close()
    print("New data saved for", patch)


# # check the saved data
# for patch in patches:

#     inpath = inpathF + patch + inpathP
#     hdf = pd.HDFStore(inpath, mode='r')
#     print("hdf built from", inpath)

#     for i in range(len(bins)-1):
    
#         key = Z + '__' + str(bins[i]) + str(bins[i+1])
#         savekey = '_' + str(bins[i]) + str(bins[i+1])

#         data = hdf.get(key=key)
#         print("Selected data from", key)

#         outpath = outpathF + patch + savekey + outpathP

#         data2 = pd.read_feather(outpath)
#         print("Loaded data from", outpath)

#         print(len(data)==len(data2))

#     data = hdf.get(key='allBins')
#     print("Selected data from allBins")

#     outpath = outpathF + patch + '_allBins' + outpathP
    
#     data2 = pd.read_feather(outpath)
#     print("Loaded data from", outpath)
#     print(len(data)==len(data2))
    
#     hdf.close()
