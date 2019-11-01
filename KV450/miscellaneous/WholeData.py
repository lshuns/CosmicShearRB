#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 15:09:03 2019

@author: ssli

saving a whole data set by adding all bins objects

need to be added into 2bins.py in the future
"""

import pandas as pd


bins = [1, 3, 5, 7, 9, 12]
patches = ["G9","G12","G15","G23","GS"]
# column used as redshift
Z = 'Z_B'

# input path
pathF = "/disks/shear15/ssli/KV450/CorrFunc/data/"
pathP = ".h5"


for patch in patches:

    path = pathF + patch + pathP
    hdf = pd.HDFStore(path, mode='a')
    print("Loaded data from", path)

    data = []
    for i in range(len(bins)-1):
    
        key = Z + '__' + str(bins[i]) + str(bins[i+1])

        tmp = hdf.get(key=key)
        print("Selected data from", key)
        data.append(tmp)

    data = pd.concat(data)

    hdf.put(key='allBins', value=data, format='table', data_columns=True)
    hdf.close()
    print("New data saved for", patch)