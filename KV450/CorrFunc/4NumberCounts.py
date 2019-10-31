#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 09:46:10 2019

@author: ssli

Count object numbers in different bins

"""

import numpy as np 
import pandas as pd

def binCountFunc(path, Z, bins):
    """
    Function for number counting for different bins
    """

    dat = pd.read_csv(path)

    res = []
    for i in range(len(bins)-1):
        key = Z + '__' + str(bins[i]) + str(bins[i+1])
        
        N = np.sum(dat[dat.bin==key]['Number'].values)
        print("Number in bin (", str(bins[i] * 1e-1), ",", str(bins[i+1] * 1e-1), "]:", N)

        res.append(N)

    N_tot = np.sum(res)
    print("Total Number in all bins:", N_tot)
    res.append(N_tot)

    return np.array(res)


if __name__ == '__main__':
    
    bins = [1, 3, 5, 7, 9, 12]
    # column used as redshift
    Z = 'Z_B'

    # path directory
    path = "/disks/shear15/ssli/KV450/CorrFunc/log/log_bins.csv"

    res = binCountFunc(path, Z, bins)

