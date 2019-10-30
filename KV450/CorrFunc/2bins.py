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

"""


import pandas as pd
import multiprocessing as mp


def BinFunc(patch, path, Z, bins, pq):
    """
    Function for redshift binning
    """

    hdf = pd.HDFStore(path, mode='a')

    dat = hdf.get(key='whole')
    print("Loaded data from", path)

    for i in range(len(bins)-1):
    
        zmin = bins[i] * 1e-1
        zmax = bins[i+1] * 1e-1
        tmp = dat[(dat[Z] > zmin) & (dat[Z] <= zmax)]

        save_key = Z + '__' + str(bins[i]) + str(bins[i+1])

        hdf.put(key=save_key, value=tmp, format='table', data_columns=True)
        
        # selected numbers
        N_bin = len(tmp)

        # data information
        logdata = {"patch": patch, 'bin': save_key, 'N': N_bin}
        pq.put(logdata)

        print(patch, "succeed in selection and save for bin (", zmin, ",", zmax, "]")

    hdf.close()

    print("Finished binning for", patch)
    

if __name__ == "__main__":
    
    bins = [1, 3, 5, 7, 9, 12]
    patches = ["G9","G12","G15","G23","GS"]
    # column used as redshift
    Z = 'Z_B'

    # path directory
    pathF = "/disks/shear15/ssli/KV450/CorrFunc/data/"
    pathP = ".h5"

    ### running information
    log = open("/disks/shear15/ssli/KV450/CorrFunc/log/log_bins.csv", "w")

    # for mp
    jobs = []
    pq = mp.Queue()

    for patch in patches:
        path = pathF + patch + pathP

        p = mp.Process(target=BinFunc, args=(patch, path, Z, bins, pq))
        jobs.append(p)
        p.start()

    for p in jobs:
        p.join()

    print("All processing done.")
    print("Start saving data information.")

    # data information
    print("patch,bin,Number", file=log)
    while not pq.empty():
        tmp = pq.get()
        print(tmp["patch"], tmp['bin'], tmp['N'], sep=',', file=log)

    print("All done.")
