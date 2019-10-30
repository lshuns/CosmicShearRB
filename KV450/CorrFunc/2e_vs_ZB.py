#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 16:21:10 2019

@author: ssli

weighted average of e for each tomographic bin and each patch
"""

import numpy as np
import pandas as pd
import multiprocessing as mp


def MeanFunc(patch, path, Z, bins, pq):
    """
    Function for weighted average calculation
    """

    hdf = pd.HDFStore(path, mode='r')
    print("Loaded data from", path)

    # boots
    nboot = 30
    e1w = np.zeros(nboot)
    e2w = np.zeros(nboot)

    for i in range(len(bins)-1):
    
        key = Z + '__' + str(bins[i]) + str(bins[i+1])

        data = hdf.select(key=key, 
            columns=["bias_corrected_e1", "bias_corrected_e2", "recal_weight"])
        print("Selected data from", key)


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

        e1_ave_b = np.mean(e1w)
        e1_err_b = np.std(e1w)
        e2_ave_b = np.mean(e2w)
        e2_err_b = np.std(e2w)

        # weight mean from direct calculation
        e1_ave = np.dot(e1, wt) / np.sum(wt)
        e2_ave = np.dot(e2, wt) / np.sum(wt)

        # output
        out = {"patch": patch, "key": key,
        'e1_ave': e1_ave, 'e2_ave': e2_ave, 
        'e1_ave_b': e1_ave_b, 'e1_err_b': e1_err_b, 
        'e2_ave_b': e2_ave_b, 'e2_err_b': e2_err_b}
        pq.put(out)

        print("Done with ", key, "for", patch)
    hdf.close()
    print("Done with", patch)


if __name__ == "__main__":

    bins = [1, 3, 5, 7, 9, 12]
    patches = ["G9","G12","G15","G23","GS"]
    # column used as redshift
    Z = 'Z_B'

    # input path
    pathF = "/disks/shear15/ssli/KV450/CorrFunc/data/"
    pathP = ".h5"

    # output path 
    log = open("/disks/shear15/ssli/KV450/CorrFunc/log/e_vs_ZB.csv", mode='w')

    # for mp
    jobs = []
    pq = mp.Queue()

    for patch in patches:
        path = pathF + patch + pathP

        p = mp.Process(target=MeanFunc, args=(patch, path, Z, bins, pq))
        jobs.append(p)
        p.start()

    for p in jobs:
        p.join()

    print("All processing done.")
    print("Start saving output.")

    # data information
    print("patch,key,e1_ave,e2_ave,e1_ave_b,e1_err_b,e2_ave_b,e2_err_b", file=log)
    while not pq.empty():
        tmp = pq.get()
        print(tmp['patch'], tmp['key'], tmp['e1_ave'], tmp['e2_ave'],
            tmp['e1_ave_b'], tmp['e1_err_b'], tmp['e2_ave_b'], tmp['e2_err_b'], sep=',', 
            file=log)

    print("All done.")
