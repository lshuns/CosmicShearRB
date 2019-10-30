#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 16:09:28 2019

@author: ssli

Using TreeCorr to calculate the 2-point shear correlation function
e is modified using weighted e from 2e_vs_ZB.py
bin_slop is set to be 0.05
"""

import treecorr
import multiprocessing as mp
import pandas as pd
import numpy as np


def CorrFunc(bins, cat1, cat2, outpath, nbins, mins, maxs, units, bin_slop, nthr):
    """
    funciton for correlation calculation (auto & cross)
    """

    print("Start correlation calculation for", bins)

    gg = treecorr.GGCorrelation(min_sep=mins, max_sep=maxs, 
                                nbins=nbins, sep_units=units, bin_slop=bin_slop)

    if cat2 == None:
        # auto-correlation
        gg.process(cat1, num_threads=nthr)
    else:
        # cross-correlation
        gg.process(cat1, cat2, num_threads=nthr)
    gg.write(outpath)
    print("Results saved in", outpath)


if __name__ == "__main__":

    bins = [1, 3, 5, 7, 9, 12]
    patches = ["G9","G12","G15","G23","GS"]
    # column used as redshift
    Z = 'Z_B'

    # input path
    # data
    inpathF = "/disks/shear15/ssli/KV450/CorrFunc/data/"
    inpathP = ".h5"
    # c-term
    pathC = "/disks/shear15/ssli/KV450/CorrFunc/log/e_vs_ZB.csv"

    # output path 
    outpathF = "/disks/shear15/ssli/KV450/CorrFunc/full/bins_"
    outpathP = ".dat"

    # General parameters
    nbins = 9
    mins = 0.5
    maxs = 300.
    units = "arcmin"
    bin_slop = 0.05
    nthr = 8
    

    # weighted e
    we = pd.read_csv(pathC)

    # data
    hdfL = []
    for patch in patches:
        path = inpathF + patch + inpathP
        hdf = pd.HDFStore(path, mode='r')
        hdfL.append(hdf)


    # build treecorr catalog
    cat = []
    for i in range(len(bins)-1):
    
        key = Z + '__' + str(bins[i]) + str(bins[i+1])
        data = []
        for j in range(len(hdfL)):
            hdf = hdfL[j]
            patch = patches[j]

            tmp = hdf.select(key=key, 
                            columns=["ALPHA_J2000", "DELTA_J2000", 
                                    "bias_corrected_e1", "bias_corrected_e2", "recal_weight"])
            tmp['bias_corrected_e1'] -= we[(we.patch==patch) & (we.key==key)]['e1_ave'].values 
            tmp['bias_corrected_e2'] -= we[(we.patch==patch) & (we.key==key)]['e2_ave'].values

            data.append(tmp)

        data = pd.concat(data)

        print("Loaded data from", key)


        cat.append(treecorr.Catalog(ra=data["ALPHA_J2000"], dec=data["DELTA_J2000"], 
                                    ra_units="deg", dec_units="deg", 
                                    w=data["recal_weight"],
                                    g1=data["bias_corrected_e1"], g2=data["bias_corrected_e2"]))
    print("Treecorr catalogues built.")
    for hdf in hdfL:
        hdf.close()

    # calculate correlation function
    jobs = []
    for i in range(len(cat)):
        for j in range(i+1):
            outpath = outpathF + str(j+1) + str(i+1) + outpathP
            if i == j:
                p = mp.Process(target=CorrFunc, args=((str(j+1) + str(i+1)), cat[j], None, 
                                outpath, nbins, mins, maxs, units, bin_slop, nthr))
            else:
                p = mp.Process(target=CorrFunc, args=((str(j+1) + str(i+1)), cat[j], cat[i], 
                                outpath, nbins, mins, maxs, units, bin_slop, nthr))

            jobs.append(p)
            p.start()
            print("Start on",  (str(j+1) + str(i+1)), "...")


    for p in jobs:
        p.join()

    print("All done.")
