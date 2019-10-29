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


## Mode 1: Using a configuration file 
##           Using original fits data
#config_file = 'corr_orig.yaml'
#config = treecorr.read_config(config_file)
#
#
#t1 = time.time()
#treecorr.corr2(config)
#t2 = time.time()
#print('Time for calculating gg correlation = ', t2-t1)

# Mode 2: Reading in an input Catalog

def autocorr(cat, outpath, nbins, mins, maxs, units, bin_slop, nthr):
    """
    funciton for auto-correlation
    """
    gg = treecorr.GGCorrelation(min_sep=mins, max_sep=maxs, 
                                nbins=nbins, sep_units=units, bin_slop=bin_slop)
    gg.process(cat, num_threads=nthr)
    gg.write(outpath)
    print("Succeed in writing", outpath)
    
def crosscorr(cat1, cat2, outpath, nbins, mins, maxs, units, bin_slop, nthr):
    """
    funciton for auto-correlation
    """
    gg = treecorr.GGCorrelation(min_sep=mins, max_sep=maxs, 
                                nbins=nbins, sep_units=units, bin_slop=bin_slop)
    gg.process(cat1, cat2, num_threads=nthr)
    gg.write(outpath)
    print("Succeed in writing", outpath)
    
if __name__ == "__main__":
    # input path
    inpathF = "/disks/shear15/ssli/KV450/selected/pre/bins/zb_"
    inpathL = ["13","35","57","79","912"]
    inpathP = ".h5"

    # weighted e
    we = np.loadtxt("/disks/shear15/ssli/KV450/CS/mine/TreeCorr/full/pre/e_vs_ZB.dat")

    # output path 
    outpathF = "/disks/shear15/ssli/KV450/CS/mine/TreeCorr/full/pre/slop/"
    outpathP = ".dat"
    
    # General parameters
    nbins = 9
    mins = 0.5
    maxs = 300.
    units = "arcmin"
    bin_slop = 0.05
    nthr = 8

    cat = []
    for i in range(len(inpathL)):
        s = inpathL[i]
        inpath = inpathF + s + inpathP
        data = pd.read_hdf(inpath, key='whole', 
                           columns=["ALPHA_J2000", "DELTA_J2000", 
                                    "bias_corrected_e1", "bias_corrected_e2", "recal_weight"])
        # weighted mean e for each bin
        e1_ave = we[i, 1]
        e2_ave = we[i, 3]
        # e1 & e2 used for TreeCorr
        e1_corr = data["bias_corrected_e1"] - e1_ave
        e2_corr = data["bias_corrected_e2"] - e2_ave

        cat.append(treecorr.Catalog(ra=data["ALPHA_J2000"], dec=data["DELTA_J2000"], 
                                    ra_units="deg", dec_units="deg", 
                                    w=data["recal_weight"],
                                    g1=e1_corr, g2=e2_corr))
        
        if s == "13":
            outpath = outpathF + "11" + outpathP
            p11 = mp.Process(target=autocorr, args=(cat[0], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p11.start()
            print("Start on 11...")

        elif s == "35":
            outpath = outpathF + "12" + outpathP
            p12 = mp.Process(target=crosscorr, args=(cat[0], cat[1], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p12.start()
            print("Start on 12...")

            outpath = outpathF + "22" + outpathP
            p22 = mp.Process(target=autocorr, args=(cat[1], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p22.start()
            print("Start on 22...")

        elif s == "57":
            outpath = outpathF + "13" + outpathP
            p13 = mp.Process(target=crosscorr, args=(cat[0], cat[2], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p13.start()
            print("Start on 13...")

            outpath = outpathF + "23" + outpathP
            p23 = mp.Process(target=crosscorr, args=(cat[1], cat[2], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p23.start()
            print("Start on 23...")

            outpath = outpathF + "33" + outpathP
            p33 = mp.Process(target=autocorr, args=(cat[2], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p33.start()
            print("Start on 33...")

        elif s == "79":
            outpath = outpathF + "14" + outpathP
            p14 = mp.Process(target=crosscorr, args=(cat[0], cat[3], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p14.start()
            print("Start on 14...")

            outpath = outpathF + "24" + outpathP
            p24 = mp.Process(target=crosscorr, args=(cat[1], cat[3], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p24.start()
            print("Start on 24...")

            outpath = outpathF + "34" + outpathP
            p34 = mp.Process(target=crosscorr, args=(cat[2], cat[3], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p34.start()
            print("Start on 34...")

            outpath = outpathF + "44" + outpathP
            p44 = mp.Process(target=autocorr, args=(cat[3], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p44.start()
            print("Start on 44...")

        elif s == "912":
            outpath = outpathF + "15" + outpathP
            p15 = mp.Process(target=crosscorr, args=(cat[0], cat[4], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p15.start()
            print("Start on 15...")

            outpath = outpathF + "25" + outpathP
            p25 = mp.Process(target=crosscorr, args=(cat[1], cat[4], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p25.start()
            print("Start on 25...")

            outpath = outpathF + "35" + outpathP
            p35 = mp.Process(target=crosscorr, args=(cat[2], cat[4], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p35.start()
            print("Start on 35...")

            outpath = outpathF + "45" + outpathP
            p45 = mp.Process(target=crosscorr, args=(cat[3], cat[4], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p45.start()
            print("Start on 45...")

            outpath = outpathF + "55" + outpathP
            p55 = mp.Process(target=autocorr, args=(cat[4], outpath, nbins, mins, maxs, units, bin_slop, nthr))
            p55.start()
            print("Start on 55...")
    
    p11.join()
    p12.join()    
    p13.join()    
    p14.join()
    p15.join()
    
    p22.join()    
    p23.join()    
    p24.join()
    p25.join()
    
    p33.join()    
    p34.join()    
    p35.join()
    
    p44.join()
    p45.join() 
    
    p55.join()
    
    print("All done.")
