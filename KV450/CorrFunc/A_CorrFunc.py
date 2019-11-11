#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 16:09:28 2019

@author: ssli

MeanFunc:
weighted average of e for a given data set

CorrFunc:
Using TreeCorr to calculate the 2-point shear correlation function
e is modified using weighted e 
bin_slop is set to be 0.05

CorrRearrangeFunc:
Function for rearrange the correlation results
"""

def MeanFunc(e1, e2, wt):
    """
    Function for weighted average calculation
    """
    print("Start calculation of weighted average calculation (MeanFunc)...")

    # boots
    nboot = 30
    e1w = np.zeros(nboot)
    e2w = np.zeros(nboot)

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
    out = {'e1_ave': e1_ave, 'e2_ave': e2_ave, 
        'e1_ave_b': e1_ave_b, 'e1_err_b': e1_err_b, 
        'e2_ave_b': e2_ave_b, 'e2_err_b': e2_err_b}

    print("Finished calculation of weighted average calculation (MeanFunc).")
    return out


def CorrFunc(cat1, cat2, outpath, nbins, mins, maxs, units, bin_slop, nthr):
    """
    funciton for correlation calculation (auto & cross)
    """
    print("Start correlation calculation (CorrFunc)...")

    gg = treecorr.GGCorrelation(min_sep=mins, max_sep=maxs, 
                                nbins=nbins, sep_units=units, bin_slop=bin_slop)

    if cat2 == None:
        # auto-correlation
        gg.process(cat1, num_threads=nthr)
    else:
        # cross-correlation
        gg.process(cat1, cat2, num_threads=nthr)
    gg.write(outpath)
    print("Correlation results saved in", outpath)
    print("Finished correlation calculation (CorrFunc).")


def CorrRearrangeFunc(Nbins, inpathF, inpathP,outpath):
    """
    Function for rearrange the correlation results
    """
    print("Start rearrange correlation results (CorrRearrangeFunc)...")

    res = pd.DataFrame(columns=['r_nom', 'meanr', 'meanlogr', 'xi_pm_real', 'xi_pm_im', 'sigma_pm', 'p1m2', 'itomo', 'jtomo'])


    for i in range(Nbins):
        j = i + 1
        k = j

        while k <= Nbins:
            inpath = inpathF + str(j) + str(k) + inpathP
            dat = np.loadtxt(inpath)

            # xi_p
            tmp = pd.DataFrame({'r_nom': dat[0:7,0],
                                'meanr': dat[0:7,1],
                                'meanlogr': dat[0:7,2],
                                'xi_pm_real': dat[0:7,3],
                                'xi_pm_im': dat[0:7,5],
                                'sigma_pm': dat[0:7,7],
                                'p1m2': np.full((7,),1),
                                'itomo': np.full((7,),1) * j,
                                'jtomo': np.full((7,),1) * k})
            res = pd.concat([res, tmp], ignore_index=True)

            # xi_m
            tmp = pd.DataFrame({'r_nom': dat[-6:,0],
                                'meanr': dat[-6:,1],
                                'meanlogr': dat[-6:,2],
                                'xi_pm_real': dat[-6:,4],
                                'xi_pm_im': dat[-6:,6],
                                'sigma_pm': dat[-6:,8],
                                'p1m2': np.full((6,),2),
                                'itomo': np.full((6,),1) * j,
                                'jtomo': np.full((6,),1) * k})
            res = pd.concat([res, tmp], ignore_index=True)

            k += 1

    res.to_csv(path_or_buf=outpath,
               header=True, index=False) 
    print("Results saved to", outpath)
    print("Finished rearrange correlation results (CorrRearrangeFunc).")


if __name__ == "__main__":

    import feather
    import treecorr
    import multiprocessing as mp
    import pandas as pd
    import numpy as np

    import time
    start_time = time.time()

    bins = [1, 3, 5, 7, 9, 12]
    patches = ["G9","G12","G15","G23","GS"]
    # column used as redshift
    Z = 'Z_B'

    # input path
    # data
    inpathF = "/disks/shear15/ssli/KV450/selected/"

    # output path 
    # c-term
    log_c = open("/disks/shear15/ssli/KV450/selected/log/e_vs_ZB.csv", 'w')
    print("patch,key,e1_ave,e2_ave,e1_ave_b,e1_err_b,e2_ave_b,e2_err_b", file=log_c)

    outpathF = "/disks/shear15/ssli/KV450/CorrFunc/full/bins_"
    outpathP = "_whole.dat"

    # General parameters
    nbins = 9
    mins = 0.5
    maxs = 300.
    units = "arcmin"
    bin_slop = 0.05
    nthr = 8

    # build treecorr catalog
    cat = []
    for i in range(len(bins)-1):
        key = Z + '__' + str(bins[i]) + str(bins[i+1])
    
        data = []
        for patch in patches:
            inpath = inpathF + patch + '__' + key + '.feather'

            tmp = feather.read_dataframe(inpath)

            # c-term           
            c_term = MeanFunc(tmp['bias_corrected_e1'], tmp['bias_corrected_e2'], tmp["recal_weight"])
            print(patch, key, 
                c_term['e1_ave'], c_term['e2_ave'],
                c_term['e1_ave_b'], c_term['e1_err_b'], c_term['e2_ave_b'], c_term['e2_err_b'], 
                sep=',', file=log_c)

            tmp['bias_corrected_e1'] -= c_term['e1_ave'] 
            tmp['bias_corrected_e2'] -= c_term['e2_ave']

            data.append(tmp)

        data = pd.concat(data)
        print("Data built for", key)

        cat.append(treecorr.Catalog(ra=data["ALPHA_J2000"], dec=data["DELTA_J2000"], 
                                    ra_units="deg", dec_units="deg", 
                                    w=data["recal_weight"],
                                    g1=data["bias_corrected_e1"], g2=data["bias_corrected_e2"]))
    print("Treecorr catalogues built.")
    log_c.close()

    # calculate correlation function
    jobs = []
    for i in range(len(cat)):
        for j in range(i+1):
            outpath = outpathF + str(j+1) + str(i+1) + outpathP
            if i == j:
                p = mp.Process(target=CorrFunc, args=(cat[j], None, outpath, 
                                nbins, mins, maxs, units, bin_slop, nthr))
            else:
                p = mp.Process(target=CorrFunc, args=(cat[j], cat[i], outpath, 
                                nbins, mins, maxs, units, bin_slop, nthr))

            jobs.append(p)
            p.start()
            print("Start on",  (str(j+1) + str(i+1)), "...")

    for p in jobs:
        p.join()


    # rearrange the result
    Nbins = 5
    inpathF = "/disks/shear15/ssli/KV450/CorrFunc/full/bins_"
    inpathP = "_whole.dat"

    outpath = "/disks/shear15/ssli/KV450/CorrFunc/results_whole.csv"
    CorrRearrangeFunc(Nbins, inpathF, inpathP, outpath)

    print("All done.")

    print("Running time", time.time()-start_time, 'seconds')
