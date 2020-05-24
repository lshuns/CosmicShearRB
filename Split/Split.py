#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 16:56:40 2019

@author: ssli

Module to perform split of sources

SplitByNumFunc (worn-out, use SplitByWgFunc instead):
    Split sample based on object numbers.

SplitByWgFunc:
    Split sample based on object weights (wg) using assigned ratio.

SplitByValFunc (mp supported):
    Split sample based on a specific value.
    less: <= cvalue
    greater: > cvalue
"""

import pandas as pd 
import numpy as np
import feather

# worn-out
# try to use weight instead
def SplitByNumFunc(data, para, ratio, outdirF, outdirP):
    """
    Function for sample split based on object numbers
    """
    data.sort_values(by=[para], ascending=True, inplace=True)
    Ntot = len(data)
    N1 = int(Ntot*ratio)
    N2 = Ntot - N1

    data1 = data.head(n=N1)
    data2 = data.tail(n=N2)

    outdir = outdirF + '_head' + outdirP
    data1 = data1.reset_index(drop=True)
    data1.to_feather(outdir)
    print("Data saved to", outdir)

    outdir = outdirF + '_tail' + outdirP
    data2 = data2.reset_index(drop=True)
    data2.to_feather(outdir)
    print("Data saved to", outdir)


def SplitByWgFunc(data, para, wg, ratio, outdirF=None, outdirP=None):
    """
    Function for sample split based on object weights (wg)
    """
    data.sort_values(by=[para], ascending=True, inplace=True)    
    Wgtot = np.sum(data[wg].values)
    print("The total weight for whole data", Wgtot)
    Wg1 = Wgtot*ratio
    print("Desired total weight for the first part", Wg1)

    Ntot = len(data)
    N1_test = int(Ntot*ratio)
    Wg1_test = np.sum(data[wg][:N1_test].values)

    print("Guessed total weight for the first part", Wg1_test)

    if Wg1_test == Wg1:
        N1 = N1_test
    elif Wg1_test < Wg1:
        while Wg1_test < Wg1:
            N1_test += 1
            Wg1_test = np.sum(data[wg][:N1_test].values)
        N1 = N1_test
    elif Wg1_test > Wg1:
        while Wg1_test > Wg1:
            N1_test -= 1
            Wg1_test = np.sum(data[wg][:N1_test].values)
        N1 = N1_test
    
    N2 = Ntot - N1
    data1 = data.head(n=N1)
    data2 = data.tail(n=N2)

    print("Resulted total weights for the first part", np.sum(data1[wg].values))
    print("Resulted total weights for the latter part", np.sum(data2[wg].values))
    print("The weight for the edge object", data[wg][N1])
    print("The {:} value for the front edge object".format(para), data1[para].values[-1])
    print("The {:} value for the end edge object".format(para), data2[para].values[0])

    if outdirF != None:
        outdir = outdirF + '_head_wg' + outdirP
        data1 = data1.reset_index(drop=True)
        data1.to_feather(outdir)
        print("Data saved to", outdir)

        outdir = outdirF + '_tail_wg' + outdirP
        data2 = data2.reset_index(drop=True)
        data2.to_feather(outdir)
        print("Data saved to", outdir)


def SplitByValFunc(data, para_col, cvalue, 
        outdir, save_name_prefix, wg_col=None, pq=None):
    """
    Split sample based on a specific value.

    Parameters
    ----------
    data : pandas.DataFrame() or numpy.recarray
        Data being split.
    
    para_col: str
        Column name used as the split parameter
    
    cvalue: float
        Threshold value.

    outdir: str
        Directory for results.

    save_name_prefix: str
        File name prefix of the saved catalogues.

    wg_col: str, optional
        Column name used as weights. Required when log is desired.

    pq: multiprocessing.Queue() object, optional
        Object used for communication 
        between parallel processes.
    """

    data1 = data[data[para_col]<=cvalue]
    data2 = data[data[para_col]>cvalue]

    outpath = outdir + save_name_prefix + '_' + para_col + '_less' + str(cvalue) + '.feather'
    data1 = data1.reset_index(drop=True)
    data1.to_feather(outpath)
    print("Low data saved to", outpath)

    outpath = outdir + save_name_prefix + '_' + para_col + '_greater' + str(cvalue) + '.feather'
    data2 = data2.reset_index(drop=True)
    data2.to_feather(outpath)
    print("High data saved to", outpath)

    # data information
    if pq != None:
        # selected numbers
        wg_tot = np.sum(data[wg_col].values)
        N_tot = len(data)
        wg1 = np.sum(data1[wg_col].values)
        N1 = len(data1)
        wg2 = np.sum(data2[wg_col].values)
        N2 = len(data2)

        logdata = {"save_name_prefix": save_name_prefix, \
            'wg_tot': wg_tot, 'wg_less': wg1, 'wg_greater': wg2, \
            'number_tot': N_tot, 'number_less': N1, 'number_greater':  N2}
        pq.put(logdata)

if __name__ == "__main__":
    import time

    # +++++++++++++++++ SplitByNumFunc (worn out)
    # bins = ["13", "35", "57", "79", "912"]
    # patches = ["G9","G12","G15","G23","GS"]

    # # input directory
    # inpathF = "/disks/shear15/ssli/KV450/CorrFunc/data/feather/"
    # inpathP = ".feather"

    # para = "T_B"
    # ratio = 0.5

    # # output directory
    # outpathF = "/disks/shear15/ssli/KV450/ColorSplit/data/halfHalf/"
    # outpathP = ".feather"


    # for Bin in bins:
    #     for patch in patches:
    #         inpath = inpathF + patch + "_" + Bin + inpathP
    #         data = feather.read_dataframe(inpath)
    #         print("Data loaded from", inpath)

    #         outdirF = outpathF + patch + "_" + Bin

    #         SplitByNumFunc(data, para, ratio, outdirF, outpathP)



    # # +++++++++++++++++ SplitByWgFunc
    # bins = ["13", "35", "57", "79", "912"]
    # patches = ["G9","G12","G15","G23","GS"]

    # # input directory
    # inpathF = "/disks/shear15/ssli/KV450/CorrFunc/data/feather/"
    # inpathP = ".feather"

    # para = "T_B"
    # wg = "recal_weight"
    # ratio = 0.5

    # # output directory
    # outpathF = "/disks/shear15/ssli/KV450/Split/data/halfHalf/"
    # outpathP = ".feather"


    # for Bin in bins:
    #     for patch in patches:
    #         inpath = inpathF + patch + "_" + Bin + inpathP
    #         data = feather.read_dataframe(inpath)
    #         print("Data loaded from", inpath)

    #         outdirF = outpathF + patch + "_" + Bin

    #         SplitByWgFunc(data, para, wg, ratio, outdirF, outpathP)
