#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 16:56:40 2019

@author: ssli

SplitByNumFunc (worn-out, use SplitByWgFunc instead):
split sample based on object numbers

SplitByWgFunc:
split sample based on object weights (wg)

SplitByValFunc:
split sample based on a specific value
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

def SplitByWgFunc(data, para, wg, ratio, outdirF, outdirP):
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

    outdir = outdirF + '_head_wg' + outdirP
    data1 = data1.reset_index(drop=True)
    data1.to_feather(outdir)
    print("Data saved to", outdir)

    outdir = outdirF + '_tail_wg' + outdirP
    data2 = data2.reset_index(drop=True)
    data2.to_feather(outdir)
    print("Data saved to", outdir)


def SplitByValFunc(data, para, cvalue, wg, outdir, savekey):
    """
    Function for sample split based on a specific value
    """
    print("Total weight for whole data", np.sum(data[wg].values))
    data1 = data[data[para]<=cvalue]
    data2 = data[data[para]>cvalue]
    print("Total weight for low part", np.sum(data1[wg].values))
    print("Total weight for high part", np.sum(data2[wg].values))

    outpath = outdir + savekey + para + '_' + str(cvalue) + '_less.feather'
    data1 = data1.reset_index(drop=True)
    data1.to_feather(outpath)
    print("Low data saved to", outpath)

    outpath = outdir + savekey + para + '_' + str(cvalue) + '_greater.feather'
    data2 = data2.reset_index(drop=True)
    data2.to_feather(outpath)
    print("High data saved to", outpath)


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

    # +++++++++++++++++++ SplitByValFunc
    # # SimIm (Bins)
    # Start = time.time()

    # indir = "/disks/shear15/ssli/SimCat/"
    # ingene = 'Bin_all_13_PSF__ZB9_in__'

    # bins = ['13', '35', '57', '79', '912']

    # para = 'TB9_in'
    # cvalue = 3
    # wg = 'LFweight'

    # # out
    # outdir = "/disks/shear15/ssli/SimCat/split/"

    # for Bin in bins:
    #     inpath = indir + ingene + Bin + '.feather'
    #     data = feather.read_dataframe(inpath)
    #     print("Data loaded from", inpath)

    #     savekey = ingene + Bin + '_'

    #     SplitByValFunc(data, para, cvalue, wg, outdir, savekey)

    # print("Finished in", time.time()-Start)
    # # Finished in 8.959846496582031


    # SimIm (whole)
    Start = time.time()

    indir = "/disks/shear15/ssli/SimCat/"
    ingene = 'SimCatSelec_all_13_PSF'

    para = 'TB9_in'
    cvalue = 3
    wg = 'LFweight'

    # out
    outdir = "/disks/shear15/ssli/SimCat/split/"

    inpath = indir + ingene + '.feather'
    data = feather.read_dataframe(inpath)
    print("Data loaded from", inpath)

    savekey = ingene + '_'

    SplitByValFunc(data, para, cvalue, wg, outdir, savekey)

    print("Finished in", time.time()-Start)
    # Finished in 12.693177700042725


    # # KV450
    # Start = time.time()

    # indir = "/disks/shear15/ssli/KV450/CorrFunc/data/feather/"

    # bins = ['13', '35', '57', '79', '912']
    # patches = ["G9","G12","G15","G23","GS"]

    # para = 'T_B'
    # cvalue = 3
    # wg = 'recal_weight'

    # # out
    # outdir = "/disks/shear15/ssli/KV450/Split/data/"

    # for Bin in bins:
    #     for patch in patches:
    #         inpath = indir + patch + "_" + Bin + '.feather'
    #         data = feather.read_dataframe(inpath)
    #         print("Data loaded from", inpath)

    #         savekey = patch + '_' + Bin + '_'

    #         SplitByValFunc(data, para, cvalue, wg, outdir, savekey)

    # print("Finished in", time.time()-Start)
    # # Finished in 




