#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 16:56:40 2019

@author: ssli

SplitByNumFunc:
split sample based on object numbers
"""

import pandas as pd 
import feather

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

# def splitByValFunc():
#     """
#     Function for sample split based on specific values
#     """

if __name__ == "__main__":

    bins = ["13", "35", "57", "79", "912"]
    patches = ["G9","G12","G15","G23","GS"]

    # input directory
    inpathF = "/disks/shear15/ssli/KV450/CorrFunc/data/feather/"
    inpathP = ".feather"

    para = "T_B"
    ratio = 0.5

    # output directory
    outpathF = "/disks/shear15/ssli/KV450/ColorSplit/data/halfHalf/"
    outpathP = ".feather"


    for Bin in bins:
        for patch in patches:
            inpath = inpathF + patch + "_" + Bin + inpathP
            data = feather.read_dataframe(inpath)
            print("Data loaded from", inpath)

            outdirF = outpathF + patch + "_" + Bin

            SplitByNumFunc(data, para, ratio, outdirF, outpathP)
