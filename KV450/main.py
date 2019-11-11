#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:05:44 2019

@author: ssli

Scripts for running the whole Cosmic Shear Analysis
"""

import Split
import CorrFunc


# ++++++++++ Data have been preprocessed & splitted 
import feather
import treecorr
import multiprocessing as mp
import pandas as pd
import numpy as np
import time 

# # ++++++++++++++++++++ CorrFunc 
# start_time = time.time()

# bins = [1, 3, 5, 7, 9, 12]
# Nbins = 5
# patches = ["G9","G12","G15","G23","GS"]

# # input path
# # splited data
# inpathF = "/disks/shear15/ssli/KV450/Split/data/halfHalf/"
# # samples
# samples = ['head', 'tail']

# # log
# log = "/disks/shear15/ssli/KV450/selected/log/"

# # output path 
# outpathF = "/disks/shear15/ssli/KV450/CorrFunc/full/bins_"
# outpathP = ".dat"
# outpathF_r = "/disks/shear15/ssli/KV450/CorrFunc/results_"

# # General parameters
# nbins = 9
# mins = 0.5
# maxs = 300.
# units = "arcmin"
# bin_slop = 0.05
# nthr = 8

# # running
# for sample in samples:
#     # c-term
#     pathc = log+"e_vs_ZB_"+sample+".csv"
#     log_c = open(pathc, 'w')    
#     print("patch,key,e1_ave,e2_ave,e1_ave_b,e1_err_b,e2_ave_b,e2_err_b", file=log_c)

#     # build treecorr catalog
#     cat = []
#     for i in range(len(bins)-1):
#         key = str(bins[i]) + str(bins[i+1]) + '_' + sample
    
#         data = []
#         for patch in patches:
#             inpath = inpathF + patch + '_' + key + '.feather'

#             tmp = feather.read_dataframe(inpath)

#             # c-term           
#             c_term = CorrFunc.MeanFunc(tmp['bias_corrected_e1'], tmp['bias_corrected_e2'], tmp["recal_weight"])
#             print(patch, key, 
#                 c_term['e1_ave'], c_term['e2_ave'],
#                 c_term['e1_ave_b'], c_term['e1_err_b'], c_term['e2_ave_b'], c_term['e2_err_b'], 
#                 sep=',', file=log_c)

#             tmp['bias_corrected_e1'] -= c_term['e1_ave'] 
#             tmp['bias_corrected_e2'] -= c_term['e2_ave']

#             data.append(tmp)

#         data = pd.concat(data)
#         print("Data built for", key)

#         cat.append(treecorr.Catalog(ra=data["ALPHA_J2000"], dec=data["DELTA_J2000"], 
#                                     ra_units="deg", dec_units="deg", 
#                                     w=data["recal_weight"],
#                                     g1=data["bias_corrected_e1"], g2=data["bias_corrected_e2"]))
#     print("Treecorr catalogues built.")
#     log_c.close()

#     # calculate correlation function
#     jobs = []
#     for i in range(len(cat)):
#         for j in range(i+1):
#             outpath = outpathF + str(j+1) + str(i+1) + '_' + sample + outpathP
#             if i == j:
#                 p = mp.Process(target=CorrFunc.CorrFunc, args=(cat[j], None, outpath, 
#                                 nbins, mins, maxs, units, bin_slop, nthr))
#             else:
#                 p = mp.Process(target=CorrFunc.CorrFunc, args=(cat[j], cat[i], outpath, 
#                                 nbins, mins, maxs, units, bin_slop, nthr))

#             jobs.append(p)
#             p.start()
#             print("Start on",  (str(j+1) + str(i+1)), "...")

#     for p in jobs:
#         p.join()

#     # rearrange the result
    
#     inpathF_r = outpathF
#     inpathP = '_' + sample + outpathP

#     outpath = "/disks/shear15/ssli/KV450/CorrFunc/results_" + sample + ".csv"
#     CorrFunc.CorrRearrangeFunc(Nbins, inpathF_r, inpathP, outpath)

#     print("Done for", sample)

# print("Running time", time.time()-start_time, 'seconds')

# # plot
# # custom settings for plot
# XLIM_P = [0.1, 300]
# XLIM_M = [5, 350]
# YLIM = [-4, 6]
# # color
# CR = ['black', 'red', 'blue'] 
# # marker
# MK = ['o', 'o', 'o']
# # marker size
# MS = 4
# # capthick
# MW = 0.5
# # linestyle
# LS = ['-', '--', '-.']
# # linewidth
# LW_H = 1.0
# LW_G = 0.5

# # save as pdf or png
# pdfORpng = 'png'

# # Number of bins 
# N_bins = 5

# # input directory
# inpath_whole = "/disks/shear15/ssli/KV450/CorrFunc/results_whole.csv"
# inpath_head = "/disks/shear15/ssli/KV450/CorrFunc/results_head.csv"
# inpath_tail = "/disks/shear15/ssli/KV450/CorrFunc/results_tail.csv"

# inpaths = [inpath_whole, inpath_head, inpath_tail]
# names = ['whole', 'small $T_B$', 'big $T_B$']

# # output directory
# outpath1 = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/split_"
# outpath2 = "/disks/shear15/ssli/KV450/CorrFunc/split_"
# outpaths = [outpath1, outpath2]

# CorrFunc.xiPlotFunc(inpaths, names, outpaths, N_bins, pdfORpng, XLIM_P, YLIM, CR, MK, MS, MW, LS, LW_H, LW_G, 'P')
# CorrFunc.xiPlotFunc(inpaths, names, outpaths, N_bins, pdfORpng, XLIM_M, YLIM, CR, MK, MS, MW, LS, LW_H, LW_G, 'M')



# ++++++++++++++++++++++++++++++++++++++++++++++ HistFunc
bins = ["13", "35", "57", "79", "912"]
patches = ["G9","G12","G15","G23","GS"]

# input directory
inpathF = "/disks/shear15/ssli/KV450/Split/data/halfHalf/"
inpathP = ".feather"

para_n = "T_B"
wg_n = "recal_weight"

# output directory
outpathF1 = "/disks/shear15/ssli/KV450/Split/hist/hist_TB_"
outpathF2 = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/Split/hist/hist_TB_"
outpathP = ".png"

# plot related
XLABEL = "$T_B$"
YLABEL = "Distribution"
# plot related
COLORS = ['red', 'blue']
NBINS = [60, 60]
LABELS = ['small $T_B$', 'big $T_B$']

for Bin in bins:
    for patch in patches:
        inpath1 = inpathF + patch + "_" + Bin + "_head" + inpathP
        sample1 = feather.read_dataframe(inpath1)
        print("Loaded data from", inpath1)

        inpath2 = inpathF + patch + "_" + Bin + "_tail" + inpathP
        sample2 = feather.read_dataframe(inpath2)
        print("Loaded data from", inpath2)

        values = [sample1["T_B"], sample2["T_B"]]
        wgs = [sample1["recal_weight"], sample2["recal_weight"]]

        # output path 
        outpath1 = outpathF1 + patch + "_" + Bin + outpathP
        outpath2 = outpathF2 + patch + "_" + Bin + outpathP
        outpaths = [outpath1, outpath2]

        Split.HistFunc(values, wgs, outpaths,
                COLORS, NBINS, LABELS,
                XLABEL, YLABEL)


# ++++++++++ Running the whole pipeline from origin data