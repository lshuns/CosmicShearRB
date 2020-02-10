#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 11:56:24 2020

@author: ssli

script to run the correlation function using TreeCorr
"""


import multiprocessing as mp
import feather
import treecorr
import pandas as pd
import numpy as np


import time

import os
import sys
# Self-defined package
sys.path.insert(0,os.path.realpath('..')) 
from CorrFunc import MeanFunc, CorrFunc, CorrCosmoFunc, CorrPlotFunc


Start = time.time()

# General parameters
nzbins = 5
patches = ["G9","G12","G15","G23","GS"]

# parameters for treecorr
theta_nbins = 9
theta_min = 0.5
theta_max = 300.
theta_unit = "arcmin"
theta_bin_slop = 0.05
nthr = 8


# m value
mss = []
#
path_m_whole = "/disks/shear15/ssli/CosmicShear/shear_bias/Summary_multiplicative_whole.csv"
path_m_red = "/disks/shear15/ssli/CosmicShear/shear_bias/Summary_multiplicative_red.csv"
path_m_blue = "/disks/shear15/ssli/CosmicShear/shear_bias/Summary_multiplicative_blue.csv"
#
tmp = pd.read_csv(path_m_whole)
mss.append(tmp['m'].values)
tmp = pd.read_csv(path_m_red)
mss.append(tmp['m'].values)
tmp = pd.read_csv(path_m_blue)
mss.append(tmp['m'].values)


# input path
# data 
# whole / red / blue
inpathFs = ["/disks/shear15/ssli/KV450/tomo/",\
            "/disks/shear15/ssli/KV450/split/",\
            "/disks/shear15/ssli/KV450/split/"]
inpathPs =['.feather', \
            '_T_B_less3.feather', \
            '_T_B_greater3.feather']

# output path 
outpathF = "/disks/shear15/ssli/CosmicShear/data_vector/treecorr/"
outpathPs =['_whole.dat', \
            '_red.dat', \
            '_blue.dat']

# rearrange catalogue path
Reoutpaths_plot_withoutK = ["/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withoutK_whole.dat", \
                        "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withoutK_red.dat", \
                        "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withoutK_blue.dat"]
#
Reoutpaths_plot_withK = ["/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_whole.dat", \
                    "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_red.dat", \
                    "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_blue.dat"]
#
ReoutpathFs_cosmo = ["/disks/shear15/ssli/CosmicShear/data_vector/for_cosmo/xi_for_cosmo_", \
                        "/disks/shear15/ssli/CosmicShear/data_vector/for_cosmo/xi_for_cosmo_", \
                        "/disks/shear15/ssli/CosmicShear/data_vector/for_cosmo/xi_for_cosmo_"]
#
ReoutpathPs_cosmo = ["_whole.dat", "_red.dat", "_blue.dat"]


# c-term
log_cW = open("/disks/shear15/ssli/CosmicShear/shear_bias/e_vs_ZB_whole.csv", 'w')
print("patch,bin,e1_ave,e2_ave,e1_ave_b,e1_err_b,e2_ave_b,e2_err_b", file=log_cW)
log_cR = open("/disks/shear15/ssli/CosmicShear/shear_bias/e_vs_ZB_red.csv", 'w')
print("patch,bin,e1_ave,e2_ave,e1_ave_b,e1_err_b,e2_ave_b,e2_err_b", file=log_cR)
log_cB = open("/disks/shear15/ssli/CosmicShear/shear_bias/e_vs_ZB_blue.csv", 'w')
print("patch,bin,e1_ave,e2_ave,e1_ave_b,e1_err_b,e2_ave_b,e2_err_b", file=log_cB)
log_cs = [log_cW, log_cR, log_cB]


# jobs = []
# for kk in range(len(inpathFs)):
#     ms = mss[kk]

#     inpathF = inpathFs[kk]
#     inpathP = inpathPs[kk]
#     outpathP = outpathPs[kk]
#     log_c = log_cs[kk]

#     # build treecorr catalog
#     cat = []
#     for i in range(nzbins):
        
#         m = ms[i]
    
#         data = []
#         for patch in patches:

#             inpath = inpathF + patch + '_tomo' + str(i+1) + inpathP

#             tmp = feather.read_dataframe(inpath)

#             # c-term           
#             c_term = MeanFunc(tmp['bias_corrected_e1'], tmp['bias_corrected_e2'], tmp["recal_weight"])
#             print(patch, i+1, 
#                 c_term['e1_ave'], c_term['e2_ave'],
#                 c_term['e1_ave_b'], c_term['e1_err_b'], c_term['e2_ave_b'], c_term['e2_err_b'], 
#                 sep=',', file=log_c)

#             tmp['bias_corrected_e1'] -= c_term['e1_ave'] 
#             tmp['bias_corrected_e2'] -= c_term['e2_ave']

#             data.append(tmp)


#         data = pd.concat(data)
            
#         print("Data built for bin", str(i+1), "in loop", kk)

#         cat.append(treecorr.Catalog(ra=data["ALPHA_J2000"], dec=data["DELTA_J2000"], 
#                                     ra_units="deg", dec_units="deg", 
#                                     w=data["recal_weight"],
#                                     g1=data["bias_corrected_e1"], g2=data["bias_corrected_e2"]))
    
#     print("Treecorr catalogues built", "in loop", kk)
#     log_c.close()


#     # calculate correlation function
#     for idx_z1 in range(nzbins):
#         for idx_z2 in range(idx_z1, nzbins):

#             outpath = outpathF + 'tomo_' + str(idx_z1+1) + '_' + str(idx_z2+1) + outpathP

#             if idx_z1 == idx_z2:
#                 p = mp.Process(target=CorrFunc, args=(cat[idx_z1], None, outpath, 
#                                 theta_nbins, theta_min, theta_max, theta_unit, theta_bin_slop, nthr))

#             else:
#                 p = mp.Process(target=CorrFunc, args=(cat[idx_z1], cat[idx_z2], outpath, 
#                                 theta_nbins, theta_min, theta_max, theta_unit, theta_bin_slop, nthr))

#             jobs.append(p)

#             p.start()

#             print("Loop", kk, "start on",  (str(idx_z1+1) + str(idx_z2+1)), "...")

# for p in jobs:
#     p.join()

# print("Running Finished.")

# Rearrange the result
print("Start rearrange the results...")

for i in range(len(Reoutpaths_plot_withoutK)):
    Reoutpath_plot_withK = Reoutpaths_plot_withK[i]
    Reoutpath_plot_withoutK = Reoutpaths_plot_withoutK[i]

    ms = mss[i]

    outpathP = outpathPs[i]

    ReoutpathF_cosmo = ReoutpathFs_cosmo[i]
    ReoutpathP_cosmo = ReoutpathPs_cosmo[i]
    
    CorrCosmoFunc(Nbins=nzbins, ntheta=theta_nbins,
                        inpathF=outpathF, inpathP=outpathP, 
                        m_list=[],
                        outpathF=ReoutpathF_cosmo, outpathP=ReoutpathP_cosmo)
    CorrCosmoFunc(Nbins=nzbins, ntheta=theta_nbins,
                        inpathF=outpathF, inpathP=outpathP, 
                        m_list=ms,
                        outpathF=ReoutpathF_cosmo, outpathP=ReoutpathP_cosmo)

    CorrPlotFunc(Nbins=nzbins, 
                        inpathF=ReoutpathF_cosmo, inpathP=ReoutpathP_cosmo, 
                        withK=False,
                        outpath=Reoutpath_plot_withoutK)

    CorrPlotFunc(Nbins=nzbins, 
                        inpathF=ReoutpathF_cosmo, inpathP=ReoutpathP_cosmo, 
                        withK=True,
                        outpath=Reoutpath_plot_withK)


print("All finished in", time.time()-Start)