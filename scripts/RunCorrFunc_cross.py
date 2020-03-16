#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 11:56:24 2020

@author: ssli

Script to run the correlation function using TreeCorr
    Calculate the two subsamples jointly (cross-correlation required by covariance estimation)

Package:
    CorrFunc: MeanFunc, CorrFunc, CorrPlotFunc, CorrCosmoFunc, CorrErrCosmoFunc, CorrErrPlotFunc

Data location:
    Input: KV450/split
    Output: CosmicShear/data_vector/cross_subsamples
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
from CorrFunc import MeanFunc, CorrFunc, CorrPlotFunc, CorrCosmoFunc, CorrErrCosmoFunc, CorrErrPlotFunc


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

# input path
# data 
# red / blue
inpathF = "/disks/shear15/ssli/KV450/split/"
inpathPs =['_T_B_greater3.feather', \
            '_T_B_less3.feather']
# number of subsamples
n_subsamples = len(inpathPs)
# number of bins used for correlation function calculation 
ncorrbins = nzbins*n_subsamples


# output path 
outpathF = "/disks/shear15/ssli/CosmicShear/data_vector/cross_subsamples/xi_blue_red_"

jobs = []

for patch in patches:
    
    # build treecorr catalogue
    cat = [] # zbins organized as blue 5 then red 5
    for inpathP in inpathPs:

        for i in range(nzbins):

            inpath = inpathF + patch + '_tomo' + str(i+1) + inpathP

            tmp = feather.read_dataframe(inpath)

            # c-term           
            c_term = MeanFunc(tmp['bias_corrected_e1'], tmp['bias_corrected_e2'], tmp["recal_weight"])

            tmp['bias_corrected_e1'] -= c_term['e1_ave'] 
            tmp['bias_corrected_e2'] -= c_term['e2_ave']

            data = tmp
            
            cat.append(treecorr.Catalog(ra=data["ALPHA_J2000"], dec=data["DELTA_J2000"], 
                                        ra_units="deg", dec_units="deg", 
                                        w=data["recal_weight"],
                                        g1=data["bias_corrected_e1"], g2=data["bias_corrected_e2"]))
    
    print("Treecorr catalogues built for patch", patch)
    print("Catalogue length", len(cat))


    # calculate correlation function
    for idx_z1 in range(ncorrbins):
        for idx_z2 in range(idx_z1, ncorrbins):

                outpath = outpathF + 'tomo_' + str(idx_z1+1) + '_' + str(idx_z2+1) + '_' + patch + '.dat'

                if idx_z1 == idx_z2:
                    p = mp.Process(target=CorrFunc, args=(cat[idx_z1], None, outpath, 
                                    theta_nbins, theta_min, theta_max, theta_unit, theta_bin_slop, nthr))

                else:
                    p = mp.Process(target=CorrFunc, args=(cat[idx_z1], cat[idx_z2], outpath, 
                                    theta_nbins, theta_min, theta_max, theta_unit, theta_bin_slop, nthr))

                jobs.append(p)

                p.start()

                print('patch', patch, "bin",  (str(idx_z1+1) + '_' + str(idx_z2+1)), "started...")

for p in jobs:
    p.join()

print("Running Finished.")


# combine all the patches
header ="   r_nom       meanr       meanlogr       xip          xim         xip_im       xim_im     sigma_xip    sigma_xim      weight       npairs"
for idx_z1 in range(ncorrbins):
    for idx_z2 in range(idx_z1, ncorrbins):

        outpath = outpathF + 'tomo_' + str(idx_z1+1) + '_' + str(idx_z2+1) + '.dat'

        xi_list = []
        for patch in patches:
            inpath = outpathF + 'tomo_' + str(idx_z1+1) + '_' + str(idx_z2+1) + '_' + patch + '.dat'

            xi_list.append(np.loadtxt(inpath))

        xi_out = np.zeros(np.shape(xi_list[0]))

        xis = np.dstack(xi_list)
        
        wg = 1./xis[:,7]**2. # the error for xi
        for id_col in range(len(xi_out[:])):
            if (id_col != 7) and (id_col != 8) and (id_col != 9) and (id_col != 10):
                # weight average for parameters
                xi_out[:, id_col] = np.average(xis[:,id_col], axis=1, weights=wg)
            elif (id_col==7) or (id_col==8):
                # for errors
                xi_out[:,id_col] = np.sqrt(1./np.sum(1./xis[:,id_col]**2, axis=1))
            else:
                # for weights and npairs
                xi_out[:,id_col] = np.sum(xis[:,id_col], axis=1)

        np.savetxt(outpath, xi_out, fmt='%.4e', header=header)


print("All finished in", time.time()-Start)
# eemmeer (March 16, 2020)
# All finished in 158.4673569202423
