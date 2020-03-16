#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 11:56:24 2020

@author: ssli

Script to run the correlation function using TreeCorr
    Analyse different samples separately

Package:
    CorrFunc: MeanFunc, CorrFunc, CorrPlotFunc, CorrCosmoFunc, CorrErrCosmoFunc, CorrErrPlotFunc

Data location:
    Input: 
        m-value: CosmicShear/shear_bias/
        data: KV450/tomo, KV450/split
    Output: CosmicShear/data_vector/
        original: treecorr/
        rearranged: for_cosmo, for_plot
        c-term: CosmicShear/shear_bias/
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


# m value
mss = []
#
path_m_whole = "/disks/shear15/ssli/CosmicShear/shear_bias/Summary_m_whole.csv"
path_m_red = "/disks/shear15/ssli/CosmicShear/shear_bias/Summary_m_less3.csv"
path_m_blue = "/disks/shear15/ssli/CosmicShear/shear_bias/Summary_m_greater3.csv"
#
tmp = pd.read_csv(path_m_whole)
ms = tmp.sort_values(by='bin')['m'].values
mss.append(ms)
tmp = pd.read_csv(path_m_red)
ms = tmp.sort_values(by='bin')['m'].values
mss.append(ms)
tmp = pd.read_csv(path_m_blue)
ms = tmp.sort_values(by='bin')['m'].values
mss.append(ms)

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
            '_less3.dat', \
            '_greater3.dat']

# rearrange catalogue path
#
ReoutpathFs_cosmo = ["/disks/shear15/ssli/CosmicShear/data_vector/for_cosmo/xi_for_cosmo_", \
                        "/disks/shear15/ssli/CosmicShear/data_vector/for_cosmo/xi_for_cosmo_", \
                        "/disks/shear15/ssli/CosmicShear/data_vector/for_cosmo/xi_for_cosmo_"]
#
ReoutpathPs_cosmo = ["_whole.dat", "_less3.dat", "_greater3.dat"]
#
Reoutpaths_plot_withoutK = ["/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withoutK_whole.dat", \
                        "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withoutK_less3.dat", \
                        "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withoutK_greater3.dat"]
#
Reoutpaths_plot_withK = ["/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_whole.dat", \
                    "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_less3.dat", \
                    "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_greater3.dat"]
#
Reoutpaths_plot_err = ["/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_whole.dat", \
                    "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_less3.dat", \
                    "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_greater3.dat"]


# c-term
log_cW = open("/disks/shear15/ssli/CosmicShear/shear_bias/e_vs_ZB_whole.csv", 'w')
print("patch,bin,e1_ave,e2_ave,e1_ave_b,e1_err_b,e2_ave_b,e2_err_b", file=log_cW)
log_cR = open("/disks/shear15/ssli/CosmicShear/shear_bias/e_vs_ZB_less3.csv", 'w')
print("patch,bin,e1_ave,e2_ave,e1_ave_b,e1_err_b,e2_ave_b,e2_err_b", file=log_cR)
log_cB = open("/disks/shear15/ssli/CosmicShear/shear_bias/e_vs_ZB_greater3.csv", 'w')
print("patch,bin,e1_ave,e2_ave,e1_ave_b,e1_err_b,e2_ave_b,e2_err_b", file=log_cB)
log_cs = [log_cW, log_cR, log_cB]


jobs = []
for kk in range(len(inpathFs)):

    inpathF = inpathFs[kk]
    inpathP = inpathPs[kk]
    outpathP = outpathPs[kk]
    log_c = log_cs[kk]

    # build treecorr catalog
    cat = [] # patches followed by zbins
    for patch in patches:
        for i in range(nzbins):

            inpath = inpathF + patch + '_tomo' + str(i+1) + inpathP

            tmp = feather.read_dataframe(inpath)

            # c-term           
            c_term = MeanFunc(tmp['bias_corrected_e1'], tmp['bias_corrected_e2'], tmp["recal_weight"])
            print(patch, i+1, 
                c_term['e1_ave'], c_term['e2_ave'],
                c_term['e1_ave_b'], c_term['e1_err_b'], c_term['e2_ave_b'], c_term['e2_err_b'], 
                sep=',', file=log_c)

            tmp['bias_corrected_e1'] -= c_term['e1_ave'] 
            tmp['bias_corrected_e2'] -= c_term['e2_ave']

            data = tmp
            
            print("Data built for patch", patch, "bin", str(i+1), "in loop", kk)

            cat.append(treecorr.Catalog(ra=data["ALPHA_J2000"], dec=data["DELTA_J2000"], 
                                        ra_units="deg", dec_units="deg", 
                                        w=data["recal_weight"],
                                        g1=data["bias_corrected_e1"], g2=data["bias_corrected_e2"]))
    
    print("Treecorr catalogues built", "in loop", kk)
    log_c.close()

    # calculate correlation function
    for id_patch in range(len(patches)):
        patch = patches[id_patch]

        for idx_z1_ori in range(nzbins):
            idx_z1 = idx_z1_ori + id_patch*nzbins
            for idx_z2_ori in range(idx_z1_ori, nzbins):
                idx_z2 = idx_z2_ori + id_patch*nzbins

                outpath = outpathF + 'tomo_' + str(idx_z1_ori+1) + '_' + str(idx_z2_ori+1) + '_' + patch + outpathP

                if idx_z1_ori == idx_z2_ori:
                    p = mp.Process(target=CorrFunc, args=(cat[idx_z1], None, outpath, 
                                    theta_nbins, theta_min, theta_max, theta_unit, theta_bin_slop, nthr))

                else:
                    p = mp.Process(target=CorrFunc, args=(cat[idx_z1], cat[idx_z2], outpath, 
                                    theta_nbins, theta_min, theta_max, theta_unit, theta_bin_slop, nthr))

                jobs.append(p)

                p.start()

                print("Loop", kk, 'patch', patch, "bin",  (str(idx_z1+1) + '_' + str(idx_z2+1)), "started...")

for p in jobs:
    p.join()

print("Running Finished.")


# combine all the patches
header ="   r_nom       meanr       meanlogr       xip          xim         xip_im       xim_im     sigma_xip    sigma_xim      weight       npairs"
for kk in range(len(inpathFs)):

    outpathP = outpathPs[kk]

    for idx_z1 in range(nzbins):
        for idx_z2 in range(idx_z1, nzbins):

            outpath = outpathF + 'tomo_' + str(idx_z1+1) + '_' + str(idx_z2+1) + outpathP

            xi_list = []
            for patch in patches:
                inpath = outpathF + 'tomo_' + str(idx_z1+1) + '_' + str(idx_z2+1) + '_' + patch + outpathP

                xi_list.append(np.loadtxt(inpath))

            print('xi_list in combining patches', xi_list)
            xi_out = np.zeros(np.shape(xi_list[0]))

            xis = np.dstack(xi_list)
            print('xis in combining patches', xis)

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


# Rearrange the result
print("Start rearrange the results...")

for i in range(len(Reoutpaths_plot_withoutK)):

    ms = mss[i]

    outpathP = outpathPs[i]

    ReoutpathF_cosmo = ReoutpathFs_cosmo[i]
    ReoutpathP_cosmo = ReoutpathPs_cosmo[i]

    Reoutpath_plot_withK = Reoutpaths_plot_withK[i]
    Reoutpath_plot_withoutK = Reoutpaths_plot_withoutK[i]
    Reoutpath_plot_err = Reoutpaths_plot_err[i]
    
    CorrCosmoFunc(Nbins=nzbins, ntheta=theta_nbins,
                        inpathF=outpathF, inpathP=outpathP, 
                        m_list=[],
                        outpathF=ReoutpathF_cosmo, outpathP=ReoutpathP_cosmo)
    CorrCosmoFunc(Nbins=nzbins, ntheta=theta_nbins,
                        inpathF=outpathF, inpathP=outpathP, 
                        m_list=ms,
                        outpathF=ReoutpathF_cosmo, outpathP=ReoutpathP_cosmo)

    CorrErrCosmoFunc(Nbins=nzbins, ntheta=theta_nbins,
                        inpathF=outpathF, inpathP=outpathP, 
                        outpathF=ReoutpathF_cosmo, outpathP=ReoutpathP_cosmo)

    CorrPlotFunc(Nbins=nzbins, 
                        inpathF=ReoutpathF_cosmo, inpathP=ReoutpathP_cosmo, 
                        withK=False,
                        outpath=Reoutpath_plot_withoutK)

    CorrPlotFunc(Nbins=nzbins, 
                        inpathF=ReoutpathF_cosmo, inpathP=ReoutpathP_cosmo, 
                        withK=True,
                        outpath=Reoutpath_plot_withK)

    CorrErrPlotFunc(Nbins=nzbins, 
                        inpathF=ReoutpathF_cosmo, inpathP=ReoutpathP_cosmo, 
                        outpath=Reoutpath_plot_err)



print("All finished in", time.time()-Start)
# eemmeer (March 16, 2020)
# All finished in 279.2406942844391
