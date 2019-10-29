#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 11:59:18 2019

@author: ssli
"""

import pandas as pd 
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

import matplotlib as mpl
#mpl.use('Agg')

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True


#
#
#
#### read the data 
#inpath = "/disks/shear15/ssli/KV450/KV450_G12_reweight_3x4x4_v2_good.cat"
### with scope making the file be closed automatically
#with fits.open(inpath) as hdul:
#	# hdul.info() # information of whole hdul
#	# ### PrimaryHDU 
#	# ## Header
#	# hdr = hdul[0].header
#	# # print(repr(hdr))
#	# hdr_key = hdr.keys() # a list of all keywords
#	# print(list(hdr_key))
#
#	# ### Objects
#	# ## Header
#	# hdr = hdul[1].header 
#	# print(repr(hdr))
#	# hdr_key = hdr.keys() # a list of all keywords
#	# print(list(hdr_key))
#	## Data
#    data = hdul[1].data
#	# cols = hdul[1].columns # column information
#	# print(cols.info())
#	# print(cols.names)
#	# ['SeqNr', 'FLUX_AUTO_THELI', 'FLUXERR_AUTO_THELI', 'MAG_AUTO', 'MAGERR_AUTO_THELI', 'KRON_RADIUS_THELI', 
#	# 'BackGr_THELI', 'Level_THELI', 'MU_THRESHOLD_THELI', 'MaxVal_THELI', 'MU_MAX_THELI', 
#	# 'ISOAREA_WORLD_THELI', 'Xpos_THELI', 'Ypos_THELI', 'ALPHA_J2000', 'DELTA_J2000', 
#	# 'A_WORLD', 'B_WORLD', 'THETA_J2000', 'ERRA_WORLD_THELI', 'ERRA_IMAGE_THELI', 'ERRB_WORLD_THELI', 
#	# 'ERRB_IMAGE_THELI', 'THETA_SKY_THELI', 'ERRTHETA_SKY_THELI', 'THETA_WORLD_THELI', 'ERRTHETA_WORLD_THELI', 
#	# 'FWHM_IMAGE_THELI', 'FWHM_WORLD_THELI', 'Flag_THELI', 'FLUX_RADIUS_THELI', 'NIMAFLAGS_ISO_THELI', 
#	# 'CLASS_STAR_THELI', 'IMAFLAGS_ISO_THELI', 'X_WORLD_THELI', 'Y_WORLD_THELI', 'XY_WORLD_THELI', 
#	# 'X2_WORLD_THELI', 'Y2_WORLD_THELI', 'ERRX2_IMAGE_THELI', 'ERRX2_WORLD_THELI', 'ERRXY_WORLD_THELI', 
#	# 'ERRY2_WORLD_THELI', 'ERRXY_IMAGE_THELI', 'ERRY2_IMAGE_THELI', 'XM2_THELI', 'YM2_THELI', 'Corr_THELI', 
#	# 'CXX_IMAGE_THELI', 'CXY_IMAGE_THELI', 'CYY_IMAGE_THELI', 'CXX_WORLD_THELI', 'CXY_WORLD_THELI', 'CYY_WORLD_THELI', 
#	# 'ERRCXX_IMAGE_THELI', 'ERRCXY_IMAGE_THELI', 'ERRCYY_IMAGE_THELI', 'ERRCXX_WORLD_THELI', 'ERRCXY_WORLD_THELI', 
#	# 'ERRCYY_WORLD_THELI', 'NPIX_THELI', 'XMIN_IMAGE_THELI', 'XMAX_IMAGE_THELI', 'YMIN_IMAGE_THELI', 'YMAX_IMAGE_THELI', 
#	# 'A_THELI', 'B_THELI', 'POSANG_THELI', 'ERRTHETA_IMAGE_THELI', 'ELLIPTICITY_THELI', 'ELONGATION_THELI', 
#	# 'ISO0_THELI', 'ISO1_THELI', 'ISO2_THELI', 'ISO3_THELI', 'ISO4_THELI', 'ISO5_THELI', 'ISO6_THELI', 'ISO7_THELI', 
#	# 'EXTINCTION_u', 'EXTINCTION_g', 'EXTINCTION_r', 'EXTINCTION_i', 
#	# 'MAG_LIM_u', 'MAG_LIM_g', 'MAG_LIM_r', 'MAG_LIM_i', 
#	# 'MAG_ISO_THELI', 'MAGERR_ISO_THELI', 'MAG_ISOCOR_THELI', 'MAGERR_ISOCOR_THELI', 
#	# 'MAG_BEST_THELI', 'MAGERR_BEST_THELI', 'MAG_APER_THELI', 'MAGERR_APER_THELI', 
#	# 'FLUX_BEST_THELI', 'FLUXERR_BEST_THELI', 'FLUX_ISO_THELI', 'FLUXERR_ISO_THELI', 
#	# 'FLUX_ISOCOR_THELI', 'FLUXERR_ISOCOR_THELI', 
#	# 'FLUX_APER_4_THELI', 'FLUX_APER_5_THELI', 'FLUX_APER_6_THELI', 'FLUX_APER_7_THELI', 'FLUX_APER_8_THELI', 'FLUX_APER_9_THELI', 'FLUX_APER_10_THELI', 'FLUX_APER_11_THELI', 'FLUX_APER_12_THELI', 'FLUX_APER_13_THELI', 'FLUX_APER_14_THELI', 'FLUX_APER_15_THELI', 'FLUX_APER_16_THELI', 'FLUX_APER_17_THELI', 'FLUX_APER_18_THELI', 'FLUX_APER_19_THELI', 'FLUX_APER_20_THELI', 'FLUX_APER_25_THELI', 'FLUX_APER_30_THELI', 'FLUX_APER_35_THELI', 'FLUX_APER_40_THELI', 'FLUX_APER_45_THELI', 'FLUX_APER_50_THELI', 'FLUX_APER_55_THELI', 
#	# 'FLUXERR_APER_4_THELI', 'FLUXERR_APER_5_THELI', 'FLUXERR_APER_6_THELI', 'FLUXERR_APER_7_THELI', 'FLUXERR_APER_8_THELI', 'FLUXERR_APER_9_THELI', 'FLUXERR_APER_10_THELI', 'FLUXERR_APER_11_THELI', 'FLUXERR_APER_12_THELI', 'FLUXERR_APER_13_THELI', 'FLUXERR_APER_14_THELI', 'FLUXERR_APER_15_THELI', 'FLUXERR_APER_16_THELI', 'FLUXERR_APER_17_THELI', 'FLUXERR_APER_18_THELI', 'FLUXERR_APER_19_THELI', 'FLUXERR_APER_20_THELI', 'FLUXERR_APER_25_THELI', 'FLUXERR_APER_30_THELI', 'FLUXERR_APER_35_THELI', 'FLUXERR_APER_40_THELI', 'FLUXERR_APER_45_THELI', 'FLUXERR_APER_50_THELI', 'FLUXERR_APER_55_THELI', 
#	# 'Agaper_u', 'Bgaper_u', 'PAgaap_u', 'MAG_GAAP_u', 'MAGERR_GAAP_u', 'FLUX_GAAP_u', 'FLUXERR_GAAP_u', 
#	# 'Agaper_g', 'Bgaper_g', 'PAgaap_g', 'MAG_GAAP_g', 'MAGERR_GAAP_g', 'FLUX_GAAP_g', 'FLUXERR_GAAP_g', 
#	# 'Agaper_r', 'Bgaper_r', 'PAgaap_r', 'MAG_GAAP_r', 'MAGERR_GAAP_r', 'FLUX_GAAP_r', 'FLUXERR_GAAP_r', 
#	# 'Agaper_i', 'Bgaper_i', 'PAgaap_i', 'MAG_GAAP_i', 'MAGERR_GAAP_i', 'FLUX_GAAP_i', 'FLUXERR_GAAP_i', 
#	# 'Z_B_ugri', 'Z_B_MIN_ugri', 'Z_B_MAX_ugri', 'T_B_ugri', 
#	# 'ODDS_ugri', 'Z_ML_ugri', 'T_ML_ugri', 'CHI_SQUARED_BPZ_ugri', 'M_0_ugri', 'FIELD_POS', 
#	# 'MAG_GAAP_Z', 'MAGERR_GAAP_Z', 'FLUX_GAAP_Z', 'FLUXERR_GAAP_Z', 'GAAP_Flag_Z', 'GAAP_nexp_Z', 
#	# 'MAG_GAAP_Y', 'MAGERR_GAAP_Y', 'FLUX_GAAP_Y', 'FLUXERR_GAAP_Y', 'GAAP_Flag_Y', 'GAAP_nexp_Y', 
#	# 'MAG_GAAP_J', 'MAGERR_GAAP_J', 'FLUX_GAAP_J', 'FLUXERR_GAAP_J', 'GAAP_Flag_J', 'GAAP_nexp_J', 
#	# 'MAG_GAAP_H', 'MAGERR_GAAP_H', 'FLUX_GAAP_H', 'FLUXERR_GAAP_H', 'GAAP_Flag_H', 'GAAP_nexp_H', 
#	# 'MAG_GAAP_Ks', 'MAGERR_GAAP_Ks', 'FLUX_GAAP_Ks', 'FLUXERR_GAAP_Ks', 'GAAP_Flag_Ks', 'GAAP_nexp_Ks', 
#	# 'EXTINCTION_Z', 'EXTINCTION_Y', 'EXTINCTION_J', 'EXTINCTION_H', 'EXTINCTION_Ks', 
#	# 'MAG_LIM_Z', 'MAG_LIM_Y', 'MAG_LIM_J', 'MAG_LIM_H', 'MAG_LIM_Ks', 
#	# 'GAAP_Flag_u', 'GAAP_Flag_g', 'GAAP_Flag_r', 'GAAP_Flag_i', 
#	# 'Z_B', 'Z_B_MIN', 'Z_B_MAX', 
#	# 'T_B', 'ODDS', 'Z_ML', 'T_ML', 'CHI_SQUARED_BPZ', 'M_0', 
#	# 'BPZ_FILT', 'NBPZ_FILT', 'BPZ_NONDETFILT', 'NBPZ_NONDETFILT', 'BPZ_FLAGFILT', 'NBPZ_FLAGFILT', 
#	# 'GAAP_Flag_ugriZYJHKs', 'wcsx', 'wcsy', 'bias_corrected_e1', 'bias_corrected_e2', 'weight', 'fitclass', 
#	# 'bias_corrected_scalelength_pixels', 'bulge_fraction', 'model_flux', 'pixel_SNratio', 'model_SNratio', 
#	# 'contamination_radius', 'PSF_e1', 'PSF_e2', 'PSF_Strehl_ratio', 'PSF_Q11', 'PSF_Q22', 'PSF_Q12', 
#	# 'star_galaxy_f_probability', 'r_correction', '2D_measurement_variance', 'mean_likelihood_e', 
#	# 'e1_correction', 'e2_correction', 'neighbour_mag', 'neighbour_distance', 
#	# 'catmag', 'n_exposures_used', 'cat_ID', 
#	# 'PSF_e1_exp1', 'PSF_e2_exp1', 'PSF_e1_exp2', 'PSF_e2_exp2', 'PSF_e1_exp3', 'PSF_e2_exp3', 'PSF_e1_exp4', 'PSF_e2_exp4', 'PSF_e1_exp5', 'PSF_e2_exp5', 
#	# 'MASK', 'SG_FLAG', 'MAG_GAAP_u_CALIB', 'MAG_GAAP_g_CALIB', 'MAG_GAAP_r_CALIB', 'MAG_GAAP_i_CALIB', 
#	# 'THELI_NAME', 'THELI_INT', 'SeqNr_field', 'weight_A', 'gmr', 'imr', 'eabs', 
#	# 'binary_cut', 'recal_weight', '2D_measurement_variance_corr']
#
#	##
#	# pdf = pd.DataFrame({'mag_u': data['MAG_GAAP_u_CALIB'], \
#	# 	'mag_g': data['MAG_GAAP_g_CALIB'], \
#	# 	'mag_r':data['MAG_GAAP_r_CALIB'],
#	# 	'mag_i': data['MAG_GAAP_i_CALIB'], \
#	# 	'T_B': data['T_B']})
#
##	mag_u = data['MAG_GAAP_u_CALIB']
#    mag_g = data['MAG_GAAP_g_CALIB']
#    mag_r = data['MAG_GAAP_r_CALIB']
#    mag_i = data['MAG_GAAP_i_CALIB']
#    mag_Z = data['MAG_GAAP_Z']
#    T_B = data['T_B']
#    flag = data['GAAP_Flag_ugriZYJHKs']
#
#
#    magg = []
#    magr = []
#    magi = []
#    magZ = []
#    TB = []
#	
#    for i in range(len(T_B)):
##        if (mag_r[i] < 25) and (mag_r[i] > 20):
#        if (flag[i] == 0) and (mag_g[i]!=99.) and (mag_r[i]!=99.) and (mag_i[i]!=99.) and (mag_Z[i]!=99.):
#            magg.append(mag_g[i])
#            magr.append(mag_r[i])
#            magi.append(mag_i[i])
#            magZ.append(mag_Z[i])
#            TB.append(T_B[i])
#	
#
#    new = pd.DataFrame({'magg':np.array(magg),
#                        'magr':np.array(magr),
#                        'magi':np.array(magi),
#                        'magZ':np.array(magZ),
#                        'TB':np.array(TB)
#                        })
#    new.to_csv("/disks/shear15/ssli/KV450/KV450_G12_light.csv",columns=['magg','magr','magi','magZ','TB'],sep=',',index=False)		





data = pd.read_csv("/disks/shear15/ssli/KV450/KV450_G12_light.csv",sep=',',header=0,index_col=False,\
                   usecols=['magg','magr','magi','magZ','TB'])

#magg = data['magg']
magr = data.loc[data.magZ>-99.,'magr']
#magi = data['magi']
magZ = data.loc[data.magZ>-99.,'magZ']

TB = data.loc[data.magZ>-99.,'TB']




#plt.figure()
#cm = plt.cm.get_cmap('RdYlBu')
#sc = plt.scatter(np.array(magr), np.array(magg)-np.array(magr), c=np.array(TB), vmin=0, vmax=np.amax(np.array(TB)), cmap=cm)
#plt.xlabel("r")
#plt.ylabel("g-r")
#plt.colorbar(sc)
#plt.savefig("/data1/surfdrive_ssli/Projects/6wl_SimIm/plot/r_gr.png")
#plt.close()
#
#
#plt.figure()
#cm = plt.cm.get_cmap('RdYlBu')
#sc = plt.scatter(np.array(magr), np.array(magr)-np.array(magi), c=np.array(TB), vmin=0, vmax=np.amax(np.array(TB)), cmap=cm)
#plt.xlabel("r")
#plt.ylabel("r-i")
#plt.colorbar(sc)
#plt.savefig("/data1/surfdrive_ssli/Projects/6wl_SimIm/plot/r_ri.png")
#plt.close()


plt.figure()
cm = plt.cm.get_cmap('RdYlBu')
sc = plt.scatter(np.array(magr), np.array(magr)-np.array(magZ), c=np.array(TB), vmin=0, vmax=np.amax(np.array(TB)), cmap=cm)
plt.xlabel("r")
plt.ylabel("r-Z")
plt.colorbar(sc)    
plt.savefig("/data1/surfdrive_ssli/Projects/6wl_SimIm/plot/r_rZ.png")
plt.close()










