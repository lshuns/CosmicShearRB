#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:46:16 2020

@author: ssli

Running the module of chi^2 test

Please modify necessary configurations in Cosmo/cosmic_shear_signal/input

"""

import time

import os
import sys
# Self-defined package
sys.path.insert(0, os.path.realpath('../Cosmo/cosmic_shear_signal')) 
import CosmicShear, initialise, Chi2Test


Start = time.time()

# path for all necessary input
paths = {}
# cosmological code path (CLASS)
paths['cosmo'] = '/net/eemmeer/data1/ssli/class_public' 
# father path for all the input and output
paths['data'] = '/disks/shear15/ssli/CosmicShear'
# parameter/configure file path
paths['param'] =  '/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/CosmicShearRB/Cosmo/cosmic_shear_signal/input'


# ++++++++++++++++++++++++++++++++++++++++++ whole
# name of parameter/configure files
name_param_file = 'kv450_cf_best.param'
name_conf_file = 'kv450_cf.conf'

# Initialisation
# class: data and cosmo created
cosmo, data_whole = initialise.initialise(paths)
# data filled with input files
# parameter file (with cosmological and nuisance parameters)
data_whole.read_file(name_param_file, 'data', field='', separate=False)
# configure file (with configure and hardly changed setting parameters)
data_whole.read_file(name_conf_file, 'data', field='', separate=False)

# cosmic shear signal calculation
xi_obs_whole, xi_theo_whole, theta_bins_whole, mask_whole = CosmicShear.CSsignalFunc(data_whole, cosmo, False)

# chi2 
# Number of bins
nzbins = len(data_whole.const['z_bins_min'])
# Number of correlation
nzcorrs = int(nzbins * (nzbins + 1) / 2)
#
chi2_whole, dof_whole = Chi2Test.Chi2SingleFunc(nzbins, nzcorrs, theta_bins_whole, mask_whole, data_whole, xi_obs_whole, xi_theo_whole)
print("chi2_whole", chi2_whole)
print("dof_whole", dof_whole)

# ++++++++++++++++++++++++++++++++++++++++++ red
# name of parameter/configure files
name_param_file = 'kv450_cf_best.param'
name_conf_file = 'kv450_cf_red.conf'

# Initialisation
# class: data and cosmo created
cosmo, data_red = initialise.initialise(paths)
# data filled with input files
# parameter file (with cosmological and nuisance parameters)
data_red.read_file(name_param_file, 'data', field='', separate=False)
# configure file (with configure and hardly changed setting parameters)
data_red.read_file(name_conf_file, 'data', field='', separate=False)

# cosmic shear signal calculation
xi_obs_red, xi_theo_red, theta_bins_red, mask_red = CosmicShear.CSsignalFunc(data_red, cosmo, False)
# xi_obs_red, xi_theo_red, theta_bins_red, mask_red = CosmicShear.CSsignalFunc(data_red, cosmo, True)

# chi2 
# Number of bins
nzbins = len(data_red.const['z_bins_min'])
# Number of correlation
nzcorrs = int(nzbins * (nzbins + 1) / 2)
#
chi2_red, dof_red = Chi2Test.Chi2SingleFunc(nzbins, nzcorrs, theta_bins_red, mask_red, data_red, xi_obs_red, xi_theo_red)
print("chi2_red", chi2_red)
print("dof_red", dof_red)

# ++++++++++++++++++++++++++++++++++++++++++ blue
# name of parameter/configure files
name_param_file = 'kv450_cf_best.param'
name_conf_file = 'kv450_cf_blue.conf'

# Initialisation
# class: data and cosmo created
cosmo, data_blue = initialise.initialise(paths)
# data filled with input files
# parameter file (with cosmological and nuisance parameters)
data_blue.read_file(name_param_file, 'data', field='', separate=False)
# configure file (with configure and hardly changed setting parameters)
data_blue.read_file(name_conf_file, 'data', field='', separate=False)

# cosmic shear signal calculation
xi_obs_blue, xi_theo_blue, theta_bins_blue, mask_blue = CosmicShear.CSsignalFunc(data_blue, cosmo, False)
# xi_obs_blue, xi_theo_blue, theta_bins_blue, mask_blue = CosmicShear.CSsignalFunc(data_blue, cosmo, True)

# chi2 
# Number of bins
nzbins = len(data_blue.const['z_bins_min'])
# Number of correlation
nzcorrs = int(nzbins * (nzbins + 1) / 2)
#
chi2_blue, dof_blue = Chi2Test.Chi2SingleFunc(nzbins, nzcorrs, theta_bins_blue, mask_blue, data_blue, xi_obs_blue, xi_theo_blue)
print("chi2_blue", chi2_blue)
print("dof_blue", dof_blue)


# # # ++++++++++++++++++++++++++++++++++++++++++++ red - blue
# # chi2_rb, dof_rb = Chi2Test.Chi2CoupleFunc(nzbins, nzcorrs, theta_bins_red, mask_red, mask_blue, data_red, data_blue, xi_obs_red, xi_obs_blue, xi_theo_red, xi_theo_blue)
# # print("chi2_red_blue", chi2_rb)
# # print("dof_red_blue", dof_rb)
# # print("chi2_reduced_red_blue", chi2_rb/(dof_rb-1))


# print("All finished in", time.time()-Start)
# # ('All finished in', 25.485830068588257) (covariance matrix in list form)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ old
################# results from mean values
# # whole
# Intitial cosmological parameters passed to CLASS code:
# {'omega_cdm': 0.118, 'ln10^{10}A_s': 3.158, 'omega_b': 0.022, 'n_s': 1.021, 'h': 0.745, 'non linear': 'hmcode', 'c_min': 2.484, 'Omega_k': 0.0, 'N_eff': 2.0328, 'N_ncdm': 1.0, 'm_ncdm': 0.06, 'T_ncdm': 0.71611, 'output': 'mPk', 'P_k_max_h/Mpc': 170.0, 'sBBN file': '/net/eemmeer/data1/ssli/class_public/bbn/sBBN.dat', 'k_pivot': 0.05, 'z_max_pk': 5.9750000000000005}
# sigma8 = 0.8905375160205157
# Omega_m = 0.25340144300336764
# h = 0.745

# 'chi2_whole', 187.248500267432
# 'dof_whole', 195

# #red
#
# A_IA = 0.981 
# ('chi2_red', 158.91163818760813)
# ('dof_red', 195)
#
# A_IA = 0
# chi2_red 159.10040204122924
# dof_red 195

# whole matrix
# A_IA = 0.981 
# ('chi2_red', 776.1795011987977
# ('dof_red', 195)
#
# A_IA = 0
# chi2_red 775.7116916240524
# dof_red 195

# # blue
#
# A_IA = 0.981
# ('chi2_blue', 186.727356105233)
# ('dof_blue', 195)
#
# A_IA = 0
# chi2_blue 200.51255150505304
# dof_blue 195
# All finished in 19.25131869316101

# whole matrix
# A_IA = 0.981
# ('chi2_blue', 586.8404108038512
# ('dof_blue', 195)
#
# A_IA = 0
# chi2_blue 610.9348264268585
# dof_blue 195
# All finished in 19.25131869316101


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ new
################# results from max-likelihood values
# Intitial cosmological parameters passed to CLASS code:
# {'omega_cdm': 0.05748817, 'ln10^{10}A_s': 4.697331, 'omega_b': 0.0219527, 'n_s': 1.128071, 'h': 0.7801301, 'non linear': 'hmcode', 'c_min': 2.188542, 'Omega_k': 0.0, 'N_eff': 2.0328, 'N_ncdm': 1.0, 'm_ncdm': 0.06, 'T_ncdm': 0.71611, 'output': 'mPk', 'P_k_max_h/Mpc': 170.0, 'sBBN file': '/net/eemmeer/data1/ssli/class_public/bbn/sBBN.dat', 'k_pivot': 0.05, 'z_max_pk': 5.9750000000000005}
# sigma8 = 1.1647685473207734
# Omega_m = 0.13158825944920308
# h = 0.7801301
# S8 = 0.7714141104481363

# ########### whole
# +++++ old covariance matrix (old theory vector)
# chi2_whole 181.289777435608
# dof_whole 195
# +++++ new covariance matrix (new theory vector)
# chi2_whole 181.28974941896647
# dof_whole 195
# +++++ IA = 0
# chi2_whole 182.18528673318514
# dof_whole 195

# ########### red
# +++++ IA = best value
# chi2_red 159.05609058659545
# dof_red 195
# +++++ IA = 0
# chi2_red 157.9054310819387
# dof_red 195


# ########### blue
# +++++ IA = best value
# chi2_blue 176.58063045837295
# dof_blue 195
# +++++ IA = 0
# chi2_blue 179.45973994197797
# dof_blue 195

