#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:46:16 2020

@author: ssli

Running the module of chi^2 test

Please modify necessary configurations in Cosmo/cosmic_shear_signal/input

"""

import time

import numpy as np

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


# name_param_file = 'kv450_cf_best.param'
# name_param_file = 'kv450_cf_null_nuisance.param'
# name_param_file = 'kv450_cf_mean.param'
name_param_file = 'kv450_cf_Planck.param'

# ++++++++++++++++++++++++++++++++++++++++++ whole
# name of parameter/configure files
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
# print("chi2_whole", chi2_whole)
# print("dof_whole", dof_whole)

# ++++++++++++++++++++++++++++++++++++++++++ red
# name of parameter/configure files
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
# print("chi2_red", chi2_red)
# print("dof_red", dof_red)

# ++++++++++++++++++++++++++++++++++++++++++ blue
# name of parameter/configure files
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
# print("chi2_blue", chi2_blue)
# print("dof_blue", dof_blue)


# ++++++++++++++++++++++++++++++++++++++++++++ blue - red
inDir_cov12 = '/disks/shear15/ssli/CosmicShear/covariance'
# file_name_cov12 = 'thps_cov_mar11_list_br.dat'
file_name_cov12 = 'thps_cov_mar11_usable_br.dat'
ntheta = 9
# check the mask is the same
if np.array_equal(mask_red, mask_blue):
    chi2_br, dof_br = Chi2Test.Chi2CoupleDiffFunc(nzbins, nzcorrs, ntheta, mask_red,
                            data_blue, xi_obs_blue, xi_theo_blue,
                            data_red, xi_obs_red, xi_theo_red,
                            inDir_cov12, file_name_cov12)
    # print("chi2_blue_red", chi2_br)
    # print("dof_blue_red", dof_br)
else:
    print("Something wrong with the mask!")


print("chi2_whole", chi2_whole)
print("chi2_red", chi2_red)
print("chi2_blue", chi2_blue)
print("chi2_blue_red", chi2_br)
print("All finished in", time.time()-Start)
# ('All finished in', 25.485830068588257) (covariance matrix in list form)


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 2020-04-02
# All finished in 24.952291250228882
# ########### whole (old covariance)
# +++ best KV450
# chi2_whole 181.289777435608
# +++ IA = 0
# chi2_whole 182.14159048043757
# +++ null nuisance
# chi2_whole 183.03711897015557
# +++ mean KV450
# chi2_whole 186.6469440889119
# +++ Planck
# chi2_whole 197.15775145270865

# ########### red
# +++ best KV450
# chi2_red 159.0857498091302
# +++ IA = 0
# chi2_red 157.91943258570063
# +++ null nuisance
# chi2_red 157.49834365474632
# +++ mean KV450
# chi2_red 158.74952682259186
# +++ Planck
# chi2_red 163.18518748760374

# ########### blue
# +++ best KV450
# chi2_blue 176.69554322681768
# +++ IA = 0
# chi2_blue 179.6342164659231
# +++ null nuisance
# chi2_blue 181.4358222272628
# +++ mean KV450
# chi2_blue 187.62772845975311
# +++ Planck
# chi2_blue 202.74030188552254

# ########### blue-red
# +++ best KV450
# chi2_blue_red 176.6212374720392
# +++ IA = 0
# chi2_blue_red 177.24543612145942
# +++ null nuisance
# chi2_blue_red 177.47529222732436
# +++ mean KV450
# chi2_blue_red 178.08327855721166
# +++ Planck
# chi2_blue_red 181.72856546228772



# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ old 2
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
# +++++ mean KV450 values
# chi2_whole 186.6469440889119
# +++++ Planck values
# chi2_whole 197.15775145270865

# ########### red
# +++++ IA = best value
# chi2_red 159.05609058659545
# dof_red 195
# +++++ IA = 0
# chi2_red 157.9054310819387
# dof_red 195
# +++++ mean KV450 values
# chi2_red 158.70514127432648
# +++++ Planck values
# chi2_red 162.94534120660316

# ########### blue
# +++++ IA = best value
# chi2_blue 176.58063045837295
# dof_blue 195
# +++++ IA = 0
# chi2_blue 179.45973994197797
# dof_blue 195
# +++++ mean KV450 values
# chi2_blue 186.7343894369092
# +++++ Planck values
# chi2_blue 200.86952530128653

# ########## blue-red
# +++++ IA = best value
# chi2_blue_red 175.57548773807528
# dof_blue_red 195
# +++++ IA = 0
# chi2_blue_red 176.13607482251632
# dof_blue_red 195
# +++++ mean KV450 values
# chi2_blue_red 176.81964188989852
# +++++ Planck values
# chi2_blue_red 179.9253175222906


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

