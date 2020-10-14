#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 12:59:15 2020

@author: ssli

script to make the data vector plot
"""


import pandas as pd
import numpy as np

from XiPlot import XiPlotFunc


# # ++++++++++++++++++++++++++++++++++++++++++ data: whole vs. red vs. blue

# # Number of bins 
# nzbins = 5

# # custom settings for plot
# # color: dimgray (solid) / red (half transparent) / blue (half transparent)
# CRs = ['dimgray', [1, 0, 0, 0.5], [0, 0, 1, 0.5]]
# # marker: circle / diamond / square
# MKs = ['o', 'd', 's']
# # marker size
# MSs = [2, 2, 2]
# # linestyle (not used for data)
# # LSs = ['none', 'none', 'none']
# LSs = ['-', '-.', '--']
# # linewidth
# LWs = [0.6, 0.6, 0.6]
# # linewidth of the errorbar lines
# ELWs = [0.6, 0.6, 0.6]
# # YTYPE
# YTYPE = 'orig'
    
# # output directory
# outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/data_whole_red_blue.png"


# # whole
# # data
# inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_whole.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_whole = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error
# inpath = "/disks/shear15/ssli/CosmicShear/KV450_COSMIC_SHEAR_DATA_RELEASE/COV_MAT/xipmcutcov_KV450_analytic_inc_m.dat"
# data = np.loadtxt(inpath)
# indx_i = data[:,0]
# indx_j = data[:,1]
# err_full = data[:,2]
# err_diagonal = err_full[indx_j==indx_i]
# para_whole['error'] = np.sqrt(err_diagonal)

# # red
# # data
# inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_red.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_red = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error
# inpath = "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_red.dat"
# data = np.loadtxt(inpath)
# err_diagonal = data[:,2]
# para_red['error'] = err_diagonal

# # blue
# # data
# inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_blue.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_blue = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error
# inpath = "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_blue.dat"
# data = np.loadtxt(inpath)
# err_diagonal = data[:,2]
# para_blue['error'] = err_diagonal

# paras = [para_whole, para_red, para_blue]
# names = ['data', 'data', 'data']


# XiPlotFunc(paras, names, nzbins,
#                 CRs, MKs, MSs, LSs, LWs, ELWs,
#                 YTYPE,
#                 outpath)

# # ++++++++++++++++++++++++++++++++++++++++++ data vs. theory: whole

# # Number of bins 
# nzbins = 5

# # custom settings for plot
# # color
# CRs = ['dimgray', 'black']
# # marker: circle
# MKs = ['o', None]
# # marker size
# MSs = [2, None]
# # linestyle (not used for data)
# LSs = ['none', '-']
# # linewidth
# LWs = [None, 1.0]
# # linewidth of the errorbar lines
# ELWs = [1.0, None]
# # YTYPE
# YTYPE = 'orig'
    
# # output directory
# outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/theory_data_whole.png"


# # whole
# # data
# inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_whole.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_data = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error
# inpath = "/disks/shear15/ssli/CosmicShear/KV450_COSMIC_SHEAR_DATA_RELEASE/COV_MAT/xipmcutcov_KV450_analytic_inc_m.dat"
# data = np.loadtxt(inpath)
# indx_i = data[:,0]
# indx_j = data[:,1]
# err_full = data[:,2]
# err_diagonal = err_full[indx_j==indx_i]
# para_data['error'] = np.sqrt(err_diagonal)

# # theory
# inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_cut_to_cut_values_5zbins_whole.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_theory = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error (not used)
# para_theory['error'] = np.zeros(len(theta))


# paras = [para_data, para_theory]
# names = ['data', 'theory']

# XiPlotFunc(paras, names, nzbins,
#                 CRs, MKs, MSs, LSs, LWs, ELWs,
#                 YTYPE,
#                 outpath)



# # ++++++++++++++++++++++++++++++++++++++++++ data vs. theory: red

# # Number of bins 
# nzbins = 5

# # custom settings for plot
# # color
# CRs = ['red', 'black']
# # marker: diamond
# MKs = ['d', None]
# # marker size
# MSs = [2, None]
# # linestyle (not used for data)
# LSs = ['none', '-']
# # linewidth
# LWs = [None, 1.0]
# # linewidth of the errorbar lines
# ELWs = [1.0, None]
# # YTYPE
# YTYPE = 'orig'
    
# # output directory
# outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/theory_data_red.png"


# # red
# # data
# inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_red.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_data = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error
# inpath = "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_red.dat"
# data = np.loadtxt(inpath)
# err_diagonal = data[:,2]
# para_data['error'] = err_diagonal

# # theory
# inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_cut_to_cut_values_5zbins_red.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_theory = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error (not used)
# para_theory['error'] = np.zeros(len(theta))


# paras = [para_data, para_theory]
# names = ['data', 'theory']

# XiPlotFunc(paras, names, nzbins,
#                 CRs, MKs, MSs, LSs, LWs, ELWs,
#                 YTYPE,
#                 outpath)



# # ++++++++++++++++++++++++++++++++++++++++++ data vs. theory: blue

# # Number of bins 
# nzbins = 5

# # custom settings for plot
# # color
# CRs = ['blue', 'black']
# # marker: square
# MKs = ['s', None]
# # marker size
# MSs = [2, None]
# # linestyle (not used for data)
# LSs = ['none', '-']
# # linewidth
# LWs = [None, 1.0]
# # linewidth of the errorbar lines
# ELWs = [1.0, None]
# # YTYPE
# YTYPE = 'orig'
    
# # output directory
# outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/theory_data_blue.png"


# # blue
# # data
# inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_blue.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_data = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error
# inpath = "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_blue.dat"
# data = np.loadtxt(inpath)
# err_diagonal = data[:,2]
# para_data['error'] = err_diagonal

# # theory
# inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_cut_to_cut_values_5zbins_blue.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_theory = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error (not used)
# para_theory['error'] = np.zeros(len(theta))


# paras = [para_data, para_theory]
# names = ['data', 'theory']

# XiPlotFunc(paras, names, nzbins,
#                 CRs, MKs, MSs, LSs, LWs, ELWs,
#                 YTYPE,
#                 outpath)



# # ++++++++++++++++++++++++++++++++++++++++++ data vs. theory (KV450 only): red-blue


# # Number of bins 
# nzbins = 5

# # custom settings for plot
# # color
# CRs = ['orange', 'black']
# # marker: square
# MKs = ['o', None]
# # marker size
# MSs = [2, None]
# # linestyle (not used for data)
# LSs = ['none', '-']
# # linewidth
# LWs = [None, 1.0]
# # linewidth of the errorbar lines
# ELWs = [1.0, None]
# # YTYPE
# YTYPE = 'diff'
    
# # output directory
# outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/theory_data_red_blue_diff.png"


# # red
# # data
# inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_red.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_data_red = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error
# inpath = "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_red.dat"
# data = np.loadtxt(inpath)
# err_diagonal = data[:,2]
# para_data_red['error'] = err_diagonal

# # theory
# inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_cut_to_cut_values_5zbins_red.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_theory_red = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error (not used)
# para_theory_red['error'] = np.zeros(len(theta))


# # blue
# # data
# inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_blue.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_data_blue = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error
# inpath = "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_blue.dat"
# data = np.loadtxt(inpath)
# err_diagonal = data[:,2]
# para_data_blue['error'] = err_diagonal

# # theory
# inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_cut_to_cut_values_5zbins_blue.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_theory_blue = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error (not used)
# para_theory_blue['error'] = np.zeros(len(theta))



# # difference
# # use avarage for theta
# para_data = pd.DataFrame({'theta': (para_data_red['theta'].values+para_data_blue['theta'].values)/2.,
#                             'xi': para_data_red['xi'].values - para_data_blue['xi'].values,
#                             'pm': para_data_red['pm'].values,
#                             'ito': para_data_red['ito'].values,
#                             'jto': para_data_red['jto'].values,
#                             'error': np.sqrt(para_data_red['error'].values**2.+para_data_blue['error'].values**2.)
#                         })

# para_theory = pd.DataFrame({'theta': (para_theory_red['theta'].values+para_theory_blue['theta'].values)/2.,
#                             'xi': para_theory_red['xi'].values - para_theory_blue['xi'].values,
#                             'pm': para_theory_red['pm'].values,
#                             'ito': para_theory_red['ito'].values,
#                             'jto': para_theory_red['jto'].values,
#                             'error': para_theory_red['error'].values
#                         })


# paras = [para_data, para_theory]
# names = ['data', 'theory']

# XiPlotFunc(paras, names, nzbins,
#                 CRs, MKs, MSs, LSs, LWs, ELWs,
#                 YTYPE,
#                 outpath)



# # ++++++++++++++++++++++++++++++++++++++++++ data (only): blue-red


# # Number of bins 
# nzbins = 5

# # custom settings for plot
# # color 
# # without K / with K
# CRs = ['gray', 'orange']
# # marker: square
# MKs = ['o', 'o']
# # marker size
# MSs = [2, 2]
# # linestyle (not used for data)
# LSs = ['none', 'none']
# # linewidth
# LWs = [None, None]
# # linewidth of the errorbar lines
# ELWs = [1.0, 1.0]
# # YTYPE
# YTYPE = 'diff'
    
# # output directory
# outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/data2_br_diff.png"


# ########## without K
# # red
# # data
# inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withoutK_less3.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_data_red = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error
# inpath = "/disks/shear15/ssli/CosmicShear/covariance/apr8/thps_cov_apr8_rr_inc_m_usable_mask.dat"
# err_r = np.loadtxt(inpath)


# # blue
# # data
# inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withoutK_greater3.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_data_blue = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# # error
# inpath = "/disks/shear15/ssli/CosmicShear/covariance/apr8/thps_cov_apr8_bb_inc_m_usable_mask.dat"
# err_b = np.loadtxt(inpath)

# # difference
# # cross error
# inpath = "/disks/shear15/ssli/CosmicShear/covariance/apr8/thps_cov_apr8_br_inc_m_usable_mask.dat"
# err_c = np.loadtxt(inpath)
# err = err_b + err_r - 2.*err_c
# err_diagonal = np.diag(err)
# # use avarage for theta
# para_data1 = pd.DataFrame({'theta': (para_data_red['theta'].values+para_data_blue['theta'].values)/2.,
#                             'xi': para_data_blue['xi'].values - para_data_red['xi'].values,
#                             'pm': para_data_red['pm'].values,
#                             'ito': para_data_red['ito'].values,
#                             'jto': para_data_red['jto'].values,
#                             'error': np.sqrt(err_diagonal)
#                         })


# ########## with K
# # red
# # data
# inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_less3.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_data_red = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})


# # blue
# # data
# inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_greater3.dat'
# data = np.loadtxt(inpath)
# theta = data[:,1]
# xi = data[:,2]
# pm = data[:,3]
# ito = data[:,4]
# jto = data[:,5]
# para_data_blue = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})

# # difference
# # use avarage for theta
# para_data2 = pd.DataFrame({'theta': (para_data_red['theta'].values+para_data_blue['theta'].values)/2.,
#                             'xi': para_data_blue['xi'].values - para_data_red['xi'].values,
#                             'pm': para_data_red['pm'].values,
#                             'ito': para_data_red['ito'].values,
#                             'jto': para_data_red['jto'].values,
#                             'error': np.sqrt(err_diagonal)
#                             })

# paras = [para_data1, para_data2]
# names = ['data', 'data']

# XiPlotFunc(paras, names, nzbins,
#                 CRs, MKs, MSs, LSs, LWs, ELWs,
#                 YTYPE,
#                 outpath)



# ++++++++++++++++++++++++++++++++++++++++++ Final plot: data vs. theory (KV450 & Planck): blue-red

# Number of bins 
nzbins = 5
    
# output directory
outpath1 = "/disks/shear15/ssli/surfdrive/Projects/6CosmicShear_RB/plot/publish/data_br_diff_KV450_Planck2.pdf"
outpath2 = "/disks/shear15/ssli/surfdrive/Projects/6CosmicShear_RB/plot/CorrFunc/data_br_diff_KV450_Planck2.png"


# red
# data
inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_less3.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_data_red = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error
inpath = "/disks/shear15/ssli/CosmicShear/covariance/apr8/thps_cov_apr8_rr_inc_m_usable_mask.dat"
err_r = np.loadtxt(inpath)

# KV450
inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_theory_cut_less3_KV450_best.dat'

data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_KV450_red = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})

# Planck
inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_theory_cut_less3_Planck.dat'

data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_Planck_red = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})


# blue
# data
inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_greater3.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_data_blue = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error
inpath = "/disks/shear15/ssli/CosmicShear/covariance/apr8/thps_cov_apr8_bb_inc_m_usable_mask.dat"
err_b = np.loadtxt(inpath)

# KV450
inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_theory_cut_greater3_KV450_best.dat'

data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_KV450_blue = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})

# Planck
inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_theory_cut_greater3_Planck.dat'

data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_Planck_blue = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})

# difference
# cross error
inpath = "/disks/shear15/ssli/CosmicShear/covariance/apr8/thps_cov_apr8_br_inc_m_usable_mask.dat"
err_c = np.loadtxt(inpath)
err = err_b + err_r - 2.*err_c
err_diagonal = np.diag(err)
# use avarage for theta
para_data = pd.DataFrame({'theta': (para_data_red['theta'].values+para_data_blue['theta'].values)/2.,
                            'xi': para_data_blue['xi'].values - para_data_red['xi'].values,
                            'pm': para_data_red['pm'].values,
                            'ito': para_data_red['ito'].values,
                            'jto': para_data_red['jto'].values,
                            'error': np.sqrt(err_diagonal)
                        })

para_KV450 = pd.DataFrame({'theta': (para_KV450_red['theta'].values+para_KV450_blue['theta'].values)/2.,
                            'xi': para_KV450_blue['xi'].values - para_KV450_red['xi'].values,
                            'pm': para_KV450_red['pm'].values,
                            'ito': para_KV450_red['ito'].values,
                            'jto': para_KV450_red['jto'].values
                        })

para_Planck = pd.DataFrame({'theta': (para_Planck_red['theta'].values+para_Planck_blue['theta'].values)/2.,
                            'xi': para_Planck_blue['xi'].values - para_Planck_red['xi'].values,
                            'pm': para_Planck_red['pm'].values,
                            'ito': para_Planck_red['ito'].values,
                            'jto': para_Planck_red['jto'].values
                        })

########## without K
# red
# data
inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withoutK_less3.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_data_red = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error
inpath = "/disks/shear15/ssli/CosmicShear/covariance/apr8/thps_cov_apr8_rr_inc_m_usable_mask.dat"
err_r = np.loadtxt(inpath)


# blue
# data
inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withoutK_greater3.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_data_blue = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error
inpath = "/disks/shear15/ssli/CosmicShear/covariance/apr8/thps_cov_apr8_bb_inc_m_usable_mask.dat"
err_b = np.loadtxt(inpath)

# difference
# cross error
inpath = "/disks/shear15/ssli/CosmicShear/covariance/apr8/thps_cov_apr8_br_inc_m_usable_mask.dat"
err_c = np.loadtxt(inpath)
err = err_b + err_r - 2.*err_c
err_diagonal = np.diag(err)
# use avarage for theta
para_data1 = pd.DataFrame({'theta': (para_data_red['theta'].values+para_data_blue['theta'].values)/2.,
                            'xi': para_data_blue['xi'].values - para_data_red['xi'].values,
                            'pm': para_data_red['pm'].values,
                            'ito': para_data_red['ito'].values,
                            'jto': para_data_red['jto'].values,
                            'error': np.sqrt(err_diagonal)
                        })

# custom settings for plot
paras = [para_data1, para_data, para_KV450, para_Planck]

### data KV450 Planck
names = ['data', 'data', 'theory', 'theory']
LABELS = ['raw measurement', 'with shear calibration', 'KV450 cosmology', 'Planck cosmology']
# color
CRs = ['gray', 'orange', 'red', 'black']
# marker: square
MKs = ['o', 'x', None, None]
# marker size
MSs = [2, 2, None, None]
# linestyle (not used for data)
LSs = ['none', 'none', '-', '--']
# linewidth
LWs = [None, None, 1.0, 1.0]
# linewidth of the errorbar lines
ELWs = [1.0, 1.0, None, None]
# YTYPE
YTYPE = 'diff'

XiPlotFunc(paras, names, nzbins,
                CRs, MKs, MSs, LSs, LWs, ELWs,
                YTYPE,
                outpath1,
                LABELS)
XiPlotFunc(paras, names, nzbins,
                CRs, MKs, MSs, LSs, LWs, ELWs,
                YTYPE,
                outpath2,
                LABELS)
