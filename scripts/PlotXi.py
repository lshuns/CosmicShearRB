#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 12:59:15 2020

@author: ssli

script to make the data vector plot
"""


import pandas as pd
import numpy as np

import os
import sys
# Self-defined package
sys.path.insert(0,os.path.realpath('..')) 
from Plot import XiPlotFunc


# ++++++++++++++++++++++++++++++++++++++++++ data: whole vs. red vs. blue

# Number of bins 
nzbins = 5

# custom settings for plot
# color: dimgray (solid) / red (half transparent) / blue (half transparent)
CRs = ['dimgray', [1, 0, 0, 0.5], [0, 0, 1, 0.5]]
# marker: circle / diamond / square
MKs = ['o', 'd', 's']
# marker size
MSs = [2, 2, 2]
# linestyle (not used for data)
# LSs = ['none', 'none', 'none']
LSs = ['-', '-.', '--']
# linewidth
LWs = [0.6, 0.6, 0.6]
# linewidth of the errorbar lines
ELWs = [0.6, 0.6, 0.6]
# YTYPE
YTYPE = 'orig'
    
# output directory
outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/data_whole_red_blue.png"


# whole
# data
inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_whole.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_whole = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error
inpath = "/disks/shear15/ssli/CosmicShear/KV450_COSMIC_SHEAR_DATA_RELEASE/COV_MAT/xipmcutcov_KV450_analytic_inc_m.dat"
data = np.loadtxt(inpath)
indx_i = data[:,0]
indx_j = data[:,1]
err_full = data[:,2]
err_diagonal = err_full[indx_j==indx_i]
para_whole['error'] = np.sqrt(err_diagonal)

# red
# data
inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_red.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_red = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error
inpath = "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_red.dat"
data = np.loadtxt(inpath)
err_diagonal = data[:,2]
para_red['error'] = err_diagonal

# blue
# data
inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_blue.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_blue = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error
inpath = "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_blue.dat"
data = np.loadtxt(inpath)
err_diagonal = data[:,2]
para_blue['error'] = err_diagonal

paras = [para_whole, para_red, para_blue]
names = ['data', 'data', 'data']


XiPlotFunc(paras, names, nzbins,
                CRs, MKs, MSs, LSs, LWs, ELWs,
                YTYPE,
                outpath)

# ++++++++++++++++++++++++++++++++++++++++++ data vs. theory: whole

# Number of bins 
nzbins = 5

# custom settings for plot
# color
CRs = ['dimgray', 'black']
# marker: circle
MKs = ['o', None]
# marker size
MSs = [2, None]
# linestyle (not used for data)
LSs = ['none', '-']
# linewidth
LWs = [None, 1.0]
# linewidth of the errorbar lines
ELWs = [1.0, None]
# YTYPE
YTYPE = 'orig'
    
# output directory
outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/theory_data_whole.png"


# whole
# data
inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_whole.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_data = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error
inpath = "/disks/shear15/ssli/CosmicShear/KV450_COSMIC_SHEAR_DATA_RELEASE/COV_MAT/xipmcutcov_KV450_analytic_inc_m.dat"
data = np.loadtxt(inpath)
indx_i = data[:,0]
indx_j = data[:,1]
err_full = data[:,2]
err_diagonal = err_full[indx_j==indx_i]
para_data['error'] = np.sqrt(err_diagonal)

# theory
inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_cut_to_cut_values_5zbins_whole.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_theory = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error (not used)
para_theory['error'] = np.zeros(len(theta))


paras = [para_data, para_theory]
names = ['data', 'theory']

XiPlotFunc(paras, names, nzbins,
                CRs, MKs, MSs, LSs, LWs, ELWs,
                YTYPE,
                outpath)



# ++++++++++++++++++++++++++++++++++++++++++ data vs. theory: red

# Number of bins 
nzbins = 5

# custom settings for plot
# color
CRs = ['red', 'black']
# marker: diamond
MKs = ['d', None]
# marker size
MSs = [2, None]
# linestyle (not used for data)
LSs = ['none', '-']
# linewidth
LWs = [None, 1.0]
# linewidth of the errorbar lines
ELWs = [1.0, None]
# YTYPE
YTYPE = 'orig'
    
# output directory
outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/theory_data_red.png"


# red
# data
inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_red.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_data = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error
inpath = "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_red.dat"
data = np.loadtxt(inpath)
err_diagonal = data[:,2]
para_data['error'] = err_diagonal

# theory
inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_cut_to_cut_values_5zbins_red.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_theory = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error (not used)
para_theory['error'] = np.zeros(len(theta))


paras = [para_data, para_theory]
names = ['data', 'theory']

XiPlotFunc(paras, names, nzbins,
                CRs, MKs, MSs, LSs, LWs, ELWs,
                YTYPE,
                outpath)



# ++++++++++++++++++++++++++++++++++++++++++ data vs. theory: blue

# Number of bins 
nzbins = 5

# custom settings for plot
# color
CRs = ['blue', 'black']
# marker: square
MKs = ['s', None]
# marker size
MSs = [2, None]
# linestyle (not used for data)
LSs = ['none', '-']
# linewidth
LWs = [None, 1.0]
# linewidth of the errorbar lines
ELWs = [1.0, None]
# YTYPE
YTYPE = 'orig'
    
# output directory
outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/theory_data_blue.png"


# blue
# data
inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_blue.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_data = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error
inpath = "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_blue.dat"
data = np.loadtxt(inpath)
err_diagonal = data[:,2]
para_data['error'] = err_diagonal

# theory
inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_cut_to_cut_values_5zbins_blue.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_theory = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error (not used)
para_theory['error'] = np.zeros(len(theta))


paras = [para_data, para_theory]
names = ['data', 'theory']

XiPlotFunc(paras, names, nzbins,
                CRs, MKs, MSs, LSs, LWs, ELWs,
                YTYPE,
                outpath)



# ++++++++++++++++++++++++++++++++++++++++++ data vs. theory: red-blue


# Number of bins 
nzbins = 5

# custom settings for plot
# color
CRs = ['orange', 'black']
# marker: square
MKs = ['o', None]
# marker size
MSs = [2, None]
# linestyle (not used for data)
LSs = ['none', '-']
# linewidth
LWs = [None, 1.0]
# linewidth of the errorbar lines
ELWs = [1.0, None]
# YTYPE
YTYPE = 'diff'
    
# output directory
outpath = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/CorrFunc/theory_data_red_blue_diff.png"


# red
# data
inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_red.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_data_red = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error
inpath = "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_red.dat"
data = np.loadtxt(inpath)
err_diagonal = data[:,2]
para_data_red['error'] = err_diagonal

# theory
inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_cut_to_cut_values_5zbins_red.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_theory_red = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error (not used)
para_theory_red['error'] = np.zeros(len(theta))


# blue
# data
inpath = '/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_withK_blue.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_data_blue = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error
inpath = "/disks/shear15/ssli/CosmicShear/data_vector/for_plot/xi_for_plot_error_blue.dat"
data = np.loadtxt(inpath)
err_diagonal = data[:,2]
para_data_blue['error'] = err_diagonal

# theory
inpath = '/disks/shear15/ssli/CosmicShear/theory_vector/xi_cut_to_cut_values_5zbins_blue.dat'
data = np.loadtxt(inpath)
theta = data[:,1]
xi = data[:,2]
pm = data[:,3]
ito = data[:,4]
jto = data[:,5]
para_theory_blue = pd.DataFrame({'theta': theta, 'xi': xi, 'pm': pm, 'ito': ito, 'jto': jto})
# error (not used)
para_theory_blue['error'] = np.zeros(len(theta))



# difference
# use avarage for theta
para_data = pd.DataFrame({'theta': (para_data_red['theta'].values+para_data_blue['theta'].values)/2.,
                            'xi': para_data_red['xi'].values - para_data_blue['xi'].values,
                            'pm': para_data_red['pm'].values,
                            'ito': para_data_red['ito'].values,
                            'jto': para_data_red['jto'].values,
                            'error': np.sqrt(para_data_red['error'].values**2.+para_data_blue['error'].values**2.)
                        })

para_theory = pd.DataFrame({'theta': (para_theory_red['theta'].values+para_theory_blue['theta'].values)/2.,
                            'xi': para_theory_red['xi'].values - para_theory_blue['xi'].values,
                            'pm': para_theory_red['pm'].values,
                            'ito': para_theory_red['ito'].values,
                            'jto': para_theory_red['jto'].values,
                            'error': para_theory_red['error'].values
                        })


paras = [para_data, para_theory]
names = ['data', 'theory']

XiPlotFunc(paras, names, nzbins,
                CRs, MKs, MSs, LSs, LWs, ELWs,
                YTYPE,
                outpath)
