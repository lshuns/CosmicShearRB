#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 13:50:48 2019

@author: ssli

make the histogram plot
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.use('Agg')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

# color
CT = 'black'
CR = 'red'
CB = 'blue'
# bins
NB = 100

# data
path = "/disks/shear15/ssli/Sim/MasterCat.h5"
store = pd.HDFStore(path, mode='r')
#print(store.info())
print(store.keys())

print("Succeed in loading data")

# datasets for histograms
setsG = '/bluered/'
setsL = ['t0', 't1', 'z13', 'z35', 'z57', 'z79', 'z912']

for s in setsL:
    key = setsG + s
    data = store.select(key=key, columns=['e1_true', 'e2_true', 'e1', 'e2', 
                                          'LFweight', 'b0r1'])
    e_true = np.sqrt(data['e1_true'].values**2. + data['e2_true'].values**2.)
    e_out = np.sqrt(data['e1'].values**2. + data['e2'].values**2.)
    wg = data['LFweight'].values
    
    # red
    e_true_R = np.sqrt(data.loc[data.b0r1==1, 'e1_true'].values**2.\
                       + data.loc[data.b0r1==1, 'e2_true'].values**2.)
    e_out_R = np.sqrt(data.loc[data.b0r1==1, 'e1'].values**2.\
                      + data.loc[data.b0r1==1, 'e2'].values**2.)
    wg_R = data.loc[data.b0r1==1, 'LFweight'].values
    # blue
    e_true_B = np.sqrt(data.loc[data.b0r1==0, 'e1_true'].values**2.\
                       + data.loc[data.b0r1==0, 'e2_true'].values**2.)
    e_out_B = np.sqrt(data.loc[data.b0r1==0, 'e1'].values**2.\
                      + data.loc[data.b0r1==0, 'e2'].values**2.)
    wg_B = data.loc[data.b0r1==0, 'LFweight'].values
    
    # plot
    # e_true
    fig = plt.figure()
    plt.hist(x=e_true, bins=NB, density=True, weights=wg, color=CT, label='Total', histtype='step')
    plt.hist(x=e_true_R, bins=NB, density=True, weights=wg_R, color=CR, label='Red', histtype='step')
    plt.hist(x=e_true_B, bins=NB, density=True, weights=wg_B, color=CB, label='Blue', histtype='step')
    plt.legend(frameon=False)
    plt.xlabel('e_true')
    plt.ylabel('Probability')
    plt.savefig("/net/raam/data1/surfdrive_ssli/Projects/6wl_SimIm/plot/hist_e_true_"+s+'.png')
    plt.close()

    print("Succeed in plot of e_true for", s)

    # e_out
    fig = plt.figure()
    plt.hist(x=e_out, bins=NB, density=True, weights=wg, color=CT, label='Total', histtype='step')
    plt.hist(x=e_out_R, bins=NB, density=True, weights=wg_R, color=CR, label='Red', histtype='step')
    plt.hist(x=e_out_B, bins=NB, density=True, weights=wg_B, color=CB, label='Blue', histtype='step')
    plt.legend(frameon=False)
    plt.xlabel('e_out')
    plt.ylabel('Probability')
    plt.savefig("/net/raam/data1/surfdrive_ssli/Projects/6wl_SimIm/plot/hist_e_out_"+s+'.png')
    plt.close()

    print("Succeed in plot of e_out for", s)

store.close()