#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 16:15:25 2019

@author: ssli

make the histograms based on the bins created by hist_hist.py 
"""

import pandas as pd 

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

### info for input bins 
pathbins = "/disks/shear15/ssli/Sim/hist_bins.h5"

### running information
log = open("./log_hist_plot.txt", "a")


### parameters for histograms
para = ['e_true','e_out','ZB9_in','bulge_fraction','N_in','snr_model_log','mag_in','mag_out','TB9_in']

## color
CR = 'red'
CB = 'blue'
CT = 'black'
## transparency
ALP = 1


for s in para:

    bins_T = pd.read_hdf(pathbins,key=s+'/tot/bins')
    counts_T = pd.read_hdf(pathbins,key=s+'/tot/counts')

    bins_R = pd.read_hdf(pathbins,key=s+'/red/bins')
    counts_R = pd.read_hdf(pathbins,key=s+'/red/counts')

    bins_B = pd.read_hdf(pathbins,key=s+'/blue/bins')
    counts_B = pd.read_hdf(pathbins,key=s+'/blue/counts')

    print("Succeed in loading data for",s,file=log)

    ### start plot
    fig = plt.figure()
    ##
    if (s == "mag_out") or (s == "mag_in"):
        plt.hist(bins_T[:-1], bins_T,log=True,weights=counts_T,facecolor=CT,edgecolor=CT,alpha=ALP,histtype='step',label=s+' (total)')
        plt.hist(bins_R[:-1], bins_R,log=True,weights=counts_R,facecolor=CR,edgecolor=CR,alpha=ALP,histtype='step',label=s+' (Red)')
        plt.hist(bins_B[:-1], bins_B,log=True,weights=counts_B,facecolor=CB,edgecolor=CB,alpha=ALP,histtype='step',label=s+' (Blue)')

    else:
        plt.hist(bins_T[:-1], bins_T,weights=counts_T,facecolor=CT,edgecolor=CT,alpha=ALP,histtype='step',label=s+' (total)')
        plt.hist(bins_R[:-1], bins_R,weights=counts_R,facecolor=CR,edgecolor=CR,alpha=ALP,histtype='step',label=s+' (Red)')
        plt.hist(bins_B[:-1], bins_B,weights=counts_B,facecolor=CB,edgecolor=CB,alpha=ALP,histtype='step',label=s+' (Blue)')

    # plt.legend(loc='upper left', frameon=False,fontsize="x-small")
    plt.legend(frameon=False)

    plt.xlabel(s)
    plt.ylabel("Counts")

#    plt.savefig("/data1/surfdrive_ssli/Projects/6wl_SimIm/plot/hist_"+s+".pdf")
    plt.savefig("/data1/surfdrive_ssli/Projects/6wl_SimIm/plot/hist_"+s+".png")    
    plt.close()

    print("Finished in plot for",s,file=log)


