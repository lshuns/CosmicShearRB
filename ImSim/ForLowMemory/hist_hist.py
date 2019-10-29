#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 14:23:17 2019

@author: ssli

make the histograms using chunk-by-chunk,
    which is friendly to low memory machines
(Based on np.histogram, weights are feasible)

"""


import numpy as np
import pandas as pd
import math

### info for input data
pathin = "/disks/shear15/ssli/Sim/MasterCat.h5"
chunksize = int(1e6)

### output the results 
## (bins for histograms)
pathbins = "/disks/shear15/ssli/Sim/hist_bins.h5"
hdfbins = pd.HDFStore(pathbins)

### running information
log = open("./log_hist.txt", "a")

### parameters for histograms
para = ['e_true','e_out','ZB9_in','bulge_fraction','N_in','snr_model_log','mag_in','mag_out','TB9_in']

### Number of bins
NUM_BINS = 100

### parameters for histograms
for s in para:
    ### load the data
    if s == 'e_true':
        chunks = pd.read_hdf(pathin,key='bluered/t0',\
            columns = ['e1_true','e2_true','LFweight','b0r1'],\
                iterator=True,chunksize=chunksize)
    elif s == 'e_out':
        chunks = pd.read_hdf(pathin,key='bluered/t0',\
            columns = ['e1','e2','LFweight','b0r1'],\
                iterator=True,chunksize=chunksize)
    elif s == 'snr_model_log':
        chunks = pd.read_hdf(pathin,key='bluered/t0',\
            columns = ['snr_model','LFweight','b0r1'],\
                iterator=True,chunksize=chunksize)
    else:
        chunks = pd.read_hdf(pathin,key='bluered/t0',\
            columns = [s,'LFweight','b0r1'],\
                iterator=True,chunksize=chunksize)
    
    ##
    print("Succeed in loading chunks for",s,file=log)
    
    ### first loop
    FIRST_T = True
    FIRST_R = True
    FIRST_B = True
    ##
    IP = 1
    for data in chunks:
        if s == 'e_true':
            tmp = np.sqrt(data['e1_true'].values**2.\
                +data['e2_true'].values**2.)
            tmp_R = np.sqrt(data.loc[data.b0r1==1,'e1_true'].values**2.\
                +data.loc[data.b0r1==1,'e2_true'].values**2.)
            tmp_B = np.sqrt(data.loc[data.b0r1==0,'e1_true'].values**2.\
                +data.loc[data.b0r1==0,'e2_true'].values**2.)
            # weight for hist
            wg = data['LFweight'].values
            wg_R = data.loc[data.b0r1==1,'LFweight'].values
            wg_B = data.loc[data.b0r1==0,'LFweight'].values
                
        elif s == 'e_out':
            tmp = np.sqrt(data['e1'].values**2.\
                +data['e2'].values**2.)
            tmp_R = np.sqrt(data.loc[data.b0r1==1,'e1'].values**2.\
                +data.loc[data.b0r1==1,'e2'].values**2.)
            tmp_B = np.sqrt(data.loc[data.b0r1==0,'e1'].values**2.\
                +data.loc[data.b0r1==0,'e2'].values**2.)
            # weight for hist
            wg = data['LFweight'].values
            wg_R = data.loc[data.b0r1==1,'LFweight'].values
            wg_B = data.loc[data.b0r1==0,'LFweight'].values
        
        elif s == 'snr_model_log':
            tmp = np.log(data.loc[data.snr_model>0,"snr_model"].values)
            tmp_R = np.log(data.loc[(data.snr_model>0) & (data.b0r1==1),"snr_model"].values)
            tmp_B = np.log(data.loc[(data.snr_model>0) & (data.b0r1==0),"snr_model"].values)
            # weight for hist
            wg = data.loc[data.snr_model>0,'LFweight'].values
            wg_R = data.loc[(data.snr_model>0) & (data.b0r1==1),'LFweight'].values
            wg_B = data.loc[(data.snr_model>0) & (data.b0r1==0),'LFweight'].values
            
        else:
            tmp = data[s].values
            tmp_R = data.loc[data.b0r1==1,s].values
            tmp_B = data.loc[data.b0r1==0,s].values
            # weight for hist
            wg = data['LFweight'].values
            wg_R = data.loc[data.b0r1==1,'LFweight'].values
            wg_B = data.loc[data.b0r1==0,'LFweight'].values

        ### build the histograms
        # total
        if tmp.size != 0:                            
            if FIRST_T:
                ## first loop
                MIN_BINS = np.amin(tmp)
                MAX_BINS = np.amax(tmp)
                # 
                bins, db = np.linspace(MIN_BINS,MAX_BINS,NUM_BINS,retstep=True)
                ### hist (All but the last (righthand-most) bin is half-open. )
                hist, bins = np.histogram(tmp, bins=bins, weights=wg, density=True)
            
                ### new bounds
                MIN_BINS = bins[0]
                MAX_BINS = bins[-1]
            
                FIRST_T = False
            
            else:
                MIN_BINS_N = np.amin(tmp)
                MAX_BINS_N = np.amax(tmp)
                if MIN_BINS_N < MIN_BINS:
                    N = math.ceil((MIN_BINS-MIN_BINS_N)/db)
                    ## new part of bins
                    bins_sig = np.linspace(MIN_BINS-N*db,MIN_BINS-db,N)
                    bins = np.concatenate((bins_sig,bins))
                    ## pad the original histogram matrices 
                    hist = np.pad(hist,(N,0),'constant',constant_values=0.)
                if MAX_BINS_N >= MAX_BINS:
                    N = math.ceil((MAX_BINS_N-MAX_BINS)/db)
                    ## new part of bins
                    bins_sig = np.linspace(MAX_BINS+db,MAX_BINS+N*db,N)
                    bins = np.concatenate((bins,bins_sig))
                    ## pad the original histogram matrices 
                    hist = np.pad(hist,(0,N),'constant',constant_values=0.)

                #
                hist_N, bins = np.histogram(tmp, bins=bins, weights=wg, density=True)
                #
                hist += hist_N
         
                ### new bounds
                MIN_BINS = bins[0]
                MAX_BINS = bins[-1]

        # red
        if tmp_R.size != 0:
            if FIRST_R:
                ## first loop
                MIN_BINS_R = np.amin(tmp_R)
                MAX_BINS_R = np.amax(tmp_R)
                #
                bins_R, db_R = np.linspace(MIN_BINS_R,MAX_BINS_R,NUM_BINS,retstep=True)
                ### hist (All but the last (righthand-most) bin is half-open. )
                hist_R, bins_R = np.histogram(tmp_R, bins=bins_R, weights=wg_R, density=True)
            
                ### new bounds
                MIN_BINS_R = bins_R[0]
                MAX_BINS_R = bins_R[-1]
            
                FIRST_R = False
            
            else:
                MIN_BINS_R_N = np.amin(tmp_R)
                MAX_BINS_R_N = np.amax(tmp_R)
                if MIN_BINS_R_N < MIN_BINS_R:
                    N = math.ceil((MIN_BINS_R-MIN_BINS_R_N)/db_R)
                    ## new part of bins
                    bins_sig = np.linspace(MIN_BINS_R-N*db_R,MIN_BINS_R-db_R,N)
                    bins_R = np.concatenate((bins_sig,bins_R))
                    ## pad the original histogram matrices 
                    hist_R = np.pad(hist_R,(N,0),'constant',constant_values=0.)
                if MAX_BINS_R_N >= MAX_BINS_R:
                    N = math.ceil((MAX_BINS_R_N-MAX_BINS_R)/db_R)
                    ## new part of bins
                    bins_sig = np.linspace(MAX_BINS_R+db_R,MAX_BINS_R+N*db_R,N)
                    bins_R = np.concatenate((bins_R,bins_sig))
                    ## pad the original histogram matrices 
                    hist_R = np.pad(hist_R,(0,N),'constant',constant_values=0.)

                #
                hist_R_N, bins_R = np.histogram(tmp_R, bins=bins_R, weights=wg_R, density=True)
            
                #
                hist_R += hist_R_N
            
                ### new bounds
                MIN_BINS_R = bins_R[0]
                MAX_BINS_R = bins_R[-1]
            
        # blue 
        if tmp_B.size != 0:
                
            if FIRST_B:
                ## first loop
                MIN_BINS_B = np.amin(tmp_B)
                MAX_BINS_B = np.amax(tmp_B)
                #
                bins_B, db_B = np.linspace(MIN_BINS_B,MAX_BINS_B,NUM_BINS,retstep=True)
                #
                hist_B, bins_B = np.histogram(tmp_B, bins=bins_B, weights=wg_B, density=True)
            
                ### new bounds
                MIN_BINS_B = bins_B[0]
                MAX_BINS_B = bins_B[-1]
  
                FIRST_B = False
            
            else:
                MIN_BINS_B_N = np.amin(tmp_B)
                MAX_BINS_B_N = np.amax(tmp_B)
                if MIN_BINS_B_N < MIN_BINS_B:
                    N = math.ceil((MIN_BINS_B-MIN_BINS_B_N)/db_B)
                    ## new part of bins
                    bins_sig = np.linspace(MIN_BINS_B-N*db_B,MIN_BINS_B-db_B,N)
                    bins_B = np.concatenate((bins_sig,bins_B))
                    ## pad the original histogram matrices 
                    hist_B = np.pad(hist_B,(N,0),'constant',constant_values=0.)
                                
                if MAX_BINS_B_N >= MAX_BINS_B:
                    N = math.ceil((MAX_BINS_B_N-MAX_BINS_B)/db_B)
                    ## new part of bins
                    bins_sig = np.linspace(MAX_BINS_B+db_B,MAX_BINS_B+N*db_B,N)
                    bins_B = np.concatenate((bins_B,bins_sig))
                    ## pad the original histogram matrices 
                    hist_B = np.pad(hist_B,(0,N),'constant',constant_values=0.)
                    
                #
                hist_B_N, bins_B = np.histogram(tmp_B, bins=bins_B, weights=wg_B, density=True)

                #
                hist_B += hist_B_N
            
                ### new bounds
                MIN_BINS_B = bins_B[0]
                MAX_BINS_B = bins_B[-1]
        ##
        print("Succeed in loop",IP,"for",s,file=log)
        IP += 1
    ##
    print("Finished in",s,file=log)
        

    ### save the results
    hdfbins.put(key=s+'/tot/bins',value=pd.Series(bins),format='fixed')
    hdfbins.put(key=s+'/red/bins',value=pd.Series(bins_R),format='fixed')
    hdfbins.put(key=s+'/blue/bins',value=pd.Series(bins_B),format='fixed')
    #
    hdfbins.put(key=s+'/tot/counts',value=pd.Series(hist),format='fixed')
    hdfbins.put(key=s+'/red/counts',value=pd.Series(hist_R),format='fixed')
    hdfbins.put(key=s+'/blue/counts',value=pd.Series(hist_B),format='fixed')
    
    ##
    print("Data Saved for",s,file=log)
    

hdfbins.close()
