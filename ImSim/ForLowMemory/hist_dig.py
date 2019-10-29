#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 13:48:05 2019

@author: ssli

make the histograms using chunk-by-chunk,
    which is friendly to low memory machines
(Based on digitize, weights are infeasible)

(worn out)
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
## (parameters for histograms)
pathpara = "/disks/shear15/ssli/Sim/hist_para.h5"
hdfpara = pd.HDFStore(pathpara)

### running information
log = open("./log.txt", "a")

### parameters for histograms
para = ['e_true','e_out','ZB9_in','bulge_fraction','N_in','snr_model','mag_in','mag_out','TB9_in']

### Number of bins
NUM_BINS = 100

### parameters for histograms
for s in para:
    ### load the data
    if s == 'e_true':
        chunks = pd.read_hdf(pathin,key='bluered/t0',\
            columns = ['e1_true','e2_true','b0r1'],\
                iterator=True,chunksize=chunksize)
    elif s == 'e_out':
        chunks = pd.read_hdf(pathin,key='bluered/t0',\
            columns = ['e1','e2','b0r1'],\
                iterator=True,chunksize=chunksize)
    else:
        chunks = pd.read_hdf(pathin,key='bluered/t0',\
            columns = [s,'b0r1'],\
                iterator=True,chunksize=chunksize)
    
    ##
    print("Succeed in loading chunks for",s,file=log)
    
    ### first loop
    FIRST = True
    FIRST_R = True
    FIRST_B = True
    ##
    IP = 1
    for data in chunks:
        if s == 'e_true':
            tmp_R = np.sqrt(data.loc[data.b0r1==1,'e1_true'].values**2.\
                +data.loc[data.b0r1==1,'e2_true'].values**2.)
            tmp_B = np.sqrt(data.loc[data.b0r1==0,'e1_true'].values**2.\
                +data.loc[data.b0r1==0,'e2_true'].values**2.)
                
        elif s == 'e_out':
            tmp_R = np.sqrt(data.loc[data.b0r1==1,'e1'].values**2.\
                +data.loc[data.b0r1==1,'e2'].values**2.)
            tmp_B = np.sqrt(data.loc[data.b0r1==0,'e1'].values**2.\
                +data.loc[data.b0r1==0,'e2'].values**2.)
        else:
            tmp_R = data.loc[data.b0r1==1,s].values
            tmp_B = data.loc[data.b0r1==0,s].values

        if FIRST:
            ### first loop
            hdfpara.put(key=s+'/red',value=pd.Series(tmp_R),format='table')
            hdfpara.put(key=s+'/blue',value=pd.Series(tmp_B),format='table')
            #
            FIRST = False
        else:
            hdfpara.append(key=s+'/red',value=pd.Series(tmp_R))
            hdfpara.append(key=s+'/blue',value=pd.Series(tmp_B))
            

        ### build the histograms
        # red
        if tmp_R.size != 0:
            if FIRST_R:
                ## first loop
                MIN_BINS_R = np.amin(tmp_R)
                MAX_BINS_R = np.amax(tmp_R)
                #
                bins_R, db_R = np.linspace(MIN_BINS_R,MAX_BINS_R,NUM_BINS,retstep=True)
                bins_R = np.append(bins_R,MAX_BINS_R+db_R)
                ### hist [min,max)
                hist_R = np.bincount(np.digitize(tmp_R,bins_R))[1:]
            
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

                    if MAX_BINS_R_N >= (MAX_BINS_R+N*db_R):
                        bins_R = np.append(bins_R,(MAX_BINS_R+(N+1)*db_R))
                        hist_R = np.append(hist_R,0)
                #
                hist_R_N = np.bincount(np.digitize(tmp_R,bins_R))[1:]
                if len(hist_R_N) != len(hist_R):
                    N = len(hist_R)-len(hist_R_N)
                    hist_R_N = np.pad(hist_R_N,(0,N),'constant',constant_values=0.)
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
                bins_B = np.append(bins_B, MAX_BINS_B+db_B)
                hist_B = np.bincount(np.digitize(tmp_B,bins_B))[1:]

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
                    if MAX_BINS_B_N >= (MAX_BINS_B+N*db_B):
                        bins_B = np.append(bins_B,(MAX_BINS_B+(N+1)*db_B))
                        hist_B = np.append(hist_B,0)
                    
                #
                hist_B_N = np.bincount(np.digitize(tmp_B,bins_B))[1:]
                if len(hist_B_N) != len(hist_B):
                    N = len(hist_B)-len(hist_B_N)
                    hist_B_N = np.pad(hist_B_N,(0,N),'constant',constant_values=0.)
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
    hdfbins.put(key=s+'/red/bins',value=pd.Series(bins_R),format='fixed')
    hdfbins.put(key=s+'/blue/bins',value=pd.Series(bins_B),format='fixed')
    #
    hdfbins.put(key=s+'/red/counts',value=pd.Series(hist_R),format='fixed')
    hdfbins.put(key=s+'/blue/counts',value=pd.Series(hist_B),format='fixed')
    
    ##
    print("Data Saved for",s,file=log)
    

hdfbins.close()
hdfpara.close()
