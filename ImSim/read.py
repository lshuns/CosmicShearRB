#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:01:10 2019

@author: ssli

load the data from the simulation catalog
set selection criteria 
save for further use
"""

import numpy as np
import pandas as pd

# +++++++++++++++++++++++++ columns in MasterCat
#Index(['oldLFweight', 'LFweight', 'e1', 'e2', 'mag_out', 'g1', 'g2', 'strehl',
#       'snr_model', 'size_out', 'mag_in', 'size_in', 'N_in', 'e1_in', 'e2_in',
#       'ZB4_in', 'Cat_ID', 'size_corr', 'e1_corr', 'e2_corr', 'nd',
#       'LS variance', 'fitclass', 'rotation', 'psf_size_in', 'psf_e1_in',
#       'psf_e2_in', 'X_IMAGE', 'Y_IMAGE', 'prior_matched', 'ZB9_in',
#       'FLUX_RADIUS', 'FWHM_IMAGE', 'bulge_fraction', 'star_gal_prob', 'f_in',
#       'contamination_radius', 'CHI2NU_HI', 'TB9_in', 'e1_true', 'e2_true',
#       'b0r1'],
#      dtype='object')

# +++++++++++++++++++++++++ Main code
indata = "/disks/shear15/ssli/Sim/MasterCat_TSTnewinputglobalRecalbluered_all_5_PSF.csv"
outdata_csv = "/disks/shear15/ssli/Sim/MasterCat_TSTnewinputglobalRecalbluered_all_5_PSF_new.csv"
outdata_hdf = "/disks/shear15/ssli/Sim/MasterCat.h5"
chunk = int(1e6)
# write mode
md = 'w'
# header
head = True

# load initial data
tmp = pd.read_csv(indata,sep=',',header=0,index_col=False,
                  iterator=True,chunksize=chunk)

for data in tmp:
    # criteria
    data_tmp = data[(data['fitclass']==0) | (data['fitclass']==-6) | (data['fitclass']==-9)].copy()
    data_new = data_tmp[(data_tmp['size_out']>0) & (data_tmp['contamination_radius']>4.25) & (data_tmp['prior_matched']==1) & (data_tmp['size_in']>0) & (data_tmp['LFweight']>0)].copy()
    # ellipticity
    x = data_new['e1_in']+data_new['g1']
    y = data_new['e2_in']+data_new['g2']
    k = 1.+data_new['g1']*data_new['e1_in']+data_new['g2']*data_new['e2_in']
    m = data_new['g1']*data_new['e2_in']-data_new['g2']*data_new['e1_in']
    #
    data_new.loc[:,'e1_true'] = (x*k+y*m)/(k*k+m*m)
    data_new.loc[:,'e2_true'] = (y*k-x*m)/(k*k+m*m)
    # split blue and red galaxies
    data_new['b0r1'] = np.where(data_new['TB9_in']<=1.9, 1, 0)
    # save new data
    data_new.to_csv(outdata_csv,sep=',',index=False,mode=md,header=head)
    data_new.to_hdf(outdata_hdf,key='bluered/t0',mode=md,format='table',append=True)
    #
    md = 'a'
    head = False