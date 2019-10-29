#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 09:45:15 2019

@author: ssli

Set new criteria to remove negative TB and ZB
add binning
"""

import pandas as pd

# data
path = "/disks/shear15/ssli/Sim/MasterCat.h5"
store = pd.HDFStore(path, mode='a')

# running information
log = open("/disks/shear15/ssli/Sim/log/log_read2.txt", "a")

# previous data
datain = store.get(key='bluered/t0')
# rename invalid Python identifier
datain.rename(columns={'LS variance': 'LS_variance',}, inplace=True)
store.put(key='bluered/t0', value=datain, format='table', data_columns=True)

print("Columns", file=log)
print(datain.columns, file=log)

print("Information for MasterCat.h5/bluered/t0", file=log)
print("Number of objects:", len(datain), file=log)
print("Number of red:", len(datain[datain.b0r1==1]), file=log)
print("Number of blue:", len(datain[datain.b0r1==0]), file=log)

print("Succeed in loading data.")

# new criteria
dataout = datain[(datain.ZB9_in>0) & (datain.TB9_in>0)]
store.put(key='bluered/t1', value=dataout, format='table', data_columns=True)

print("Information for MasterCat.h5/bluered/t1", file=log)
print("Number of objects:", len(dataout), file=log)
print("Number of red:", len(dataout[dataout.b0r1==1]), file=log)
print("Number of blue:", len(dataout[dataout.b0r1==0]), file=log)

print("Succeed in saving t1.")

# binning
keys = ['z13', 'z35', 'z57', 'z79', 'z912']
zmin = 0.1
dz = 0.2
for i in range(5):
    zmax = zmin + dz

    if i < 4:
        tmp = dataout[(dataout.ZB9_in > zmin) & (dataout.ZB9_in <= zmax)]
    else:
        tmp = dataout[(dataout.ZB9_in > zmin) & (dataout.ZB9_in <= 1.2)]

    store.put(key='bluered/'+keys[i], value=tmp, format='table', data_columns=True)
    
    print("Information for MaterCat.h5/bluered/"+keys[i], file=log)
    print("Number of objects:", len(tmp), file=log)
    print("Number of red:", len(tmp[tmp.b0r1==1]), file=log)
    print("Number of blue:", len(tmp[tmp.b0r1==0]), file=log)
    
    print("Succeed in saving", keys[i])
    
    zmin = zmax

store.close()