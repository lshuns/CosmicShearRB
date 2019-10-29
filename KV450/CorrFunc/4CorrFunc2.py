#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 12:03:16 2019

@author: ssli

Rearrange the correlation function from CorrFunc.py
Make it comparable to the results from Hildebrandt 2018 
"""

import numpy as np
import pandas as pd

# input directory
# for old result from 3oldCorrFunc.py
pathinF = "/disks/shear15/ssli/KV450/CS/mine/TreeCorr/full/mine/old/" 
# for new result from 3newCorrFunc.py
# pathinF = "/disks/shear15/ssli/KV450/CS/mine/TreeCorr/full/mine/new/" 

pathinP = ".dat"

# output directory
# for old
pathout = "/disks/shear15/ssli/KV450/CS/mine/TreeCorr/xi_mine_old.csv"
# for new
# pathout = "/disks/shear15/ssli/KV450/CS/mine/TreeCorr/xi_mine_new.csv"

res = pd.DataFrame(columns=['r_nom', 'meanr', 'meanlogr', 'xi_p/m (real)', 'xi_p/m (abs)', 'sigma_p/m', '(p=1, m=2)', 'itomo', 'jtomo'])
for i in range(5):
    j = i + 1
    k = j
    
    while k <= 5:
        pathin = pathinF + str(j) + str(k) + pathinP 
        dat = np.loadtxt(pathin)
        
        # xi_p
        tmp = pd.DataFrame({'r_nom': dat[0:7,0],
                            'meanr': dat[0:7,1],
                            'meanlogr': dat[0:7,2],
                            'xi_p/m (real)': dat[0:7,3],
                            'xi_p/m (abs)': np.sqrt(dat[0:7,3]**2. + dat[0:7,5]**2.),
                            'sigma_p/m': dat[0:7,7],
                            '(p=1, m=2)': np.full((7,),1),
                            'itomo': np.full((7,),1) * j,
                            'jtomo': np.full((7,),1) * k
                            
                })
        res = pd.concat([res, tmp], ignore_index=True)
        
        # xi_m
        tmp = pd.DataFrame({'r_nom': dat[-6:,0],
                            'meanr': dat[-6:,1],
                            'meanlogr': dat[-6:,2],
                            'xi_p/m (real)': dat[-6:,4],
                            'xi_p/m (abs)': np.sqrt(dat[-6:,4]**2. + dat[-6:,6]**2.),
                            'sigma_p/m': dat[-6:,8],
                            '(p=1, m=2)': np.full((6,),2),
                            'itomo': np.full((6,),1) * j,
                            'jtomo': np.full((6,),1) * k
                })
        res = pd.concat([res, tmp], ignore_index=True)

        k += 1

res.to_csv(path_or_buf=pathout,
           header=True, index=False)