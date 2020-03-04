#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 13:21:00 2020

@author: ssli

Script to run Covariance Prepare
"""

import feather

import os
import sys
# Self-defined package
sys.path.insert(0,os.path.realpath('..')) 
from Covariance import NeffSigmaeFunc

# input directory
indir = "/disks/shear15/ssli/KV450/split/all_tomo"
# input postfix
inP_r = "_T_B_less3"
inP_b = "_T_B_greater3"


area = 341.3 * 3600. # 1/arcmin^2

outdir = "/disks/shear15/ssli/CosmicShear/covariance/prepare/neff_sigmae"


for inP in [inP_r, inP_b]:
    WorA = 'w' 
    for i in range(5):
        inpath = indir + str(i+1) + inP + '.feather'
        indata = feather.read_dataframe(inpath)
        outpath = outdir + inP + '.txt'

        id_zbin = i + 1
        e1 = indata['bias_corrected_e1']
        e2 = indata['bias_corrected_e2']
        wg = indata['recal_weight']

        NeffSigmaeFunc(id_zbin, e1, e2, wg, area, outpath, WorA)
        WorA = 'a'
        print("Finished in", id_zbin, inP)
