#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 15:49:41 2020

@author: ssli

Running the module of cosmic shear signal prediction

Please modify necessary configurations in Cosmo/cosmic_shear_signal/input

"""



import time

import os
import sys
# Self-defined package
sys.path.insert(0,os.path.realpath('../Cosmo/cosmic_shear_signal')) 
import CosmicShear, initialise


Start = time.time()

# path for all necessary input
paths = {}
# cosmological code path (CLASS)
paths['cosmo'] = '/net/eemmeer/data1/ssli/class_public' 
# father path for all the input and output
paths['data'] = '/disks/shear15/ssli/CosmicShear'
# parameter/configure file path
paths['param'] =  '/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/CosmicShearRB/Cosmo/cosmic_shear_signal/input'
# Initialisation
# class: data and cosmo created
cosmo, data = initialise.initialise(paths)


# ++++++++++++++++++++++++++++++++++++++++++ whole
# name of parameter/configure files
name_param_file = 'kv450_cf_best.param'
name_conf_file = 'kv450_cf.conf'

# data filled with input files
# parameter file (with cosmological and nuisance parameters)
data.read_file(name_param_file, 'data', field='', separate=False)
# configure file (with configure and hardly changed setting parameters)
data.read_file(name_conf_file, 'data', field='', separate=False)

# cosmic shear signal calculation
CosmicShear.CSsignalFunc(data, cosmo)


# ++++++++++++++++++++++++++++++++++++++++++ red
# name of parameter/configure files
name_param_file = 'kv450_cf_best_red.param'
name_conf_file = 'kv450_cf_red.conf'

# data filled with input files
# parameter file (with cosmological and nuisance parameters)
data.read_file(name_param_file, 'data', field='', separate=False)
# configure file (with configure and hardly changed setting parameters)
data.read_file(name_conf_file, 'data', field='', separate=False)

# cosmic shear signal calculation
CosmicShear.CSsignalFunc(data, cosmo)


# ++++++++++++++++++++++++++++++++++++++++++ blue
# name of parameter/configure files
name_param_file = 'kv450_cf_best_blue.param'
name_conf_file = 'kv450_cf_blue.conf'

# data filled with input files
# parameter file (with cosmological and nuisance parameters)
data.read_file(name_param_file, 'data', field='', separate=False)
# configure file (with configure and hardly changed setting parameters)
data.read_file(name_conf_file, 'data', field='', separate=False)

# cosmic shear signal calculation
CosmicShear.CSsignalFunc(data, cosmo)



print("All finished in", time.time()-Start)