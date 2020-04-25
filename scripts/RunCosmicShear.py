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


# +++++++++++++++++++ General setting
# parameter files
name_param_files = ['kv450_cf_best.param', 'kv450_cf_Planck.param']

# Set the output folder
out_folder = 'theory_vector'
out_suffix_s = ['KV450_best', 'Planck']


# path for all necessary input
paths = {}
# cosmological code path (CLASS)
paths['cosmo'] = '/net/eemmeer/data1/ssli/class_public' 
# father path for all the input and output
paths['data'] = '/disks/shear15/ssli/CosmicShear'
# parameter/configure file path
paths['param'] =  '/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/CosmicShearRB/Cosmo/cosmic_shear_signal/input'

# +++++++++++++++++++ Running scripts
Start = time.time()

# Initialisation
# class: data and cosmo created
cosmo, data = initialise.initialise(paths)

for i in range(len(name_param_files)):

    # parameter file (with cosmological and nuisance parameters)
    name_param_file = name_param_files[i]
    data.read_file(name_param_file, 'data', field='', separate=False)
    # output folder
    data.conf['out_folder'] = out_folder
    data.conf['out_suffix'] = out_suffix_s[i]

    # # ++++++++++++++++++++++++++++++++++++++++++ whole
    # # name of configure files
    # name_conf_file = 'kv450_cf.conf'

    # # data filled with input files
    # # configure file (with configure and hardly changed setting parameters)
    # data.read_file(name_conf_file, 'data', field='', separate=False)

    # # cosmic shear signal calculation
    # CosmicShear.CSsignalFunc(data, cosmo, save_theory_vector=True)

    # ++++++++++++++++++++++++++++++++++++++++++ red
    # name of parameter/configure files
    name_conf_file = 'kv450_cf_red.conf'

    # data filled with input files
    # configure file (with configure and hardly changed setting parameters)
    data.read_file(name_conf_file, 'data', field='', separate=False)

    # cosmic shear signal calculation
    CosmicShear.CSsignalFunc(data, cosmo, save_theory_vector=True)


    # ++++++++++++++++++++++++++++++++++++++++++ blue
    # name of parameter/configure files
    name_conf_file = 'kv450_cf_blue.conf'

    # data filled with input files
    # configure file (with configure and hardly changed setting parameters)
    data.read_file(name_conf_file, 'data', field='', separate=False)

    # cosmic shear signal calculation
    CosmicShear.CSsignalFunc(data, cosmo, save_theory_vector=True)

print("All finished in", time.time()-Start)
# eemmeer (2020-04-23)
# ('All finished in', 23.075681924819946)