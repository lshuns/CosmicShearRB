#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 20:39:36 2020

@author: ssli

Running the whole module of cosmic_shear_signal

"""

import os
import math

import initialise
import io_cs
import cosmic_shear


# +++++++++++++++++++++++++++++++++++++++++++ whole data
# path for all necessary input
paths = {}
# cosmological code path (CLASS)
paths['cosmo'] = '/net/eemmeer.strw.leidenuniv.nl/data1/ssli/class_public' 
# data path
paths['data'] = '/disks/shear15/ssli/KV450/KV450_COSMIC_SHEAR_DATA_RELEASE' 
# parameter/configure file path (same as code path)
paths['param'] = os.path.dirname(os.path.realpath(__file__))

# name of parameter/configure files
name_param_file = 'kv450_cf_best.param'
name_conf_file = 'kv450_cf.conf'

# Initialisation
# class: data and cosmo created
cosmo, data = initialise.initialise(paths)
# print(data.paths)


# data filled with input files
# parameter file (with cosmological and nuisance parameters)
data.read_file(name_param_file, 'data', field='', separate=False)
# print(data.cosmo_arguments)
# print(data.nuisance_parameters)
# configure file (with configure and hardly changed setting parameters)
data.read_file(name_conf_file, 'data', field='', separate=False)
# print(data.const)
# print(data.conf)

# create output folder
io_cs.CreateOutpathFunc(data)

# # cosmic shear signal calculation
cosmic_shear.CSsignalFunc(data, cosmo)




