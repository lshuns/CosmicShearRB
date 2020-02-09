#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 19:24:06 2020

@author: ssli

Module dealing with the simulated catalogues

SelecFunc:
    Apply mask and re-save the original catalogues.

TomoBinFunc:
    Build and save the tomographic catalogues.

    The fiducial KV450 cosmic shear bins:
    Bin1: (0.1, 0.3]
    Bin2: (0.3, 0.5]
    Bin3: (0.5, 0.7]
    Bin4: (0.7, 0.9]
    Bin5: (0.9, 1.2]

"""

from SimCat.DataProcess import SelecFunc, TomoBinFunc
