#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:14:43 2019

@author: ssli

Module to calculate the data vector (correlation function)

MeanFunc:
    weighted average of e for a given data set (for c-term determination)

CorrFunc:
    Using TreeCorr to calculate the 2-point shear correlation function
    e is modified using weighted e 
    bin_slop is set to be 0.05

CorrCosmoFunc:
    Rearrange the results to be used in cosmo analysis.

CorrErrCosmoFunc:
    Save treecorr reported errors.

CorrPlotFunc:
    Rearrange the results to be used in plot.

CorrErrPlotFunc:
    Rearrange the errors to be used in plot.

"""

from CorrFunc.CorrFunc import MeanFunc, CorrFunc, CorrPlotFunc, CorrCosmoFunc, CorrErrCosmoFunc, CorrErrPlotFunc
