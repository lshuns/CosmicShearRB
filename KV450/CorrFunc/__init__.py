#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 14:14:43 2019

@author: ssli
"""

# standard modules
import treecorr
import feather

import numpy as np
import pandas as pd
import multiprocessing as mp
import matplotlib as mpl
import matplotlib.pyplot as plt

from astropy.io import fits


# defined modules
from A_A_DataPreprocess import SelecFunc, BinFunc, BinCountFunc
from A_CorrFunc import MeanFunc, CorrFunc, CorrRearrangeFunc
from B_plot import xiPlotFunc
