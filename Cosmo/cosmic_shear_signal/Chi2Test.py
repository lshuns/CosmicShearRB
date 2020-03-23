#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 12:51:35 2020

@author: ssli

Module for chi^2 test
"""

import numpy as np
from scipy.linalg import cholesky, solve_triangular

import io_cs
import CosmicShear

OUTPATH = '/disks/shear15/ssli/CosmicShear/data_vector/test/output/'

def Chi2SingleFunc(nzbins, nzcorrs, theta_bins, mask, data, xi_obs, xi_theo):
    """
    Estimate chi^2 for one data vector
    """

    # load the full covariance matrix:
    covmat = io_cs.LoadCovarianceFunc(data, nzbins, nzcorrs, xi_theo)

    # apply mask also to covariance matrix
    mask_indices = np.where(mask == 1)[0]
    covmat = covmat[np.ix_(mask_indices, mask_indices)]

    # precompute Cholesky transform for chi^2 calculation:
    # don't invert that matrix...
    # use the Cholesky decomposition instead:
    cholesky_transform = cholesky(covmat, lower=True)

    vec = xi_theo[mask_indices] - xi_obs[mask_indices]

    yt = solve_triangular(cholesky_transform, vec, lower=True)
    chi2 = yt.dot(yt)

    return chi2, len(vec)



# cannot work due to cholesky problem
def Chi2CoupleFunc(nzbins, nzcorrs, theta_bins, mask1, mask2, data1, data2, xi_obs_1, xi_obs_2, xi_theo_1, xi_theo_2):
    """
    Estimate chi^2 for coupled data vectors
    Note: this assumes two data vectors have two separated covariance matrices
        the cross-correlation between two data vectors is ignored
    """

    # load the full covariance matrix:
    covmat_block_1 = io_cs.LoadCovarianceFunc(data1, nzbins, nzcorrs, xi_theo_1)
    covmat_block_2 = io_cs.LoadCovarianceFunc(data2, nzbins, nzcorrs, xi_theo_2)

    covmat_block_zero = np.ones_like(covmat_block_1)

    # build a combined cov-mat, for that to work we assume, that the cov-mat dimension fits
    # to the size of the *uncut*, single data-vector and is ordered in the same way as the
    # *final* data-vector created (i.e. vec = [xi+(1,1), xi-(1,1), xi+(1,2), xi-(1,2),...]!
    covmat = np.asarray(np.bmat('covmat_block_1, covmat_block_zero; covmat_block_2, covmat_block_zero'))

    # trim covariance matrix to chosen scales:
    mask = np.concatenate((mask1, mask2))
    mask_indices = np.where(mask == 1)[0]
    covmat = covmat[np.ix_(mask_indices, mask_indices)]

    # precompute Cholesky transform for chi^2 calculation:
    # don't invert that matrix...
    # use the Cholesky decomposition instead:
    cholesky_transform = cholesky(covmat, lower=True)

    vec = np.concatenate((xi_theo_1, xi_theo_2))[mask_indices] - np.concatenate((xi_obs_1, xi_obs_2))[mask_indices]
    
    yt = solve_triangular(cholesky_transform, vec, lower=True)
    chi2 = yt.dot(yt)

    return chi2, len(vec)


