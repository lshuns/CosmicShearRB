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

def Chi2SingleFunc(nzbins, nzcorrs, theta_bins, mask, data, xi_obs, xi_theo, DISCARD_POINT=False):
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

    # discard one point at a time
    if DISCARD_POINT:
        chi2s = []
        for i in range(len(vec)):
            vec_dis = np.copy(vec)
            vec_dis[i] = 0
            yt = solve_triangular(cholesky_transform, vec_dis, lower=True)
            chi2s.append(yt.dot(yt))
        return np.array(chi2s)
    yt = solve_triangular(cholesky_transform, vec, lower=True)
    chi2 = yt.dot(yt)

    return chi2, len(vec)



def Chi2CoupleDiffFunc(nzbins, nzcorrs, ntheta, mask,
                        data1, xi_obs_1, xi_theo_1,
                        data2, xi_obs_2, xi_theo_2,
                        inDir_cov12, file_name_cov12):
    """
    Estimate chi^2 for difference between two data vectors
    Note: this assumes two data vectors have two separated covariance matrices
        the cross-correlation between two data vectors is also desired
        the masks for two data vector need to be identical
    """

    # load the full covariance matrix:
    covmat_block_1 = io_cs.LoadCovarianceFunc(data1, nzbins, nzcorrs, xi_theo_1)
    covmat_block_2 = io_cs.LoadCovarianceFunc(data2, nzbins, nzcorrs, xi_theo_2)

    covmat_block_12 = io_cs.LoadCrossCovarianceFunc(inDir_cov12, file_name_cov12, ntheta, nzbins, nzcorrs, xi_theo_1, xi_theo_2)
    

    # build a combined cov-mat
    covmat = covmat_block_1 + covmat_block_2 - covmat_block_12 - covmat_block_12.transpose()

    # trim covariance matrix to chosen scales:
    mask_indices = np.where(mask == 1)[0]
    covmat = covmat[np.ix_(mask_indices, mask_indices)]

    # precompute Cholesky transform for chi^2 calculation:
    # don't invert that matrix...
    # use the Cholesky decomposition instead:
    cholesky_transform = cholesky(covmat, lower=True)

    vec = (xi_theo_1[mask_indices] - xi_obs_1[mask_indices]) - (xi_theo_2[mask_indices] - xi_obs_2[mask_indices])
    
    yt = solve_triangular(cholesky_transform, vec, lower=True)
    chi2 = yt.dot(yt)

    return chi2, len(vec)


