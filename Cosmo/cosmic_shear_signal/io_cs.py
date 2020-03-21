#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 14:22:36 2020

@author: ssli

Input-Output handling
"""

from __future__ import print_function
import os
import sys
import numpy as np


def ReadDataVectorFunc(data, nzbins, nzcorrs):
    """
    Read data vector and bring it into the desired format 
    """
    data_path = data.paths['data']
    data_sample = data.conf['sample']

    ntheta = data.const['ntheta']

    # plus one for theta-column
    data_xip = np.zeros((ntheta, nzcorrs + 1))
    data_xim = np.zeros((ntheta, nzcorrs + 1))
    idx_corr = 0
    for zbin1 in range(nzbins):
        for zbin2 in range(zbin1, nzbins):

            data_vector_path = os.path.join(data_path, 'data_vector/for_cosmo/{:}/xi_for_cosmo_tomo_{:}_{:}_withK_{:}.dat'.format(data_sample, zbin1+1, zbin2+1, data_sample))
            theta, xip, xim = np.loadtxt(data_vector_path, unpack=True)

            # this assumes theta is the same for every tomographic bin and
            # for both xi_p and xi_m!
            if idx_corr == 0:
                data_xip[:, 0] = theta
                data_xim[:, 0] = theta

            data_xip[:, idx_corr + 1] = xip
            data_xim[:, idx_corr + 1] = xim

            idx_corr += 1

    data = np.concatenate((data_xip, data_xim))

    print('Loaded data vectors from: \n', os.path.join(data_path, 'data_vector/for_cosmo/{:}'.format(data_sample)), '\n')

    return data


def ReadCutValueFunc(data, nzbins, nzcorrs, theta_bins):
    """
    Read cut values and convert into mask
    """
    data_path = data.paths['data']
    cutvalues_file = data.conf['cutvalues_file']

    ntheta = data.const['ntheta']

    cutvalues_file_path = os.path.join(data_path, 'SUPPLEMENTARY_FILES/CUT_VALUES/'+cutvalues_file)
    if os.path.exists(cutvalues_file_path):
        cut_values = np.loadtxt(cutvalues_file_path)
    else:
        raise Exception('File not found:\n {:} \n \
            Check that requested file exists in the following folder: \
            \n {:}'.format(cutvalues_file_path))

    # create the mask
    mask = np.zeros(2 * nzcorrs * ntheta)
    iz = 0
    for izl in range(nzbins):
        for izh in range(izl, nzbins):
            # this counts the bin combinations
            # iz=1 =>(1,1), iz=2 =>(1,2) etc
            iz = iz + 1
            for i in range(ntheta):
                j = (iz-1)*2*ntheta
                xi_plus_cut_low = max(cut_values[izl, 0], cut_values[izh, 0])
                xi_plus_cut_high = max(cut_values[izl, 1], cut_values[izh, 1])
                xi_minus_cut_low = max(cut_values[izl, 2], cut_values[izh, 2])
                xi_minus_cut_high = max(cut_values[izl, 3], cut_values[izh, 3])
                if ((theta_bins[i] < xi_plus_cut_high) and (theta_bins[i]>xi_plus_cut_low)):
                    mask[j+i] = 1
                if ((theta_bins[i] < xi_minus_cut_high) and (theta_bins[i]>xi_minus_cut_low)):
                    mask[ntheta + j+i] = 1

    mask_indices = np.where(mask == 1)[0]
    mask_suffix = cutvalues_file[:-4]

    return mask, mask_indices, mask_suffix


def ReadZdistribution(data, nzbins):
    """
    Read the dn/dz files
    """
    z_bins_min = data.const['z_bins_min']
    z_bins_max = data.const['z_bins_max']

    data_path = data.paths['data']
    sample_name = data.conf['sample']
    nz_method = data.conf['nz_method']

    # Create labels for loading of dn/dz-files:
    zbin_labels = []
    for i in range(nzbins):
        zbin_labels += ['{:.1f}t{:.1f}'.format(z_bins_min[i], z_bins_max[i])]

    # Read fiducial dn_dz from window files:
    z_samples = []
    hist_samples = []
    for zbin in range(nzbins):

        window_file_path = os.path.join(
            data_path, 'redshift/' + sample_name + '/Nz_{0:}/Nz_{0:}_Mean/Nz_{0:}_z{1:}.asc'.format(nz_method, zbin_labels[zbin]))

        if os.path.exists(window_file_path):
            zptemp, hist_pz = np.loadtxt(window_file_path, usecols=[0, 1], unpack=True)
            shift_to_midpoint = np.diff(zptemp)[0] / 2.
            if zbin == 0:
                zpcheck = zptemp
            if np.sum((zptemp - zpcheck)**2) > 1e-6:
                raise Exception('The redshift values for the window files at different bins do not match.')

            # we add a zero as first element because we want to integrate down to z = 0!
            z_samples += [np.concatenate((np.zeros(1), zptemp + shift_to_midpoint))]
            hist_samples += [np.concatenate((np.zeros(1), hist_pz))]

        else:
            raise Exception("dn/dz file not found:\n %s"%window_file_path)
    print('Loaded redshift distributions from: \n', os.path.join(
            data_path, 'redshift/' + sample_name + '/Nz_{0:}/Nz_{0:}_Mean/'.format(nz_method)), '\n')

    return np.asarray(z_samples), np.asarray(hist_samples), len(zptemp)


def WriteVectorFunc(data, nzbins, theta_bins, mask_indices, mask_suffix, vec, fname_suffix):
    """
    Write out the cosmic shear vector in list-format with mask
    """
    ntheta = data.const['ntheta']

    data_path = data.paths['data']
    out_folder = data.conf['out_folder']
    cutvalues_file = data.conf['cutvalues_file']

    # Initialise all output colunms
    thetas_all = []
    idx_pm = []
    idx_tomo_z1 = []
    idx_tomo_z2 = []
    for idx_z1 in range(nzbins):
        for idx_z2 in range(idx_z1, nzbins):
            thetas_all.append(theta_bins)
            idx_pm.append(np.concatenate((np.ones(ntheta), np.ones(ntheta) + 1)))
            idx_tomo_z1.append(np.ones(2 * ntheta) * (idx_z1 + 1))
            idx_tomo_z2.append(np.ones(2 * ntheta) * (idx_z2 + 1))
    thetas_all = np.concatenate(thetas_all)
    idx_pm = np.concatenate(idx_pm)
    idx_tomo_z1 = np.concatenate(idx_tomo_z1)
    idx_tomo_z2 = np.concatenate(idx_tomo_z2)

    # now apply correct masking:
    thetas_all = thetas_all[mask_indices]
    idx_pm = idx_pm[mask_indices]
    idx_tomo_z1 = idx_tomo_z1[mask_indices]
    idx_tomo_z2 = idx_tomo_z2[mask_indices]
    vec = vec[mask_indices]

    idx_run = np.arange(len(idx_pm)) + 1
    savedata = np.column_stack((idx_run, thetas_all, vec, idx_pm, idx_tomo_z1, idx_tomo_z2))
    header = ' i    theta(i)\'        xi_p/m(i)  (p=1, m=2)  itomo   jtomo'
    fname = os.path.join(data_path, '{:}/xi_cut_to_{:}_{:}.dat'.format(out_folder, mask_suffix, fname_suffix))
    np.savetxt(fname, savedata, header=header, delimiter='\t', fmt=['%4i', '%.5e', '%12.5e', '%i', '%i', '%i'])
    print('Saved vector in list format cut down to scales as specified in {:}: \n'.format(cutvalues_file), fname, '\n')


def LoadCovarianceFunc(data, nzbins, nzcorrs, xi_theo):
    """
    Read in the full covariance matrix and to bring it into format of self.xi_obs.
    """

    data_path = data.paths['data']
    list_file = data.conf['list_covariance']
    usable_file = data.conf['usable_covariance']

    ntheta = data.const['ntheta']

    try:
        fname = os.path.join(data_path, 'covariance/'+usable_file)
        matrix = np.loadtxt(fname)
        print('Loaded covariance matrix (incl. shear calibration uncertainty) in a format usable from: \n', fname, '\n')

    except:
        fname = os.path.join(data_path, 'covariance/'+list_file)
        tmp_raw = np.loadtxt(fname)

        print('Loaded covariance matrix in list format from: \n', fname)
        print('Now we construct the covariance matrix in a format usable for the first time.')

        indices = np.column_stack((tmp_raw[:, :3], tmp_raw[:, 4:7]))

        # we need to add both components for full covariance
        values = tmp_raw[:, 8] + tmp_raw[:, 9]

        dim = 2 * ntheta * nzcorrs
        matrix = np.zeros((dim, dim))

        # make sure the list covariance is in right order:
        # theta -> PorM -> iz2 -> iz1
        index_lin = 0
        # this creates the correctly ordered (i.e. like self.xi_obs) full
        # 270 x 270 covariance matrix:
        for index1 in range(dim):
            for index2 in range(index1, dim):
                matrix[index1, index2] = values[index_lin]
                matrix[index2, index1] = matrix[index1, index2]
                index_lin += 1

        # apply propagation of m-correction uncertainty following
        # equation 12 from Hildebrandt et al. 2017 (arXiv:1606.05338):
        err_multiplicative_bias = 0.02
        matrix_m_corr = np.matrix(xi_theo).T * np.matrix(xi_theo) * 4. * err_multiplicative_bias**2
        matrix = matrix + np.asarray(matrix_m_corr)

        fname = fname = os.path.join(data_path, 'covariance/'+usable_file)
        if not os.path.isfile(fname):
            np.savetxt(fname, matrix)
            print('Saved covariance matrix (incl. shear calibration uncertainty) in format usable with this likelihood to: \n', fname, '\n')

    return matrix