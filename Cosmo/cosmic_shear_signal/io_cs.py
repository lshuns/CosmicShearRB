"""
Input-Output handling

Handles all the input/output of the code (at least most of it). If something is
printed that does not satisfy you (number of decimals, for instance, in the
output files), you only have to find the called function and change a number.

Whenever the arguments of the functions are :code:`command_line` or
:code:`data`, no mention of them will be done - as it is now clear. On the
contrary, if there are more arguments, they will be detailled.

This module also defines a new class :class:`File`, that extends
:py:class:`file`, which provides a tail function. It is used in
:func:`sampler.read_args_from_chain`.

Finally, the way the error messages are displayed is set there, along with
ascii-art for the exclamation mark sign.
"""

from __future__ import print_function
import os
import sys
import numpy as np


try:
    xrange
except NameError:
    xrange = range


def CreateOutpathFunc(data):
    """
    Creat RESULTS folder for results output. 
    """

    # output file
    folder_name = data.conf['out_folder']
    # output path
    folder_path = os.path.join(data.paths['data'], folder_name)

    if not os.path.isdir(folder_path):
        os.makedirs(folder_path)
        print('Created folder for results: \n', folder_path, '\n')
    else:
        print('folder', folder_path, 'already exists. Using it for results output.')


def ReadDataVectorFunc(data, nzbins, nzcorrs):
    """
    Read data vector
    """
    ntheta = data.const['ntheta']
    data_path = data.paths['data']

    # plus one for theta-column
    data_xip = np.zeros((ntheta, nzcorrs + 1))
    data_xim = np.zeros((ntheta, nzcorrs + 1))
    idx_corr = 0
    for zbin1 in xrange(nzbins):
        for zbin2 in xrange(zbin1, nzbins):
            fname = os.path.join(data_path, 'DATA_VECTOR/KV450_xi_pm_files/KV450_xi_pm_tomo_{:}_{:}_logbin_mcor.dat'.format(zbin1 + 1, zbin2 + 1))
            theta, xip, xim = np.loadtxt(fname, unpack=True)

            # this assumes theta is the same for every tomographic bin and
            # for both xi_p and xi_m!
            if idx_corr == 0:
                data_xip[:, 0] = theta
                data_xim[:, 0] = theta

            data_xip[:, idx_corr + 1] = xip
            data_xim[:, idx_corr + 1] = xim

            idx_corr += 1

    data_xipm = np.concatenate((data_xip, data_xim))

    print('Loaded data vectors from: \n', os.path.join(data_path, 'DATA_VECTOR/KV450_xi_pm_files/'), '\n')

    return data_xipm


def ReadCutValueFunc(data, nzbins, nzcorrs, theta_bins):
    """
    Read cut values and convert into mask
    """
    data_path = data.paths['data']
    cutvalues_file = data.conf['cutvalues_file']

    ntheta = data.const['ntheta']

    cutvalues_file_path = os.path.join(data_path, 'SUPPLEMENTARY_FILES/CUT_VALUES/' + cutvalues_file)
    if os.path.exists(cutvalues_file_path):
        cut_values = np.loadtxt(cutvalues_file_path)
    else:
        raise Exception('File not found:\n {:} \n \
            Check that requested file exists in the following folder: \
            \n {:}'.format(cutvalues_file_path))

    # create the mask
    mask = np.zeros(2 * nzcorrs * ntheta)
    iz = 0
    for izl in xrange(nzbins):
        for izh in xrange(izl, nzbins):
            # this counts the bin combinations
            # iz=1 =>(1,1), iz=2 =>(1,2) etc
            iz = iz + 1
            for i in xrange(ntheta):
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
    nz_method = data.conf['nz_method']


    # Create labels for loading of dn/dz-files:
    zbin_labels = []
    for i in xrange(nzbins):
        zbin_labels += ['{:.1f}t{:.1f}'.format(z_bins_min[i], z_bins_max[i])]

    # Read fiducial dn_dz from window files:
    z_samples = []
    hist_samples = []
    for zbin in xrange(nzbins):
        window_file_path = os.path.join(
            data.paths['data'], 'REDSHIFT_DISTRIBUTIONS/Nz_{0:}/Nz_{0:}_Mean/Nz_{0:}_z{1:}.asc'.format(nz_method, zbin_labels[zbin]))
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
            data.paths['data'], 'REDSHIFT_DISTRIBUTIONS/Nz_{0:}/Nz_{0:}_Mean/'.format(nz_method)), '\n')

    return np.asarray(z_samples), np.asarray(hist_samples), len(zptemp)



def WriteVectorFunc(data, nzbins, theta_bins, mask_indices, mask_suffix, vec, fname_prefix='your_filename_prefix_here'):
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
    for idx_z1 in xrange(nzbins):
        for idx_z2 in xrange(idx_z1, nzbins):
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
    fname = os.path.join(data_path, '{:}/{:}_cut_to_{:}.dat'.format(out_folder,fname_prefix, mask_suffix))
    np.savetxt(fname, savedata, header=header, delimiter='\t', fmt=['%4i', '%.5e', '%12.5e', '%i', '%i', '%i'])
    print('Saved vector in list format cut down to scales as specified in {:}: \n'.format(cutvalues_file), fname, '\n')