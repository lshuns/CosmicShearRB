#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 16:41:11 2020
@author: ssli

Post process of covariance calculation
    re-arrange the results to desired forms
"""

import numpy as np
import pandas as pd

# cut matrix into three different parts:
#       bb, rr, rb (=br)

# +++++++++++++++++++++++++++++++++++++++ general setting
# covariance type
#   mar11 for old one
#   apr8 for updated one
# cov_tag = 'mar11'
cov_tag = 'apr8'

# parent directory for input & output covariance
ParDir = "/disks/shear15/ssli/CosmicShear/covariance/"

# path to the theory vector
P_xi_theo_r = "/disks/shear15/ssli/CosmicShear/theory_vector/xi_theory_full_less3_KV450_best.dat"
P_xi_theo_b = "/disks/shear15/ssli/CosmicShear/theory_vector/xi_theory_full_greater3_KV450_best.dat"

# +++++++++++++++++++++++++++++++++++++++ to list form
inpath = ParDir + cov_tag + "/original/thps_cov_{:}_list.dat".format(cov_tag)

tmp_raw = np.loadtxt(inpath)

# discard undesired columns
tmp_raw = np.delete(tmp_raw, [2, 3, 8, 9], axis=1)

# build dataframe
df = pd.DataFrame(tmp_raw, \
    columns=['s1_bin1','s1_bin2','s1_xip0_xim1', 's1_theta', \
            's2_bin1','s2_bin2','s2_xip0_xim1', 's2_theta', \
            'Gaussian', 'non_Gaussian'])

df = df.astype(dtype= {"s1_bin1":"int32",\
                        "s1_bin2":"int32",\
                        "s1_xip0_xim1":"int32",\
                        "s1_theta":"float64",\
                        "s2_bin1":"int32",\
                        "s2_bin2":"int32",\
                        "s2_xip0_xim1":"int32",\
                        "s2_theta":"float64",\
                        "Gaussian":"float64",\
                        "non_Gaussian":"float64",\
                        })

# get blue-blue part
mask_bb = (df.s1_bin1<=5) & (df.s1_bin2<=5) & (df.s2_bin1<=5) & (df.s2_bin2<=5)
df_bb = df[mask_bb]

# get red-red part
mask_rr = (df.s1_bin1>5) & (df.s1_bin2>5) & (df.s2_bin1>5) & (df.s2_bin2>5)
df_rr = df[mask_rr]

# get blue-red part
mask_br = (df.s1_bin1<=5) & (df.s1_bin2<=5) & (df.s2_bin1>5) & (df.s2_bin2>5)
df_br = df[mask_br]


# output
outpath1 = ParDir + cov_tag + '/thps_cov_{:}_bb_list.dat'.format(cov_tag)
df_bb.to_csv(outpath1, sep=' ', index=False, header=False)
#
outpath2 = ParDir + cov_tag + '/thps_cov_{:}_rr_list.dat'.format(cov_tag)
df_rr.to_csv(outpath2, sep=' ', index=False, header=False)
#
outpath3 = ParDir + cov_tag + '/thps_cov_{:}_br_list.dat'.format(cov_tag)
df_br.to_csv(outpath3, sep=' ', index=False, header=False)
print("list form of covariance saved to: \n", outpath1, '\n', outpath2, '\n', outpath3, '\n')

# +++++++++++++++++++++++++++++++++++++++ to matrix form
print('Now we construct the covariance matrix in a format usable for cosmological calculation.')

def List2UsableFunc(tmp_raw, xi_theo1=None, xi_theo2=None, CROSS=False):

        ntheta = 9
        nzcorrs = int(5 * (5 + 1) / 2)

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
        if CROSS:
            # for cross covariance
            for index1 in range(dim):
                for index2 in range(dim):
                    matrix[index1, index2] = values[index_lin]
                    index_lin += 1
        else:
            # for auto covariance
            for index1 in range(dim):
                for index2 in range(index1, dim):
                    matrix[index1, index2] = values[index_lin]
                    matrix[index2, index1] = matrix[index1, index2]
                    index_lin += 1


        # apply propagation of m-correction uncertainty following
        # equation 12 from Hildebrandt et al. 2017 (arXiv:1606.05338):
        err_multiplicative_bias = 0.02
        matrix_m_corr = np.matrix(xi_theo1).T * np.matrix(xi_theo2) * 4. * err_multiplicative_bias**2
        matrix = matrix + np.asarray(matrix_m_corr)

        return matrix

# load xi_theo
xi_theo_r = np.loadtxt(P_xi_theo_r)[:,2]
xi_theo_b = np.loadtxt(P_xi_theo_b)[:,2]

# list form covariance
df_bb = df_bb.to_numpy()
df_rr = df_rr.to_numpy()
df_br = df_br.to_numpy()

# transfer to matrix form
df_bb = List2UsableFunc(df_bb, xi_theo_b, xi_theo_b, False)
outpath1 = ParDir + cov_tag + '/thps_cov_{:}_bb_inc_m_usable.dat'.format(cov_tag)
np.savetxt(outpath1, df_bb)

df_rr = List2UsableFunc(df_rr, xi_theo_r, xi_theo_r, False)
outpath2 = ParDir + cov_tag + '/thps_cov_{:}_rr_inc_m_usable.dat'.format(cov_tag)
np.savetxt(outpath2, df_rr)

df_br = List2UsableFunc(df_br, xi_theo_b, xi_theo_r, True)
outpath3 = ParDir + cov_tag + '/thps_cov_{:}_br_inc_m_usable.dat'.format(cov_tag)
np.savetxt(outpath3, df_br)

print('Saved covariance matrix (incl. shear calibration uncertainty) in format usable with this likelihood to: \n', outpath1, '\n', outpath2, '\n', outpath3, '\n')
