# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-05-05 21:17:39
# @Last Modified by:   lshuns
# @Last Modified time: 2020-05-25 15:50:28

# covariance matrix plot with heatmap

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import LogNorm


# +++++++++++++++++++++++++++++++++++++++ general setting
# covariance type
#   mar11 for old one
#   apr8 for updated one
# cov_tag = 'mar11'
cov_tag = 'apr8'

# parent directory to the input covariance
ParDir = "/disks/shear15/ssli/CosmicShear/covariance/"

# path to save the plot
outfile_pdf = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/cov_comb_br.pdf"
outfile_png = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/Cov/cov_comb_br.png"

# covariance matrix
# blue part
inpath = ParDir + cov_tag + '/thps_cov_{:}_bb_inc_m_usable.dat'.format(cov_tag)
cov_bb = np.loadtxt(inpath)
# print("sum of bb", np.trace(cov_bb))
# print("abs sum of bb", np.trace(np.absolute(cov_bb)))

# red part
inpath = ParDir + cov_tag + '/thps_cov_{:}_rr_inc_m_usable.dat'.format(cov_tag)
cov_rr = np.loadtxt(inpath)
# print("sum of rr", np.trace(cov_rr))
# print("abs sum of rr", np.trace(np.absolute(cov_rr)))

# cross part
inpath = ParDir + cov_tag + '/thps_cov_{:}_br_inc_m_usable.dat'.format(cov_tag)
cov_br = np.loadtxt(inpath)
cov_rb = cov_br.transpose()

# combine to joint covariance
cov_comb = np.asarray(np.bmat('cov_bb, cov_br; cov_rb, cov_rr'))

# normalize using diagonal
cov_comb_dia = np.diag(cov_comb)
cov_comb_norm = np.sqrt(np.outer(cov_comb_dia, cov_comb_dia))

#
cov_comb_plot = cov_comb/cov_comb_norm

# # +++++++++++++++++++++++++++++++++++++++++ plot
# general settings for plot
# mpl.use('Agg')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
plt.rc('font', size=12)

### normalized 
fig, ax = plt.subplots(1, 1)
# cmap = sns.diverging_palette(250, 0, sep=10, as_cmap=True)
sns.heatmap(
    cov_comb_plot, ax=ax, cmap='Reds',
    xticklabels=False, yticklabels=False,
    square=True,
    cbar_kws={'label': r'$C_{ij}/\sqrt{C_{ii}C_{jj}}$', 'ticks': [0.00, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90]}
)

midp = len(cov_comb_plot[0,:])/2
ax.axvline(x=midp, color='gray', ls='--', lw=1)
ax.axhline(y=midp, color='gray', ls='--', lw=1)

fig.text(0.27, 0.06, r'$\xi_{\rm blue}$')
fig.text(0.57, 0.06, r'$\xi_{\rm red}$')
fig.text(0.13, 0.28, r'$\xi_{\rm red}$', rotation=90)
fig.text(0.13, 0.72, r'$\xi_{\rm blue}$', rotation=90)


fig.savefig(outfile_pdf, dpi=300)
print("Plot saved in", outfile_pdf)
fig.savefig(outfile_png, dpi=300)
print("Plot saved in", outfile_png)
plt.close()
