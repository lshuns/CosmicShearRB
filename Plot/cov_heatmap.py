# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-05-05 21:17:39
# @Last Modified by:   lshuns
# @Last Modified time: 2020-07-16 13:10:26

# covariance matrix plot with heatmap

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import LogNorm
plt.rcParams["text.usetex"] =True
plt.rc('font', size=16)


# +++++++++++++++++++++++++++++++++++++++ general setting
# covariance type
#   mar11 for old one
#   apr8 for updated one
# cov_tag = 'mar11'
cov_tag = 'apr8'

# parent directory to the input covariance
ParDir = "/disks/shear15/ssli/CosmicShear/covariance/"

# path to save the plot
# outfile_pdf = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/cov_comb_br.pdf"
outfile_png = "/net/raam/data1/surfdrive_ssli/Projects/6CosmicShear_RB/plot/publish/cov_comb_br.png"

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
# plt.rc('font', size=14)

### normalized 
fig, ax = plt.subplots()
plt.subplots_adjust(left=0, bottom=0.06, right=0.99, top=0.98, wspace=0, hspace=0)
# cmap = sns.light_palette((360, 100, 0), input="husl", as_cmap=True)
cmap = 'Reds'
sns.heatmap(
    cov_comb_plot, ax=ax, cmap=cmap, vmin=0, vmax=1.0,
    xticklabels=False, yticklabels=False,
    square=True,
    cbar_kws={'label': r'$C_{ij}/\sqrt{C_{ii}C_{jj}}$', 'ticks': [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]}
)

midp = len(cov_comb_plot[0,:])/2
ax.axvline(x=midp, color='gray', ls='--', lw=1)
ax.axhline(y=midp, color='gray', ls='--', lw=1)

fig.text(0.27, 0.01, r'$\xi_{\rm blue}$')
fig.text(0.57, 0.01, r'$\xi_{\rm red}$')
fig.text(0.06, 0.28, r'$\xi_{\rm red}$', rotation=90)
fig.text(0.06, 0.72, r'$\xi_{\rm blue}$', rotation=90)


# fig.savefig(outfile_pdf, dpi=300)
# print("Plot saved in", outfile_pdf)
fig.savefig(outfile_png, dpi=300)
print("Plot saved in", outfile_png)
plt.close()
