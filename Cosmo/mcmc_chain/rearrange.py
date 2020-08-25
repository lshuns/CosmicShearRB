# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2020-07-11 17:10:08
# @Last Modified by:   lshuns
# @Last Modified time: 2020-07-11 17:30:48

### re-arrange the order of chains (mainly A_IA,s)

import numpy as np

## original order
# weights, mloglkl, A_IA_m, A_IA_s, D_z1_m, D_z2_m, D_z3_m, D_z4_m, D_z5_m, D_z1_s, D_z2_s, D_z3_s, D_z4_s, D_z5_s

## new order
header = 'weights, mloglkl, A_IA_m, D_z1_m, D_z2_m, D_z3_m, D_z4_m, D_z5_m, A_IA_s, D_z1_s, D_z2_s, D_z3_s, D_z4_s, D_z5_s'

inpaths = ['./KV450_H1_Dz_IA/KV450_H1_Dz_IA.txt', './Planck_H1_Dz_IA/Planck_H1_Dz_IA.txt']

for inpath in inpaths:
    data = np.loadtxt(inpath)

    # new order
    new_idx = [0, 1, 2, 4, 5, 6, 7, 8, 3, 9, 10, 11, 12, 13]

    # new order applied
    data_new = data[:, new_idx]

    # save

    np.savetxt(inpath+'2', data_new, header=header)