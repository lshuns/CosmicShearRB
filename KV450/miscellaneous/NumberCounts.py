#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:27:44 2019

@author: ssli

Count object numbers in different scenarios

Bad structure due to the bad log.txt :(
"""

import numpy as np 

# pre (scenario without magnitude cut)
# mine (scenario without magnitude cut)

# patch total_number selected_number_pre selected_number_mine
G9 = [1771080, 1729203, 1192193]
GS = [2123705, 2023890, 1394009]
G12 = [3575185, 3224753, 2224676]
G15 = [3437991, 3152237, 2144398]
G23 = [3313121, 3111724, 2226579]
# path_bin number
G12_1 = 254499
G12_3 = 427916
G12_5 = 900781
G12_7 = 631899
G12_9 = 669513
G12_b_pre = [G12_1, G12_3, G12_5, G12_7, G12_9]
G15_1 = 246973
G15_3 = 450479
G15_5 = 911409
G15_7 = 615506
G15_9 = 623776
G15_b_pre = [G15_1, G15_3, G15_5, G15_7, G15_9]
GS_1 = 144917
GS_3 = 254886
GS_5 = 507620
GS_7 = 366507
GS_9 = 493295
GS_b_pre = [GS_1, GS_3, GS_5, GS_7, GS_9]
G9_1 = 132439
G9_3 = 249964
G9_5 = 490090
G9_7 = 302001
G9_9 = 346097
G9_b_pre = [G9_1, G9_3, G9_5, G9_7, G9_9]
G23_1 = 248676
G23_3 = 415585
G23_5 = 828908
G23_7 = 571211
G23_9 = 648995
G23_b_pre = [G23_1, G23_3, G23_5, G23_7, G23_9]

G15_b_mine = [228172, 342013, 593893, 390994, 391927]
GS_b_mine = [135848, 200389, 339151, 233655, 317025]
G9_b_mine = [123364, 190640, 325517, 196313, 225800]
G23_b_mine = [233362, 330171, 567525, 387907, 441086]
G12_b_mine = [237498, 331140, 598375, 408560, 426184]

# band -99 99
out_r = [0, 703]
out_u = [0, 2674739]
out_g = [0, 127355]
out_i = [0, 610885]
out_Z = [100, 107673]
out_Y = [44, 330918]
out_J = [0, 112839]
out_H = [36, 741275]
out_Ks = [12, 639986]
out_tot = [192, 5346373]


# Total Number difference between different data sets
N_tot = G9[0] + GS[0] + G12[0] + G15[0] + G23[0]
N_pre = G9[1] + GS[1] + G12[1] + G15[1] + G23[1]
N_mine = G9[2] + GS[2] + G12[2] + G15[2] + G23[2]
print("N_tot", N_tot)
print("N_pre", N_pre)
print("N_mine", N_mine)
print("N_pre-N_mine", N_pre - N_mine)

# Number in bins
Nb_pre = 0
Nb_mine = 0
for i in range(5):
    N_pre = G12_b_pre[i] + G23_b_pre[i] + G9_b_pre[i] + GS_b_pre[i] + G15_b_pre[i]
    N_mine = G12_b_mine[i] + G23_b_mine[i] + G9_b_mine[i] + GS_b_mine[i] + G15_b_mine[i]

    print("N_pre in bin", i, ":", N_pre)
    print("N_mine in bin", i, ":", N_mine)

    Nb_pre += N_pre
    Nb_mine += N_mine

print("N_pre in all bins:", Nb_pre)
print("N_mine in all bins:", Nb_mine)