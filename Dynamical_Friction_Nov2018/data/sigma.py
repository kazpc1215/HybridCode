#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
# import matplotlib.animation as animation


# from matplotlib.font_manager import FontProperties
# fp = FontProperties(fname='/Users/isoya.kazuhide/Library/Fonts/ipag.ttf');


######################################################################
# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc1E-2_frag_acc/"
# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc3E-2_frag_acc/"
# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc5E-2_frag_acc/"

# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc1E-2_frag_acc/"
# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc3E-2_frag_acc/"
directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc5E-2_frag_acc/"

N_p = 3
N_tr = 3000

if(N_p == 1):
    LINE = 42  # 10000yr
    SUBDIR_NUM = 40
elif(N_p == 3):
    LINE = 34  # 1000yr
    SUBDIR_NUM = 13
    # SUBDIR_NUM = 12  # Mmax5E-15_ecc1E-2_frag_acc
    # SUBDIR_NUM = 13  # Mmax5E-15_ecc3E-2_frag_acc
    # SUBDIR_NUM = 19  # Mmax5E-15_ecc5E-2_frag_acc


time = np.empty([SUBDIR_NUM, LINE], dtype=float)  # (ファイル番号,行数)
mass_all = np.empty([SUBDIR_NUM, LINE], dtype=float)
mass_inner = np.empty([SUBDIR_NUM, LINE], dtype=float)
mass_center = np.empty([SUBDIR_NUM, LINE], dtype=float)
mass_outer = np.empty([SUBDIR_NUM, LINE], dtype=float)
sigma_inner = np.empty([SUBDIR_NUM, LINE], dtype=float)
sigma_center = np.empty([SUBDIR_NUM, LINE], dtype=float)
sigma_outer = np.empty([SUBDIR_NUM, LINE], dtype=float)
n_inner = np.empty([SUBDIR_NUM, LINE], dtype=float)
n_center = np.empty([SUBDIR_NUM, LINE], dtype=float)
n_outer = np.empty([SUBDIR_NUM, LINE], dtype=float)
#####


for subnum in range(1, SUBDIR_NUM+1):
    subdirectory = "rand%02d/" % subnum

    arr = np.genfromtxt(directory + subdirectory + "Sigma_dep.dat", dtype=np.float, delimiter="\t")
    print(arr.shape)

    time[subnum-1, :] = arr[:, 0]
    mass_all[subnum-1, :] = arr[:, 1] / arr[0, 1]  # 初期値で規格化
    mass_inner[subnum-1, :] = arr[:, 2] / arr[0, 2]
    mass_center[subnum-1, :] = arr[:, 3] / arr[0, 3]
    mass_outer[subnum-1, :] = arr[:, 4] / arr[0, 4]
    sigma_inner[subnum-1, :] = arr[:, 5] / arr[0, 5]
    sigma_center[subnum-1, :] = arr[:, 6] / arr[0, 6]
    sigma_outer[subnum-1, :] = arr[:, 7] / arr[0, 7]
    n_inner[subnum-1, :] = arr[:, 8] / arr[0, 8]
    n_center[subnum-1, :] = arr[:, 9] / arr[0, 9]
    n_outer[subnum-1, :] = arr[:, 10] / arr[0, 10]

    print(subnum, time[subnum-1, :], mass_center[subnum-1, :], sigma_center[subnum-1, :])


######################################################################

rand_error = np.zeros([9, LINE], dtype=float)
rand_error[0, :] = time[0, :]
rand_error[1, :] = mass_all.mean(axis=0)
rand_error[2, :] = mass_all.std(axis=0, ddof=1)
rand_error[3, :] = mass_center.mean(axis=0)
rand_error[4, :] = mass_center.std(axis=0, ddof=1)
rand_error[5, :] = sigma_center.mean(axis=0)
rand_error[6, :] = sigma_center.std(axis=0, ddof=1)
rand_error[7, :] = n_center.mean(axis=0)
rand_error[8, :] = n_center.std(axis=0, ddof=1)

######################################################################

plt.figure(figsize=(10, 8), dpi=100)
plt.plot(time[0, :], rand_error[5, :], color="r", label=r"$\Sigma/\Sigma_0$")
plt.errorbar(time[0, :], rand_error[5, :], yerr=rand_error[6, :], fmt=".", color="r")
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.show()
plt.close()

#####################################################################

np.savetxt(directory + "Sigma_randall.dat", rand_error.T, fmt="%.15e", delimiter="\t", newline="\n",
           header="time\tmass_tot_all\tmass_tot_all_error\tmass_tot_center\tmass_tot_center_error\tsigma_center\tsigma_center_error\tn_center\tn_center_error\t(normalized)")
