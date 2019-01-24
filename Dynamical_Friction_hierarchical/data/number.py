#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
# import matplotlib.animation as animation


# from matplotlib.font_manager import FontProperties
# fp = FontProperties(fname='/Users/isoya.kazuhide/Library/Fonts/ipag.ttf');


######################################################################
# directory = "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-2_nofrag_acc/"
directory = "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-1_nofrag_acc/"



# LINE = 42  # 10000yr
# LINE = 314  # 1E8yr
LINE = 283  # 5.62E7yr
# SUBDIR_NUM = 40
SUBDIR_NUM = 1


time = np.empty([SUBDIR_NUM, LINE], dtype=float)  # (ファイル番号,行数)
n_inner = np.empty([SUBDIR_NUM, LINE], dtype=int)
n_center = np.empty([SUBDIR_NUM, LINE], dtype=int)
n_outer = np.empty([SUBDIR_NUM, LINE], dtype=int)

for subnum in range(1, SUBDIR_NUM+1):
    subdirectory = "rand%02d/" % subnum

    arr = np.genfromtxt(directory + subdirectory + "tracerlistnumber.dat", dtype=np.float, delimiter="\t")
    print(arr.shape)

    time[subnum-1, :] = arr[:, 0]
    n_inner[subnum-1, :] = arr[:, 1]
    n_center[subnum-1, :] = arr[:, 2]
    n_outer[subnum-1, :] = arr[:, 3]

    print(subnum, time[subnum-1, :], n_center[subnum-1, :])

rand_error = np.zeros([7, LINE], dtype=float)
rand_error[0, :] = time[0, :]
rand_error[1, :] = n_inner.mean(axis=0)
rand_error[2, :] = n_inner.std(axis=0, ddof=1)
rand_error[3, :] = n_center.mean(axis=0)
rand_error[4, :] = n_center.std(axis=0, ddof=1)
rand_error[5, :] = n_outer.mean(axis=0)
rand_error[6, :] = n_outer.std(axis=0, ddof=1)

np.savetxt(directory + subdirectory + "tracerlistnumber_2.dat", rand_error[:, [i for i in range(LINE) if i < 58 or (i >= 58 and i % 16 == 9)]].T, fmt="%.15e", delimiter="\t", newline="\n",
           header="1:t[yr]\t2:n_inner\t3:n_inner_error\t4:n_center\t5:n_center_error\t6:n_outer\t7:n_outer_error")
