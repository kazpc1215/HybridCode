#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
# import matplotlib.animation as animation


# from matplotlib.font_manager import FontProperties
# fp = FontProperties(fname='/Users/isoya.kazuhide/Library/Fonts/ipag.ttf');


######################################################################
directory = "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-2_nofrag_acc/"


N_p = 1
# N_tr = 1000

if(N_p == 1):
    # LINE = 42  # 10000yr
    LINE = 314  # 1E8yr
    # SUBDIR_NUM = 40
    SUBDIR_NUM = 1
elif(N_p == 3):
    LINE = 34  # 1000yr
    SUBDIR_NUM = 13
    # SUBDIR_NUM = 12  # Mmax5E-15_ecc1E-2_frag_acc
    # SUBDIR_NUM = 13  # Mmax5E-15_ecc3E-2_frag_acc
    # SUBDIR_NUM = 19  # Mmax5E-15_ecc5E-2_frag_acc

time = np.empty([SUBDIR_NUM, LINE], dtype=float)  # (ファイル番号,行数)
axis_1 = np.empty([SUBDIR_NUM, LINE], dtype=float)
axis_2 = np.empty([SUBDIR_NUM, LINE], dtype=float)
axis_3 = np.empty([SUBDIR_NUM, LINE], dtype=float)
r_h_1 = np.empty([SUBDIR_NUM, LINE], dtype=float)
r_h_2 = np.empty([SUBDIR_NUM, LINE], dtype=float)
r_h_3 = np.empty([SUBDIR_NUM, LINE], dtype=float)
axis_1_dif = np.empty([SUBDIR_NUM, LINE], dtype=float)
axis_2_dif = np.empty([SUBDIR_NUM, LINE], dtype=float)
axis_3_dif = np.empty([SUBDIR_NUM, LINE], dtype=float)
#####

for subnum in range(1, SUBDIR_NUM+1):
    subdirectory = "rand%02d/" % subnum

    arr = np.genfromtxt(directory + subdirectory + "Planet01.dat", dtype=np.float, delimiter="\t")

    time[subnum-1, :] = arr[:, 0]
    # ecc[subnum-1, :] = arr[:, 1]
    axis_1[subnum-1, :] = arr[:, 2]
    # u[subnum-1, :] = arr[:, 3]
    # inc[subnum-1, :] = arr[:, 4]
    # Omega[subnum-1, :] = arr[:, 5]
    # omega[subnum-1, :] = arr[:, 6]
    r_h_1[subnum-1, :] = arr[:, 7]
    # radius[subnum-1, :] = arr[:, 8]
    # mass[subnum-1, :] = arr[:, 9]
    print("Planet01", subnum, time[subnum-1, 0], axis_1[subnum-1, 0], r_h_1[subnum-1, 0])
    axis_1_dif[subnum-1, :] = (axis_1[subnum-1, :] - axis_1[subnum-1, 0])

    if (N_p == 3):
        arr = np.genfromtxt(directory + subdirectory + "Planet02.dat", dtype=np.float, delimiter="\t")

        # time[subnum-1, :] = arr[:, 0]
        # ecc[subnum-1, :] = arr[:, 1]
        axis_2[subnum-1, :] = arr[:, 2]
        # u[subnum-1, :] = arr[:, 3]
        # inc[subnum-1, :] = arr[:, 4]
        # Omega[subnum-1, :] = arr[:, 5]
        # omega[subnum-1, :] = arr[:, 6]
        r_h_2[subnum-1, :] = arr[:, 7]
        # radius[subnum-1, :] = arr[:, 8]
        # mass[subnum-1, :] = arr[:, 9]
        print("Planet02", subnum, time[subnum-1, 0], axis_2[subnum-1, 0], r_h_2[subnum-1, 0])
        axis_2_dif[subnum-1, :] = (axis_2[subnum-1, :] - axis_2[subnum-1, 0])

        arr = np.genfromtxt(directory + subdirectory + "Planet03.dat", dtype=np.float, delimiter="\t")

        # time[subnum-1, :] = arr[:, 0]
        # ecc[subnum-1, :] = arr[:, 1]
        axis_3[subnum-1, :] = arr[:, 2]
        # u[subnum-1, :] = arr[:, 3]
        # inc[subnum-1, :] = arr[:, 4]
        # Omega[subnum-1, :] = arr[:, 5]
        # omega[subnum-1, :] = arr[:, 6]
        r_h_3[subnum-1, :] = arr[:, 7]
        # radius[subnum-1, :] = arr[:, 8]
        # mass[subnum-1, :] = arr[:, 9]
        print("Planet03", subnum, time[subnum-1, 0], axis_3[subnum-1, 0], r_h_3[subnum-1, 0])
        axis_3_dif[subnum-1, :] = (axis_3[subnum-1, :] - axis_3[subnum-1, 0])

if (N_p == 1):
    rand_error = np.zeros([7, LINE], dtype=float)
    rand_error[0, :] = time[0, :]
    rand_error[1, :] = axis_1.mean(axis=0)
    rand_error[2, :] = axis_1.std(axis=0, ddof=1)
    rand_error[3, :] = axis_1_dif.mean(axis=0)
    rand_error[4, :] = axis_1_dif.std(axis=0, ddof=1)
    rand_error[5, :] = r_h_1.mean(axis=0)
    rand_error[6, :] = r_h_1.std(axis=0, ddof=1)
elif (N_p == 3):
    rand_error = np.zeros([19, LINE], dtype=float)
    rand_error[0, :] = time[0, :]
    rand_error[1, :] = axis_1.mean(axis=0)
    rand_error[2, :] = axis_1.std(axis=0, ddof=1)
    rand_error[3, :] = axis_2.mean(axis=0)
    rand_error[4, :] = axis_2.std(axis=0, ddof=1)
    rand_error[5, :] = axis_3.mean(axis=0)
    rand_error[6, :] = axis_3.std(axis=0, ddof=1)
    rand_error[7, :] = axis_1_dif.mean(axis=0)
    rand_error[8, :] = axis_1_dif.std(axis=0, ddof=1)
    rand_error[9, :] = axis_2_dif.mean(axis=0)
    rand_error[10, :] = axis_2_dif.std(axis=0, ddof=1)
    rand_error[11, :] = axis_3_dif.mean(axis=0)
    rand_error[12, :] = axis_3_dif.std(axis=0, ddof=1)
    rand_error[13, :] = r_h_1.mean(axis=0)
    rand_error[14, :] = r_h_1.std(axis=0, ddof=1)
    rand_error[15, :] = r_h_2.mean(axis=0)
    rand_error[16, :] = r_h_2.std(axis=0, ddof=1)
    rand_error[17, :] = r_h_3.mean(axis=0)
    rand_error[18, :] = r_h_3.std(axis=0, ddof=1)

    print(rand_error)


plt.figure(figsize=(10, 8), dpi=100)
plt.plot(time[0, :], rand_error[1, :], color="r", label=r"$<a_1>$")
plt.errorbar(time[0, :], rand_error[1, :], yerr=rand_error[2, :], fmt=".", color="r")
if (N_p == 3):
    plt.plot(time[0, :], rand_error[3, :], color="g", label=r"$<a_2>$")
    plt.errorbar(time[0, :], rand_error[3, :], yerr=rand_error[4, :], fmt=".", color="g")
    plt.plot(time[0, :], rand_error[5, :], color="b", label=r"$<a_3>$")
    plt.errorbar(time[0, :], rand_error[5, :], yerr=rand_error[6, :], fmt=".", color="b")
plt.xscale("log")
plt.legend()
plt.show()
plt.close()



plt.figure(figsize=(10, 8), dpi=100)
plt.plot(time[0, :], axis_1_dif.mean(axis=0), color="r", label=r"$<a_1>$")
plt.errorbar(time[0, :], axis_1_dif.mean(axis=0), yerr=axis_1.std(axis=0, ddof=1), fmt=".", color="r")
if (N_p == 3):
    plt.plot(time[0, :], axis_2_dif.mean(axis=0), color="g", label=r"$<a_2>$")
    plt.errorbar(time[0, :], axis_2_dif.mean(axis=0), yerr=axis_2.std(axis=0, ddof=1), fmt=".", color="g")
    plt.plot(time[0, :], axis_3_dif.mean(axis=0), color="b", label=r"$<a_3>$")
    plt.errorbar(time[0, :], axis_3_dif.mean(axis=0), yerr=axis_3.std(axis=0, ddof=1), fmt=".", color="b")
plt.xscale("log")
plt.legend()
plt.show()
plt.close()


"""
plt.figure(figsize=(10, 8), dpi=100)
plt.plot(time[0, :], axis_2.mean(axis=0), color="r", label="mean")
plt.errorbar(time[0, :], axis_2.mean(axis=0), yerr=axis_2.std(axis=0, ddof=1), fmt=".", color="r")
plt.plot(time[0, :], axis_2[0, :], color="b", label="01")
plt.plot(time[0, :], axis_2[1, :], color="b", label="02")
plt.plot(time[0, :], axis_2[2, :], color="b", label="03")
plt.plot(time[0, :], axis_2[3, :], color="b", label="04")
plt.plot(time[0, :], axis_2[4, :], color="b", label="05")
plt.plot(time[0, :], axis_2[5, :], color="b", label="06")
plt.plot(time[0, :], axis_2[6, :], color="b", label="07")
plt.plot(time[0, :], axis_2[7, :], color="b", label="08")
plt.plot(time[0, :], axis_2[8, :], color="b", label="09")
plt.plot(time[0, :], axis_2[9, :], color="b", label="10")
plt.plot(time[0, :], axis_2[10, :], color="b", label="11")
plt.plot(time[0, :], axis_2[11, :], color="b", label="12")
plt.plot(time[0, :], axis_2[12, :], color="b", label="13")

plt.xscale("log")
plt.legend()
plt.show()
plt.close()
"""

if (N_p == 1):
    np.savetxt(directory + "axis_evo.dat", rand_error[:, [i for i in range(LINE) if i < 58 or (i >= 58 and i % 16 == 9)]].T, fmt="%.15e", delimiter="\t", newline="\n", header="time\taxis_1_mean\taxis_1_error\taxis_1_dif_mean\taxis_1_dif_error\tr_h_1_mean\tr_h_1_error")
elif (N_p == 3):
    np.savetxt(directory + "axis_evo.dat", rand_error.T, fmt="%.15e", delimiter="\t", newline="\n", header="time\taxis_1_mean\taxis_1_error\taxis_2_mean\taxis_2_error\taxis_3_mean\taxis_3_error\taxis_1_dif_mean\taxis_1_dif_error\taxis_2_dif_mean\taxis_2_dif_error\taxis_3_dif_mean\taxis_3_dif_error\tr_h_1_mean\tr_h_1_error\tr_h_2_mean\tr_h_2_error\tr_h_3_mean\tr_h_3_error")
