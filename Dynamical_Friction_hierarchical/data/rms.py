#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
# import matplotlib.animation as animation


# from matplotlib.font_manager import FontProperties
# fp = FontProperties(fname='/Users/isoya.kazuhide/Library/Fonts/ipag.ttf');


######################################################################
directory = "Ntr1E2_t1E8_dtlog_Mtot6E-7_Mmax3E-9_ecc5E-2_frag_acc/"
# directory = "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax3E-9_ecc5E-2_frag_acc/"

N_p = 1
N_tr = 100

if(N_p == 1):
    # LINE = 42  # 10000yr
    # LINE = 314  # 1E8yr
    LINE = 279  # 5.23
    # SUBDIR_NUM = 40
    SUBDIR_NUM = 1
elif(N_p == 3):
    LINE = 34  # 1000yr
    SUBDIR_NUM = 13
    # SUBDIR_NUM = 12  # Mmax5E-15_ecc1E-2_frag_acc
    # SUBDIR_NUM = 13  # Mmax5E-15_ecc3E-2_frag_acc
    # SUBDIR_NUM = 19  # Mmax5E-15_ecc5E-2_frag_acc


time = np.empty([SUBDIR_NUM, LINE], dtype=float)  # (ファイル番号,行数)
ecc_1 = np.empty([SUBDIR_NUM, LINE], dtype=float)
ecc_p_ms = np.empty([SUBDIR_NUM, LINE], dtype=float)
ecc_tr_ms = np.empty([SUBDIR_NUM, LINE], dtype=float)
inc_1 = np.empty([SUBDIR_NUM, LINE], dtype=float)
inc_p_ms = np.empty([SUBDIR_NUM, LINE], dtype=float)
inc_tr_ms = np.empty([SUBDIR_NUM, LINE], dtype=float)
if(N_p == 3):
    ecc_2 = np.empty([SUBDIR_NUM, LINE], dtype=float)
    ecc_3 = np.empty([SUBDIR_NUM, LINE], dtype=float)
    inc_2 = np.empty([SUBDIR_NUM, LINE], dtype=float)
    inc_3 = np.empty([SUBDIR_NUM, LINE], dtype=float)
#####


for subnum in range(1, SUBDIR_NUM+1):
    subdirectory = "rand%02d/" % subnum

    arr = np.genfromtxt(directory + subdirectory + "RMS_Ntr%3d.dat" % N_tr, dtype=np.float, delimiter="\t")
    print(arr.shape)

    if(N_p == 1):
        time[subnum-1, :] = arr[:, 0]
        ecc_1[subnum-1, :] = arr[:, 1] ** 2
        ecc_p_ms[subnum-1, :] = arr[:, 2] ** 2
        ecc_tr_ms[subnum-1, :] = arr[:, 3] ** 2
        inc_1[subnum-1, :] = arr[:, 4] ** 2
        inc_p_ms[subnum-1, :] = arr[:, 5] ** 2
        inc_tr_ms[subnum-1, :] = arr[:, 6] ** 2
    elif(N_p == 3):
        time[subnum-1, :] = arr[:, 0]
        ecc_1[subnum-1, :] = arr[:, 1] ** 2
        ecc_2[subnum-1, :] = arr[:, 2] ** 2
        ecc_3[subnum-1, :] = arr[:, 3] ** 2
        ecc_p_ms[subnum-1, :] = arr[:, 4] ** 2
        ecc_tr_ms[subnum-1, :] = arr[:, 5] ** 2
        inc_1[subnum-1, :] = arr[:, 6] ** 2
        inc_2[subnum-1, :] = arr[:, 7] ** 2
        inc_3[subnum-1, :] = arr[:, 8] ** 2
        inc_p_ms[subnum-1, :] = arr[:, 9] ** 2
        inc_tr_ms[subnum-1, :] = arr[:, 10] ** 2

    print(subnum, time[subnum-1, 1], np.sqrt(ecc_p_ms[subnum-1, 1]), np.sqrt(inc_p_ms[subnum-1, 1]))


######################################################################
if(N_p == 1):
    rand_rms = np.zeros([13, LINE], dtype=float)
    rand_rms[1, :] = np.sqrt(ecc_1.mean(axis=0))
    rand_rms[2, :] = 0.5 * ecc_1.std(axis=0, ddof=1) / np.sqrt(ecc_1.mean(axis=0))
    rand_rms[3, :] = np.sqrt(ecc_p_ms.mean(axis=0))
    rand_rms[4, :] = 0.5 * ecc_p_ms.std(axis=0, ddof=1) / np.sqrt(ecc_p_ms.mean(axis=0))
    rand_rms[5, :] = np.sqrt(ecc_tr_ms.mean(axis=0))
    rand_rms[6, :] = 0.5 * ecc_tr_ms.std(axis=0, ddof=1) / np.sqrt(ecc_tr_ms.mean(axis=0))
    rand_rms[7, :] = np.sqrt(inc_1.mean(axis=0))
    rand_rms[8, :] = 0.5 * inc_1.std(axis=0, ddof=1) / np.sqrt(inc_1.mean(axis=0))
    rand_rms[9, :] = np.sqrt(inc_p_ms.mean(axis=0))
    rand_rms[10, :] = 0.5 * inc_p_ms.std(axis=0, ddof=1) / np.sqrt(inc_p_ms.mean(axis=0))
    rand_rms[11, :] = np.sqrt(inc_tr_ms.mean(axis=0))
    rand_rms[12, :] = 0.5 * inc_tr_ms.std(axis=0, ddof=1) / np.sqrt(inc_tr_ms.mean(axis=0))
elif(N_p == 3):
    rand_rms = np.zeros([21, LINE], dtype=float)
    rand_rms[1, :] = np.sqrt(ecc_1.mean(axis=0))
    rand_rms[2, :] = 0.5 * ecc_1.std(axis=0, ddof=1) / np.sqrt(ecc_1.mean(axis=0))
    rand_rms[3, :] = np.sqrt(ecc_2.mean(axis=0))
    rand_rms[4, :] = 0.5 * ecc_2.std(axis=0, ddof=1) / np.sqrt(ecc_2.mean(axis=0))
    rand_rms[5, :] = np.sqrt(ecc_3.mean(axis=0))
    rand_rms[6, :] = 0.5 * ecc_3.std(axis=0, ddof=1) / np.sqrt(ecc_3.mean(axis=0))
    rand_rms[7, :] = np.sqrt(ecc_p_ms.mean(axis=0))
    rand_rms[8, :] = 0.5 * ecc_p_ms.std(axis=0, ddof=1) / np.sqrt(ecc_p_ms.mean(axis=0))
    rand_rms[9, :] = np.sqrt(ecc_tr_ms.mean(axis=0))
    rand_rms[10, :] = 0.5 * ecc_tr_ms.std(axis=0, ddof=1) / np.sqrt(ecc_tr_ms.mean(axis=0))
    rand_rms[11, :] = np.sqrt(inc_1.mean(axis=0))
    rand_rms[12, :] = 0.5 * inc_1.std(axis=0, ddof=1) / np.sqrt(inc_1.mean(axis=0))
    rand_rms[13, :] = np.sqrt(inc_2.mean(axis=0))
    rand_rms[14, :] = 0.5 * inc_2.std(axis=0, ddof=1) / np.sqrt(inc_2.mean(axis=0))
    rand_rms[15, :] = np.sqrt(inc_3.mean(axis=0))
    rand_rms[16, :] = 0.5 * inc_3.std(axis=0, ddof=1) / np.sqrt(inc_3.mean(axis=0))
    rand_rms[17, :] = np.sqrt(inc_p_ms.mean(axis=0))
    rand_rms[18, :] = 0.5 * inc_p_ms.std(axis=0, ddof=1) / np.sqrt(inc_p_ms.mean(axis=0))
    rand_rms[19, :] = np.sqrt(inc_tr_ms.mean(axis=0))
    rand_rms[20, :] = 0.5 * inc_tr_ms.std(axis=0, ddof=1) / np.sqrt(inc_tr_ms.mean(axis=0))

rand_rms[0, :] = time[0, :]

# print(rand_rms[4, :])
######################################################################
plt.figure(figsize=(10, 8), dpi=100)

if(N_p == 1):
    plt.title(r"$N_{\rm tr}=1000,M_{\rm tot}=10 {\rm M_{\oplus}}$", fontsize=18)
    plt.plot(time[0, :], rand_rms[3, :], color="r", label=r"$e_{\rm p}$")
    plt.errorbar(time[0, :], rand_rms[3, :], yerr=rand_rms[4, :], fmt=".")
    plt.plot(time[0, :], rand_rms[5, :], color="b", label=r"$e_{\rm tr,rms}$")
    plt.errorbar(time[0, :], rand_rms[5, :], yerr=rand_rms[6, :], fmt=".")
elif(N_p == 3):
    plt.title(r"$N_{\rm tr}=3000,M_{\rm tot}=30 {\rm M_{\oplus}}$", fontsize=18)
    plt.plot(time[0, :], rand_rms[1, :], color="c", linestyle="dashed", linewidth=0.5, label=r"$e_{\rm p,1}$")
    plt.errorbar(time[0, :], rand_rms[1, :], yerr=rand_rms[2, :], fmt=".")
    plt.plot(time[0, :], rand_rms[3, :], color="m", linestyle="dashed", linewidth=0.5, label=r"$e_{\rm p,2}$")
    plt.errorbar(time[0, :], rand_rms[3, :], yerr=rand_rms[4, :], fmt=".")
    plt.plot(time[0, :], rand_rms[5, :], color="y", linestyle="dashed", linewidth=0.5, label=r"$e_{\rm p,3}$")
    plt.errorbar(time[0, :], rand_rms[5, :], yerr=rand_rms[6, :], fmt=".")
    plt.plot(time[0, :], rand_rms[7, :], color="r", label=r"$e_{\rm p,rms}$")
    plt.errorbar(time[0, :], rand_rms[7, :], yerr=rand_rms[8, :], fmt=".")
    plt.plot(time[0, :], rand_rms[9, :], color="b", label=r"$e_{\rm tr,rms}$")
    plt.errorbar(time[0, :], rand_rms[9, :], yerr=rand_rms[10, :], fmt=".")

plt.xscale("log")
plt.yscale("log")
plt.xlim([0.1, 1E8+1])
# plt.ylim([0.03, 0.12])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.xlabel('time[yr]', fontsize=25)
plt.ylabel('ecc', fontsize=25)
plt.grid(True)
plt.legend(loc="upper left", fontsize=15)

plt.tight_layout()

plt.show()

plt.close()

#####################################################################
plt.figure(figsize=(10, 8), dpi=100)

if(N_p == 1):
    plt.title(r"$N_{\rm tr}=1000,M_{\rm tot}=10 {\rm M_{\oplus}}$", fontsize=18)
    plt.plot(time[0, :], rand_rms[9, :], color="r", label=r"$i_{\rm p}$")
    plt.errorbar(time[0, :], rand_rms[9, :], yerr=rand_rms[10, :], fmt=".")
    plt.plot(time[0, :], rand_rms[11, :], color="b", label=r"$i_{\rm tr,rms}$")
    plt.errorbar(time[0, :], rand_rms[11, :], yerr=rand_rms[12, :], fmt=".")
elif(N_p == 3):
    plt.title(r"$N_{\rm tr}=3000,M_{\rm tot}=30 {\rm M_{\oplus}}$", fontsize=18)
    plt.plot(time[0, :], rand_rms[11, :], color="c", linestyle="dashed", linewidth=0.5, label=r"$i_{\rm p,1}$")
    plt.errorbar(time[0, :], rand_rms[11, :], yerr=rand_rms[12, :], fmt=".")
    plt.plot(time[0, :], rand_rms[13, :], color="m", linestyle="dashed", linewidth=0.5, label=r"$i_{\rm p,2}$")
    plt.errorbar(time[0, :], rand_rms[13, :], yerr=rand_rms[14, :], fmt=".")
    plt.plot(time[0, :], rand_rms[15, :], color="y", linestyle="dashed", linewidth=0.5, label=r"$i_{\rm p,3}$")
    plt.errorbar(time[0, :], rand_rms[15, :], yerr=rand_rms[16, :], fmt=".")
    plt.plot(time[0, :], rand_rms[17, :], color="r", label=r"$i_{\rm p,rms}$")
    plt.errorbar(time[0, :], rand_rms[17, :], yerr=rand_rms[18, :], fmt=".")
    plt.plot(time[0, :], rand_rms[19, :], color="b", label=r"$i_{\rm tr,rms}$")
    plt.errorbar(time[0, :], rand_rms[19, :], yerr=rand_rms[20, :], fmt=".")


plt.xscale("log")
plt.yscale("log")
plt.xlim([0.1, 1E8+1])
# plt.ylim([0.03, 0.06])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.xlabel('time[yr]', fontsize=25)
plt.ylabel('inc [rad]', fontsize=25)
plt.grid(True)
plt.legend(loc="upper left", fontsize=15)


plt.tight_layout()

plt.show()

plt.close()

#####################################################################


if (N_p == 1):
    np.savetxt(directory + "RMS_randall.dat", rand_rms[:, [i for i in range(LINE) if i < 58 or (i >= 58 and i % 16 == 9)]].T, fmt="%.15e", delimiter="\t", newline="\n",
               header="time\tecc_1\tecc_1_error\tecc_p_rms\tecc_p_rms_error\tecc_tr_rms\tecc_tr_rms_error\tinc_1\tinc_1_error\tinc_p_rms\tinc_p_rms_error\tinc_tr_rms\tinc_tr_rms_error")
elif (N_p == 3):
    np.savetxt(directory + "RMS_randall.dat", rand_rms.T, fmt="%.15e", delimiter="\t", newline="\n",
               header="time\tecc_1\tecc_1_error\tecc_2\tecc_2_error\tecc_3\tecc_3_error\tecc_p_rms\tecc_p_rms_error\tecc_tr_rms\tecc_tr_rms_error\tinc_1\tinc_1_error\tinc_2\tinc_2_error\tinc_3\tinc_3_error\tinc_p_rms\tinc_p_rms_error\tinc_tr_rms\tinc_tr_rms_error")
