#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
# import matplotlib.animation as animation


# from matplotlib.font_manager import FontProperties
# fp = FontProperties(fname='/Users/isoya.kazuhide/Library/Fonts/ipag.ttf');


######################################################################
directory = "S8E2_t1E3_dtlog_M1E24g_ecc1E-4/"
# directory = "S1E3_t1E3_dtlog_M1E24g_ecc1E-4/"


N = 1000


LINE = 34  # 1000yr

SUBDIR_NUM = 35


time = np.empty([SUBDIR_NUM, LINE], dtype=float)  # (ファイル番号,行数)
ecc_S_ms = np.empty([SUBDIR_NUM, LINE], dtype=float)
ecc_L_ms = np.empty([SUBDIR_NUM, LINE], dtype=float)
inc_S_ms = np.empty([SUBDIR_NUM, LINE], dtype=float)
inc_L_ms = np.empty([SUBDIR_NUM, LINE], dtype=float)
#####


for subnum in range(1, SUBDIR_NUM+1):
    subdirectory = "rand%02d/" % subnum

    arr = np.genfromtxt(directory + subdirectory + "RMS_Ntr%4d.dat" % N, dtype=np.float, delimiter="\t")
    print(arr.shape)

    time[subnum-1, :] = arr[:, 0]
    ecc_S_ms[subnum-1, :] = arr[:, 1] ** 2
    ecc_L_ms[subnum-1, :] = arr[:, 2] ** 2
    inc_S_ms[subnum-1, :] = arr[:, 3] ** 2
    inc_L_ms[subnum-1, :] = arr[:, 4] ** 2

    print(subnum, time[subnum-1, 0], np.sqrt(ecc_S_ms[subnum-1, 0]), np.sqrt(inc_S_ms[subnum-1, 0]))


######################################################################
rand_rms = np.zeros([13, LINE], dtype=float)
rand_rms[0, :] = time[0, :]
rand_rms[1, :] = np.sqrt(ecc_S_ms.mean(axis=0))
rand_rms[2, :] = 0.5 * ecc_S_ms.std(axis=0, ddof=1) / np.sqrt(ecc_S_ms.mean(axis=0))
rand_rms[3, :] = np.sqrt(ecc_L_ms.mean(axis=0))
rand_rms[4, :] = 0.5 * ecc_L_ms.std(axis=0, ddof=1) / np.sqrt(ecc_L_ms.mean(axis=0))
rand_rms[5, :] = np.sqrt(inc_S_ms.mean(axis=0))
rand_rms[6, :] = 0.5 * inc_S_ms.std(axis=0, ddof=1) / np.sqrt(inc_S_ms.mean(axis=0))
rand_rms[7, :] = np.sqrt(inc_L_ms.mean(axis=0))
rand_rms[8, :] = 0.5 * inc_L_ms.std(axis=0, ddof=1) / np.sqrt(inc_L_ms.mean(axis=0))


print(rand_rms[2, :])
######################################################################
plt.figure(figsize=(10, 8), dpi=100)

plt.plot(time[1, :], rand_rms[1, :], color="r", label=r"$<e_S^2>^{1/2}$")
plt.errorbar(time[1, :], rand_rms[1, :], yerr=rand_rms[2, :], fmt=".")
plt.plot(time[1, :], rand_rms[3, :], color="g", label=r"$<e_L^2>^{1/2}$")
plt.errorbar(time[1, :], rand_rms[3, :], yerr=rand_rms[4, :], fmt=".")
plt.plot(time[1, :], rand_rms[5, :], color="b", label=r"$<i_S^2>^{1/2}$")
plt.errorbar(time[1, :], rand_rms[5, :], yerr=rand_rms[6, :], fmt=".")
plt.plot(time[1, :], rand_rms[7, :], color="c", label=r"$<i_L^2>^{1/2}$")
plt.errorbar(time[1, :], rand_rms[7, :], yerr=rand_rms[8, :], fmt=".")

# plt.xscale("log")
# plt.yscale("log")
plt.xlim([0.1, 1E3])
plt.ylim([0.0, 0.006])
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


np.savetxt(directory + "RMS_randall.dat", rand_rms[:, [i for i in range(LINE) if i < 59 or (i >= 59 and i % 16 == 10)]].T, fmt="%.15e", delimiter="\t", newline="\n",
               header="time\tecc_S_rms\tecc_S_rms_error\tecc_L_rms\tecc_L_rms_error\tinc_S_rms\tinc_S_rms_error\tinc_L_rms\tinc_L_rms_error")
