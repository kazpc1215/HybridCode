#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np


# from matplotlib.font_manager import FontProperties
# fp = FontProperties(fname='/Users/isoya.kazuhide/Library/Fonts/ipag.ttf');

def Jacobi(axis):
    return(axis / x + 2.0 * np.sqrt(x / axis) * np.sqrt(1.0 - y * y))


######################################################################

directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc1E-2_nofrag_acc/"
# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc3E-2_nofrag_acc/"
# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc5E-2_nofrag_acc/"

# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc1E-2_frag_acc/"
# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc3E-2_frag_acc/"
# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-15_ecc5E-2_frag_acc/"

# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc1E-2_frag_acc/"
# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc3E-2_frag_acc/"
# directory = "Ntr3E3_t1E3_dtlog_Mtot3E-5_Mmax5E-18_ecc5E-2_frag_acc/"

N_p = 3
N_tr = 3000


if(N_p == 1):
    # LINE = 42  # 10000yr
    LINE = 299  # 1E8yr
    SUBDIR_NUM = 40
elif(N_p == 3):
    LINE = 34  # 1000yr
    SUBDIR_NUM = 13
    # SUBDIR_NUM = 12  # Mmax5E-15_ecc1E-2_frag_acc
    # SUBDIR_NUM = 13  # Mmax5E-15_ecc3E-2_frag_acc
    # SUBDIR_NUM = 19  # Mmax5E-15_ecc5E-2_frag_acc

X_MESH = 2000
Y_MESH = 5000

x = [i/X_MESH for i in range(int(0.5*X_MESH), int(1.5*X_MESH)+1)]
y = [i/Y_MESH for i in range(int(0.2*Y_MESH)+1)]
x, y = np.meshgrid(x, y)


for subnum in range(1, SUBDIR_NUM+1):
    subdirectory = "rand%02d/" % subnum

    time = np.empty([N_p+N_tr+1, LINE], dtype=float)  # (ファイル番号,行数)
    ecc = np.empty([N_p+N_tr+1, LINE], dtype=float)
    axis = np.empty([N_p+N_tr+1, LINE], dtype=float)
    u = np.empty([N_p+N_tr+1, LINE], dtype=float)
    inc = np.empty([N_p+N_tr+1, LINE], dtype=float)
    Omega = np.empty([N_p+N_tr+1, LINE], dtype=float)
    omega = np.empty([N_p+N_tr+1, LINE], dtype=float)
    r_h = np.empty([N_p+N_tr+1, LINE], dtype=float)
    radius = np.empty([N_p+N_tr+1, LINE], dtype=float)
    mass = np.empty([N_p+N_tr+1, LINE], dtype=float)
    #####
    for n in range(1, N_p+1):
        arr = np.genfromtxt(directory + subdirectory + "Planet%02d.dat" % n, dtype=np.float, delimiter="\t")
        if(LINE - arr.shape[0] > 0):
            print(subnum, n, "padding")
            arr = np.pad(arr, [(0, LINE - arr.shape[0]), (0, 0)], 'constant', constant_values=np.nan)
            # print(arr)

        time[n, :] = arr[:, 0]
        ecc[n, :] = arr[:, 1]
        axis[n, :] = arr[:, 2]
        # u[n, :] = arr[:, 3]
        inc[n, :] = arr[:, 4]
        # Omega[n, :] = arr[:, 5]
        # omega[n, :] = arr[:, 6]
        # r_h[n, :] = arr[:, 7]
        # radius[n, :] = arr[:, 8]
        # mass[n, :] = arr[:, 9]
        # print(n, time[n, 1], axis[n, 1], ecc[n, 1], inc[n, 1])

    #####
    for n in range(N_p+1, N_p+N_tr+1):
        arr = np.genfromtxt(directory + subdirectory + "tracer%06d.dat" % (n-N_p), dtype=np.float, delimiter="\t")
        if(LINE - arr.shape[0] > 0):
            print(subnum, n, "padding")
            arr = np.pad(arr, [(0, LINE - arr.shape[0]), (0, 0)], 'constant', constant_values=np.nan)
            # print(arr)

        time[n, :] = arr[:, 0]
        ecc[n, :] = arr[:, 1]
        axis[n, :] = arr[:, 2]
        # u[n, :] = arr[:, 3]
        inc[n, :] = arr[:, 4]
        # Omega[n, :] = arr[:, 5]
        # omega[n, :] = arr[:, 6]
        # r_h[n, :] = arr[:, 7]
        # radius[n, :] = arr[:, 8]
        # mass[n, :] = arr[:, 9]
        # print(n, time[n, 1], axis[n, 1], ecc[n, 1], inc[n, 1])

    #####

    LIST = [i for i in range(LINE) if i < 59 or (i >= 59 and i % 16 == 10)]
    step = 1

    for T in LIST:
        #####
        plt.figure(figsize=(10, 8), dpi=100)
        plt.title(r"${\rm time}: %.3e {\rm yr}$" % time[1, T], fontsize=25)
        if(N_p == 1):
            # plt.title(r"$N_{\rm tr}=1000,M_{\rm tot}=10 {\rm M_{\oplus}},{\rm time}: %.3e {\rm yr}$" % time[1, T], fontsize=18)
            plt.xlim([0.8, 1.2])
        elif(N_p == 3):
            # plt.title(r"$N_{\rm tr}=3000,M_{\rm tot}=30 {\rm M_{\oplus}},{\rm time}: %.3e {\rm yr}$" % time[1, T], fontsize=18)
            plt.xlim([0.75, 1.3])

        # plt.ylim([0, 0.35])
        plt.ylim([0, 0.2])
        plt.xlabel('semi-major axis [AU]', fontsize=30)
        plt.ylabel('ecc', fontsize=30)

        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid(zorder=1)

        plt.contour(x, y, Jacobi(axis[1, T]), colors=["g"], levels=[3.0], linestyles="dashed", linewidths=2, zorder=2)
        if(N_p == 3):
            plt.contour(x, y, Jacobi(axis[2, T]), colors=["g"], levels=[3.0], linestyles="dashed", linewidths=2, zorder=2)
            plt.contour(x, y, Jacobi(axis[3, T]), colors=["g"], levels=[3.0], linestyles="dashed", linewidths=2, zorder=2)
        plt.scatter(axis[N_p+1:, T], ecc[N_p+1:, T], color="b", s=10, label="Tracer", alpha=0.5, zorder=3)
        plt.scatter(axis[1:N_p+1, T], ecc[1:N_p+1, T], color="r", s=50, label="Planet", zorder=4)
        plt.legend(loc="upper left", fontsize=25)
        plt.tight_layout()
        filename = "../image/" + directory + subdirectory + "axis_ecc_T%02d.png" % step
        plt.savefig(filename)
        plt.close()
        #####

        #####
        plt.figure(figsize=(10, 8), dpi=100)
        plt.title(r"${\rm time}: %.3e {\rm yr}$" % time[1, T], fontsize=25)
        if(N_p == 1):
            # plt.title(r"$N_{\rm tr}=1000,M_{\rm tot}=10 {\rm M_{\oplus}},{\rm time}: %.3e {\rm yr}$" % time[1, T], fontsize=18)
            plt.title(r"${\rm time}: %.3e {\rm yr}$" % time[1, T], fontsize=25)
            plt.xlim([0.8, 1.2])
        elif(N_p == 3):
            # plt.title(r"$N_{\rm tr}=3000,M_{\rm tot}=30 {\rm M_{\oplus}},{\rm time}: %.3e {\rm yr}$" % time[1, T], fontsize=18)
            plt.xlim([0.75, 1.3])

        # plt.ylim([0, 0.16])
        plt.ylim([0, 0.1])
        plt.xlabel('semi-major axis [AU]', fontsize=30)
        plt.ylabel('inc [rad]', fontsize=30)

        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.grid(zorder=1)

        plt.scatter(axis[N_p+1:, T], inc[N_p+1:, T], color="b", s=10, label="Tracer", alpha=0.5, zorder=2)
        plt.scatter(axis[1:N_p+1, T], inc[1:N_p+1, T], color="r", s=50, label="Planet", zorder=3)
        plt.legend(loc="upper left", fontsize=25)
        plt.tight_layout()
        filename = "../image/" + directory + subdirectory + "axis_inc_T%02d.png" % step
        plt.savefig(filename)
        plt.close()
        #####
        print(subnum, time[1, T])
        # plt.show()

        step = step + 1
