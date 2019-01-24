#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
import numpy as np



# directory = "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-2_nofrag_acc/"
directory = "Ntr1E2_t1E8_dtlog_Mtot3E-7_Mmax5E-15_ecc1E-1_nofrag_acc/"

N_p = 1
N_tr = 100

if(N_p == 1):
    # LINE = 42  # 10000yr
    # LINE = 299  # 1E8yr
    SUBDIR_NUM = 1
elif(N_p == 3):
    # LINE = 34  # 1000yr
    SUBDIR_NUM = 13


for subnum in range(1, SUBDIR_NUM+1):
    subdirectory = "rand%02d/" % subnum
    mass_p = 3.0E-6
    mass_tr = 3.0E-9
    total_mass = 1.0 + 3.0E-6 + 3.0E-7  # 中心星込みの総質量
    dE_correct = 0.0
    dLx = 0.0
    dLy = 0.0
    dLz = 0.0


    """
    with open(directory + subdirectory + "ENERGY.dat", "r") as f:
        collision_time = [float(line.split("\t")[1][:12]) for line in f.readlines() if line.find("collision") >= 0]  # collision を含む行を検索し、TABで区切り、2列目の文字列の前から12文字目までを小数で表示
        print(collision_time)
    """

    with open(directory + subdirectory + "Collision_time.dat", "w") as f:
        print("#1:time",
              "2:N_tr",
              "3:total_energy",
              "4:rel_total_energy",
              "5:dE_correct",
              "6:rel_dE_correct",
              "7:energy_p",
              "8:rel_energy_p",
              "9:total_angmom.",
              "10:rel_total_angmom.",
              "11:abs_L_p",
              "12:rel_abs_L_p",
              "13:total_Lx",
              "14:rel_total_Lx",
              "15:dLx",
              "16:rel_dLx",
              "17:total_Ly",
              "18:rel_total_Ly",
              "19:dLy",
              "20:rel_dLy",
              "21:total_Lz",
              "22:rel_total_Lz",
              "23:dLz",
              "24:rel_dLz",
              file=f, sep="\t")

        print(0.0, N_tr, file=f, sep="\t")


    for n_col in range(1, N_tr+1):


        datafile = directory + subdirectory + "Collision_%03d.dat" % n_col
        if os.path.exists(datafile):
            arr = np.genfromtxt(datafile, dtype=np.float, delimiter="\t")
            # print(arr.shape)


            time = np.empty([SUBDIR_NUM, N_p+N_tr+1-n_col], dtype=float)  # (ファイル番号,行数)
            i_tr = np.empty([SUBDIR_NUM, N_p+N_tr+1-n_col], dtype=int)
            x = np.empty([SUBDIR_NUM, N_p+N_tr+1-n_col], dtype=float)
            y = np.empty([SUBDIR_NUM, N_p+N_tr+1-n_col], dtype=float)
            z = np.empty([SUBDIR_NUM, N_p+N_tr+1-n_col], dtype=float)
            abs_r = np.empty([SUBDIR_NUM, N_p+N_tr+1-n_col], dtype=float)
            vx = np.empty([SUBDIR_NUM, N_p+N_tr+1-n_col], dtype=float)
            vy = np.empty([SUBDIR_NUM, N_p+N_tr+1-n_col], dtype=float)
            vz = np.empty([SUBDIR_NUM, N_p+N_tr+1-n_col], dtype=float)
            abs_v = np.empty([SUBDIR_NUM, N_p+N_tr+1-n_col], dtype=float)
            radius = np.empty([SUBDIR_NUM, N_p+N_tr+1-n_col], dtype=float)
            temp = np.zeros([N_p+N_tr-n_col], dtype=float)


            time[subnum-1, :] = arr[:, 0]
            i_tr[subnum-1, :] = arr[:, 1]
            x[subnum-1, :] = arr[:, 2]
            y[subnum-1, :] = arr[:, 3]
            z[subnum-1, :] = arr[:, 4]
            abs_r[subnum-1, :] = arr[:, 5]
            vx[subnum-1, :] = arr[:, 7]
            vy[subnum-1, :] = arr[:, 8]
            vz[subnum-1, :] = arr[:, 9]
            abs_v[subnum-1, :] = arr[:, 10]
            radius[subnum-1, :] = arr[:, 12]


            temp = mass_tr * x[subnum-1, 1:]
            x_G = (mass_p * x[subnum-1, 0] + temp.sum()) / total_mass
            temp = np.zeros([N_p+N_tr-n_col], dtype=float)

            temp = mass_tr * y[subnum-1, 1:]
            y_G = (mass_p * y[subnum-1, 0] + temp.sum()) / total_mass
            temp = np.zeros([N_p+N_tr-n_col], dtype=float)

            temp = mass_tr * z[subnum-1, 1:]
            z_G = (mass_p * z[subnum-1, 0] + temp.sum()) / total_mass
            temp = np.zeros([N_p+N_tr-n_col], dtype=float)


            temp = mass_tr * vx[subnum-1, 1:]
            vx_G = (mass_p * vx[subnum-1, 0] + temp.sum()) / total_mass
            temp = np.zeros([N_p+N_tr-n_col], dtype=float)

            temp = mass_tr * vy[subnum-1, 1:]
            vy_G = (mass_p * vy[subnum-1, 0] + temp.sum()) / total_mass
            temp = np.zeros([N_p+N_tr-n_col], dtype=float)

            temp = mass_tr * vz[subnum-1, 1:]
            vz_G = (mass_p * vz[subnum-1, 0] + temp.sum()) / total_mass
            temp = np.zeros([N_p+N_tr-n_col], dtype=float)

            ##########
            energy_p = 0.5 * mass_p * ((vx[subnum-1, 0] - vx_G)**2 + (vy[subnum-1, 0] - vy_G)**2 + (vz[subnum-1, 0] - vz_G)**2) - mass_p / abs_r[subnum-1, 0]
            energy_i = 0.5 * mass_tr * ((vx[subnum-1, 1:] - vx_G)**2 + (vy[subnum-1, 1:] - vy_G)**2 + (vz[subnum-1, 1:] - vz_G)**2)  # indirect term 有りの tracer の運動エネルギー
            energy_i += - mass_tr / abs_r[subnum-1, 1:]  # planet（0）を除く、中心星からの重力エネルギー
            energy_rel = - mass_p * mass_tr / np.sqrt((x[subnum-1, 1:] - x[subnum-1, 0])**2 + (y[subnum-1, 1:] - y[subnum-1, 0])**2 + (z[subnum-1, 1:] - z[subnum-1, 0])**2)  # tracer と planet（0）との相互重力エネルギー
            total_energy = 0.5 * (vx_G**2 + vy_G**2 + vz_G**2) + energy_p + energy_i.sum() + energy_rel.sum() + dE_correct
            ##########


            ##########
            Lx_p = mass_p * (y[subnum-1, 0] * vz[subnum-1, 0] - z[subnum-1, 0] * vy[subnum-1, 0])
            Ly_p = mass_p * (z[subnum-1, 0] * vx[subnum-1, 0] - x[subnum-1, 0] * vz[subnum-1, 0])
            Lz_p = mass_p * (x[subnum-1, 0] * vy[subnum-1, 0] - y[subnum-1, 0] * vx[subnum-1, 0])
            Lx_i = mass_tr * (y[subnum-1, 1:] * vz[subnum-1, 1:] - z[subnum-1, 1:] * vy[subnum-1, 1:])
            Ly_i = mass_tr * (z[subnum-1, 1:] * vx[subnum-1, 1:] - x[subnum-1, 1:] * vz[subnum-1, 1:])
            Lz_i = mass_tr * (x[subnum-1, 1:] * vy[subnum-1, 1:] - y[subnum-1, 1:] * vx[subnum-1, 1:])
            total_Lx = Lx_p + Lx_i.sum() - dLx
            total_Ly = Ly_p + Ly_i.sum() - dLy
            total_Lz = Lz_p + Lz_i.sum() - dLz
            abs_L_p = np.sqrt(Lx_p**2 + Ly_p**2 + Lz_p**2)
            total_angularmomentum = np.sqrt(total_Lx**2 + total_Ly**2 + total_Lz**2)
            ##########

            if n_col == 1:
                total_energy_0 = total_energy
                energy_p_0 = energy_p
                total_Lx_0 = total_Lx
                total_Ly_0 = total_Ly
                total_Lz_0 = total_Lz
                total_angularmomentum_0 = total_angularmomentum
                abs_L_p_0 = abs_L_p

            # print(x_G, y_G, z_G, vx_G, vy_G, vz_G)
            # print(n_col, temp, temp.sum())

            print(time[subnum-1, 0])
            # print(n_col, total_energy, abs((total_energy - total_energy_0)/total_energy_0), energy_p, abs((energy_p - energy_p_0)/energy_p_0))
            print(n_col, total_angularmomentum, abs((total_angularmomentum - total_angularmomentum_0)/total_angularmomentum_0), abs_L_p, abs((abs_L_p - abs_L_p_0)/abs_L_p_0))


            with open(datafile, "r") as f:
                j_col = int(f.readline().split("\t")[1][6:])
                # print(j_col)


            ##########
            x_g12 = (mass_p * x[subnum-1, 0] + mass_tr * x[subnum-1, j_col-1]) / (mass_p + mass_tr)
            y_g12 = (mass_p * y[subnum-1, 0] + mass_tr * y[subnum-1, j_col-1]) / (mass_p + mass_tr)
            z_g12 = (mass_p * z[subnum-1, 0] + mass_tr * z[subnum-1, j_col-1]) / (mass_p + mass_tr)

            r_g12 = np.sqrt(x_g12**2 + y_g12**2 + z_g12**2)  # 中心星から衝突した2天体の重心までの距離

            dE_heat = - 0.5 * mass_p * mass_tr / (mass_p + mass_tr) * ((vx[subnum-1, j_col-1] - vx[subnum-1, 0])**2 + (vy[subnum-1, j_col-1] - vy[subnum-1, 0])**2 + (vz[subnum-1, j_col-1] - vz[subnum-1, 0])**2)

            dE_grav = mass_p * mass_tr / np.sqrt((x[subnum-1, j_col-1] - x[subnum-1, 0])**2 + (y[subnum-1, j_col-1] - y[subnum-1, 0])**2 + (z[subnum-1, j_col-1] - z[subnum-1, 0])**2)

            dE_c = - (mass_p + mass_tr) / r_g12 + mass_p / abs_r[subnum-1, 0] + mass_tr / abs_r[subnum-1, j_col-1]

            dE_correct += - dE_heat - dE_grav - dE_c

            # print(n_col, dE_heat / total_energy_0, dE_grav / total_energy_0, dE_c / total_energy_0, dE_correct / total_energy_0)
            ##########


            ##########
            vx_g12 = (mass_p * vx[subnum-1, 0] + mass_tr * vx[subnum-1, j_col-1]) / (mass_p + mass_tr)
            vy_g12 = (mass_p * vy[subnum-1, 0] + mass_tr * vy[subnum-1, j_col-1]) / (mass_p + mass_tr)
            vz_g12 = (mass_p * vz[subnum-1, 0] + mass_tr * vz[subnum-1, j_col-1]) / (mass_p + mass_tr)

            dLx += (mass_p + mass_tr) * (y_g12 * vz_g12 - z_g12 * vy_g12) - mass_p * (y[subnum-1, 0] * vz[subnum-1, 0] - z[subnum-1, 0] * vy[subnum-1, 0]) - mass_tr * (y[subnum-1, j_col-1] * vz[subnum-1, j_col-1] - z[subnum-1, j_col-1] * vy[subnum-1, j_col-1])
            dLy += (mass_p + mass_tr) * (z_g12 * vx_g12 - x_g12 * vz_g12) - mass_p * (z[subnum-1, 0] * vx[subnum-1, 0] - x[subnum-1, 0] * vz[subnum-1, 0]) - mass_tr * (z[subnum-1, j_col-1] * vx[subnum-1, j_col-1] - x[subnum-1, j_col-1] * vz[subnum-1, j_col-1])
            dLz += (mass_p + mass_tr) * (x_g12 * vy_g12 - y_g12 * vx_g12) - mass_p * (x[subnum-1, 0] * vy[subnum-1, 0] - y[subnum-1, 0] * vx[subnum-1, 0]) - mass_tr * (x[subnum-1, j_col-1] * vy[subnum-1, j_col-1] - y[subnum-1, j_col-1] * vx[subnum-1, j_col-1])
            print(n_col, dLx, dLy, dLz)
            print("")
            ##########


            with open(directory + subdirectory + "Collision_time.dat", "a") as f:
                print(time[subnum-1, 0],
                      N_tr - n_col,
                      total_energy,
                      (total_energy - total_energy_0) / abs(total_energy_0),
                      dE_correct,
                      dE_correct / abs(total_energy_0),
                      energy_p,
                      (energy_p - energy_p_0) / abs(energy_p_0),
                      total_angularmomentum,
                      (total_angularmomentum - total_angularmomentum_0) / abs(total_angularmomentum_0),
                      abs_L_p,
                      (abs_L_p - abs_L_p_0) / abs(abs_L_p_0),
                      total_Lx,
                      (total_Lx - total_Lx_0) / abs(total_Lx_0),
                      dLx,
                      dLx / abs(total_Lx_0),
                      total_Ly,
                      (total_Ly - total_Ly_0) / abs(total_Ly_0),
                      dLy,
                      dLy / abs(total_Ly_0),
                      total_Lz,
                      (total_Lz - total_Lz_0) / abs(total_Lz_0),
                      dLz,
                      dLz / abs(total_Lz_0),
                      file=f, sep="\t")


            mass_p += mass_tr
