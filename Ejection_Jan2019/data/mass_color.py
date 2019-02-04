#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from scipy.interpolate import griddata

def delta_axis(ratio):
    return (1.0/ratio + 0.5*np.cbrt(2.0E-6))/(1.0/ratio - 0.5*np.cbrt(2.0E-6))


######################################################################
# directory = "N25_t1E5_dtlog_R0.1to0.25_Theta30_EscVelo1.1_Col001_Ntr100_Mmax3.5E-10/"
# N_p = 24
# PLANET_LIST = [1, 2, 3, 4, 5, 6, 25, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
# COLLISION_PLANET = 6
# LINE = 42  # 1E4yr


# directory = "N18_t1E5_dtlog_R0.1to0.25_Theta30_EscVelo1.1_Col001_Ntr100_Mmax1.2E-9/"
# N_p = 17
# PLANET_LIST = [1, 2, 3, 4, 5, 6, 7, 18, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
# COLLISION_PLANET = 7
# LINE = 42  # 1E4yr
# LINE = 44  # 1.77E4yr


# directory = "N25_t1E5_dtlog_R0.1to0.25_Theta30_EscVelo1.1_Col020_Ntr100_Mmax7.2E-10/"
# N_p = 5
# PLANET_LIST = [1, 2, 23, 21, 25]
# COLLISION_PLANET = 2
# LINE = 50  # 1E5yr


directory = "N18_t1E5_dtlog_R0.1to0.25_Theta30_EscVelo1.1_Col014_Ntr100_Mmax1.5E-9/"
N_p = 4
PLANET_LIST = [1, 17, 16, 14]
COLLISION_PLANET = 4
LINE = 50  # 1E5yr



N_tr = 100




SUBDIR_NUM = 1



for subnum in range(1, SUBDIR_NUM+1):
# for subnum in range(1, 2):
    subdirectory = "rand%02d/" % subnum

    time = np.empty([N_p+N_tr, LINE], dtype=float)  # (ファイル番号,行数)
    num = np.empty([N_p+N_tr, LINE], dtype=float)
    x = np.empty([N_p+N_tr, LINE], dtype=float)
    y = np.empty([N_p+N_tr, LINE], dtype=float)
    z = np.empty([N_p+N_tr, LINE], dtype=float)
    r_3d = np.empty([N_p+N_tr, LINE], dtype=float)
    r_2d = np.empty([N_p+N_tr, LINE], dtype=float)
    mass = np.empty([N_p+N_tr, LINE], dtype=float)
    sigma = np.empty([N_p+N_tr, LINE], dtype=float)

    ecc = np.empty([N_p+1, LINE], dtype=float)
    axis = np.empty([N_p+1, LINE], dtype=float)
    u = np.empty([N_p+1, LINE], dtype=float)
    inc = np.empty([N_p+1, LINE], dtype=float)
    Omega = np.empty([N_p+1, LINE], dtype=float)
    omega = np.empty([N_p+1, LINE], dtype=float)
    Px = np.empty([N_p+1, LINE], dtype=float)
    Py = np.empty([N_p+1, LINE], dtype=float)
    Pz = np.empty([N_p+1, LINE], dtype=float)
    Qx = np.empty([N_p+1, LINE], dtype=float)
    Qy = np.empty([N_p+1, LINE], dtype=float)
    Qz = np.empty([N_p+1, LINE], dtype=float)
    orbit_X = np.empty([N_p+1, LINE, 100], dtype=float)
    orbit_Y = np.empty([N_p+1, LINE, 100], dtype=float)
    orbit_Z = np.empty([N_p+1, LINE, 100], dtype=float)

    # radius = np.empty([N_p+N_tr], dtype=float)
    # theta = np.empty([N_p+N_tr], dtype=float)

    for T in range(LINE):
        arr = np.genfromtxt(directory + subdirectory + "Posi_Mass_%03d.dat" % T, dtype=np.float, delimiter="\t")
        if(N_p+N_tr - arr.shape[0] > 0):
            print(T, "padding")
            arr = np.pad(arr, [(0, N_p+N_tr - arr.shape[0]), (0, 0)], 'constant', constant_values=np.NaN)

        time[:, T] = arr[:, 0]
        num[:, T] = arr[:, 1]
        x[:, T] = arr[:, 2]
        y[:, T] = arr[:, 3]
        z[:, T] = arr[:, 4]
        r_3d[:, T] = arr[:, 5]
        r_2d[:, T] = arr[:, 6]
        mass[:, T] = arr[:, 7]
        sigma[:, T] = arr[:, 10]
        print(T, int(num[0, T]), time[0, T], x[0, T])


    for n in range(1, N_p+1):
        arr2 = np.genfromtxt(directory + subdirectory + "Planet%02d.dat" % PLANET_LIST[n-1], dtype=np.float, delimiter="\t")

        ecc[n, :] = arr2[:, 1]
        axis[n, :] = arr2[:, 2]
        u[n, :] = arr2[:, 3]
        inc[n, :] = arr2[:, 4]
        Omega[n, :] = arr2[:, 5]
        omega[n, :] = arr2[:, 6]
        print(n, PLANET_LIST[n-1], time[n, 1], axis[n, 1], ecc[n, 1], inc[n, 1])

        Px[n, :] = np.cos(omega[n, :]) * np.cos(Omega[n, :]) - np.sin(omega[n, :]) * np.sin(Omega[n, :]) * np.cos(inc[n, :])
        Qx[n, :] = - np.sin(omega[n, :]) * np.cos(Omega[n, :]) - np.cos(omega[n, :]) * np.sin(Omega[n, :]) * np.cos(inc[n, :])

        Py[n, :] = np.cos(omega[n, :]) * np.sin(Omega[n, :]) + np.sin(omega[n, :]) * np.cos(Omega[n, :]) * np.cos(inc[n, :])
        Qy[n, :] = - np.sin(omega[n, :]) * np.sin(Omega[n, :]) + np.cos(omega[n, :]) * np.cos(Omega[n, :]) * np.cos(inc[n, :])

        Pz[n, :] = np.sin(omega[n, :]) * np.sin(inc[n, :])
        Qz[n, :] = np.cos(omega[n, :]) * np.sin(inc[n, :])


        #########


    for T in range(LINE):

        fig = plt.figure(figsize=(8, 6), dpi=100)
        ax = fig.add_subplot(1, 1, 1, aspect='equal')
        # ax.set_xlim([-0.6, 0.6])
        # ax.set_ylim([-0.6, 0.6])
        # ax.set_xlim([-1.0, 1.0])
        # ax.set_ylim([-1.0, 1.0])
        # ax.set_xlim([-1.5, 1.5])
        # ax.set_ylim([-1.5, 1.5])
        ax.set_xlim([-2.0, 2.0])
        ax.set_ylim([-2.0, 2.0])

        ax.set_xlabel('x [AU]', fontsize=20)
        ax.set_ylabel('y [AU]', fontsize=20)
        ax.title.set_text('%.3e [yr]' % time[1, T])



        for n in range(1, N_p+1):
            orbit_X[n, T, :] = axis[n, T] * Px[n, T] * (np.cos(np.linspace(0.0, 2.0*np.pi, 100)) - ecc[n, T]) + axis[n, T] * np.sqrt(1.0 - ecc[n, T]*ecc[n, T]) * Qx[n, T] * np.sin(np.linspace(0.0, 2.0*np.pi, 100))
            orbit_Y[n, T, :] = axis[n, T] * Py[n, T] * (np.cos(np.linspace(0.0, 2.0*np.pi, 100)) - ecc[n, T]) + axis[n, T] * np.sqrt(1.0 - ecc[n, T]*ecc[n, T]) * Qy[n, T] * np.sin(np.linspace(0.0, 2.0*np.pi, 100))
            orbit_Z[n, T, :] = axis[n, T] * Pz[n, T] * (np.cos(np.linspace(0.0, 2.0*np.pi, 100)) - ecc[n, T]) + axis[n, T] * np.sqrt(1.0 - ecc[n, T]*ecc[n, T]) * Qz[n, T] * np.sin(np.linspace(0.0, 2.0*np.pi, 100))

            if(n==COLLISION_PLANET):
                ax.plot(orbit_X[COLLISION_PLANET, T, :], orbit_Y[COLLISION_PLANET, T, :], color="k", alpha=0.7, linestyle="dashed")
                ax.scatter(x[COLLISION_PLANET-1, T], y[COLLISION_PLANET-1, T], color="k", marker="+", s=50)
            else:
                ax.plot(orbit_X[n, T, :], orbit_Y[n, T, :], color="k", alpha=0.2)
                ax.scatter(x[n-1, T], y[n-1, T], color="k", marker="+", alpha=0.2, s=50)




        fig.subplots_adjust(right=0.85)

        cm = plt.cm.get_cmap("rainbow")
        cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])

        # im = ax.scatter(x[N_p:, T], y[N_p:, T], s=10, c=mass[N_p:, T] * 2E33 / 6E27, cmap=cm, vmax=M_MAX * 2E33 / 6E27)

        im = ax.scatter(x[N_p:, T], y[N_p:, T], s=10, c=sigma[N_p:, T] * 2.0E33/1.5E13/1.5E13, cmap=cm, norm=LogNorm(0.1, sigma[N_p:, 0].max() * 2.0E33/1.5E13/1.5E13))

        ax.scatter(0, 0, c='k', marker='*', s=100)

        fig.colorbar(im, cax=cbar_ax).set_label(r'$\Sigma\ [\rm g/cm^2]$', fontsize=20)
        # fig.colorbar(im, cax=cbar_ax).set_label(r'mass$\ [\rm M_{\oplus}]$', fontsize=20)


        filename = "../image/" + directory + subdirectory + "Sigma_T%02d.png" % T
        plt.savefig(filename, format="png", dpi=100)
        # plt.show()
        plt.close()



        """
        ax = fig.add_subplot(1, 1, 1, aspect='equal')
        ax.set_xlim([-2, 2])
        ax.set_ylim([-2, 2])
        ax.set_xlabel('x [AU]', fontsize=20)
        ax.set_ylabel('y [AU]', fontsize=20)
        """

        """
        ax = fig.add_subplot(111, projection="polar")
        # ax = fig.add_subplot(111)
        ax.title.set_text('%.3e [yr]' % time[1, T])
        ax.title.set_fontsize(15)

        fig.subplots_adjust(right=0.8)

        # cm = plt.cm.get_cmap('rainbow')

        RADIUS = np.linspace(0.0, 1.5, 500)
        THETA = np.linspace(-np.pi, np.pi, 500)

        radius = np.sqrt(x[:, T]**2 + y[:, T]**2)

        if(N_p == 1):
            radius = np.append(radius, np.random.uniform(0.0, 1.0/delta_axis(5.0), 10000))
            radius = np.append(radius, np.random.uniform(delta_axis(5.0), 1.5, 5000))
        elif(N_p == 3):
            radius = np.append(radius, np.random.uniform(0.0, 1.0/delta_axis(10.0)/delta_axis(5.0), 10000))
            radius = np.append(radius, np.random.uniform(delta_axis(10.0)*delta_axis(5.0), 1.5, 5000))

        theta = np.arctan2(y[:, T], x[:, T])
        theta = np.append(theta, np.random.uniform(-np.pi-0.1, np.pi+0.1, 15000))

        sigma_grid = sigma[:, T]*8888888.88888889
        sigma_grid = np.append(sigma_grid, [0.0 for i in range(15000)])

        RADIUS, THETA = np.meshgrid(RADIUS, THETA)

        Z = griddata((radius[N_p:], theta[N_p:]), sigma_grid[N_p:], (RADIUS, THETA), method="linear")
        Z[np.isnan(Z)] = 0.0
        Z_convert = np.where(Z < 0.0, 0.0, Z)

        # im = ax.scatter(theta[3:], radius[3:], s=1, c=mass[3:, T] * 2E8, vmin=1.0, vmax=6.0, cmap=cm, label='Tracer')

        # im = ax.contourf(THETA, RADIUS, Z, levels=np.linspace(0.0, 3.0, 100), cmap=cm)
        # im = ax.contourf(THETA, RADIUS, Z_convert, cmap="jet")
        im = ax.contourf(THETA, RADIUS, Z_convert, levels=np.linspace(0.0, np.amax(Z), 100), cmap="jet")
        # im = ax.contourf(THETA, RADIUS, Z_convert, norm=LogNorm(), levels=[1E-1, 1E0, 1E1, 1E2, 1E3], cmap="jet")

        # ax.scatter(0, 0, c='k', marker='*', s=50, label='Star')
        if(N_p == 4):
            ax.scatter(theta[0], radius[0], c='k', marker='1', s=50, label='Planet1')
            ax.scatter(theta[1], radius[1], c='k', marker='2', s=50, label='Planet2')
            ax.scatter(theta[2], radius[2], c='k', marker='3', s=50, label='Planet3')
            ax.scatter(theta[3], radius[3], c='k', marker='4', s=50, label='Planet4')
        elif(N_p == 5):
            ax.scatter(theta[0], radius[0], c='k', marker='1', s=50, label='Planet1')
            ax.scatter(theta[1], radius[1], c='k', marker='2', s=50, label='Planet2')
            ax.scatter(theta[2], radius[2], c='k', marker='3', s=50, label='Planet3')
            ax.scatter(theta[3], radius[3], c='k', marker='4', s=50, label='Planet4')
            ax.scatter(theta[4], radius[4], c='k', marker='5', s=50, label='Planet5')

        ax.set_rlim(0, 1.5)
        ax.set_rticks([0.5, 1.0, 1.5])
        ax.set_thetalim(0, 2*np.pi)
        ax.set_xticks(np.pi/180. * np.linspace(0, 360, 4, endpoint=False))
        cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
        fig.colorbar(im, cax=cbar_ax, ticks=np.linspace(0.0, np.amax(Z), 10))
        # fig.colorbar(im, cax=cbar_ax)
        # cbar_ax.set_ylabel(r'mass $[10^{25} \rm g]$', fontsize=20)
        cbar_ax.set_ylabel(r'$\Sigma \ \rm [g/cm^2]$', fontsize=20)

        # plt.tight_layout()
        ax.legend()
        filename = "../image/" + directory + subdirectory + "Sigma_T%02d.png" % T
        plt.savefig(filename, format="png", dpi=100)
        # plt.show()
        plt.close()
        """
