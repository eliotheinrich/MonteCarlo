import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
from matplotlib import cm
import sys
import os


kB = 0.0816

def load_magnetization_data(filename):
    with open(filename) as f:
        s = f.readline().split('\t')
        res = int(s[0])
    
    T = np.zeros(res)
    M = np.zeros(res)
    dM = np.zeros(res)
    E = np.zeros(res)
    dE = np.zeros(res)

    with open(filename) as f:
        line = f.readline()
        for i in range(res):
            line = f.readline()

            data = np.array([x.strip('\n') for x in line.split('\t')])
            T[i] = float(data[0])/kB
            (M[i], dM[i]) = (float(x) for x in data[1].split(','))
            (E[i], dE[i]) = (float(x) for x in data[2].split(','))

    return T, M, dM, E, dE


def curie_fit(T, M):
    def curie_weiss(T, C, Tc):
        return (T - Tc)/C

    N = len(M)//2

    popt, _ = spo.curve_fit(curie_weiss, T[N:], 1/M[N:])
    (C, Tc) = popt

    return C, Tc

def plot_magnetization_curve(T, M, dM, ax, **kwargs):
    ax.errorbar(T, M, yerr = dM, **kwargs)

def plot_energy_curve(T, E, dE, ax, *args, **kwargs):
    ax.errorbar(T, E, yerr = dE, *args, **kwargs)

def plot_heat_capacity(T, dE, nmin, ax, *args, **kwargs):
    ax.plot(T[nmin:], ((dE/T)**2)[nmin:]/kB, *args, **kwargs)

if __name__ == "__main__":
    os.chdir(sys.argv[1])

    if len(os.listdir()) == 2:
        T, M, dM, E, dE = load_magnetization_data("MagnetizationCurve.txt")

        fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True)
        plot_magnetization_curve(T, M, dM, axs[0], color='k', linestyle='-')
        axs[0].set_ylabel('M(T)', fontsize=16)
        plot_energy_curve(T, E, dE, axs[1], color='k', linestyle='-')
        axs[1].set_ylabel('u(T) (meV)', fontsize=16)
        plot_heat_capacity(T, dE, 4, axs[2], 'k-')
        axs[2].set_ylabel('c(T)/$k_B$', fontsize=16)
        axs[2].set_xlabel('T (meV)', fontsize=16)

        plt.subplots_adjust(wspace=0, hspace=0)
        plt.show()

    else:
        T, M1, dM1, E1, dE1 = load_magnetization_data("MagnetizationCurve1.txt")
        T, M2, dM2, E2, dE2 = load_magnetization_data("MagnetizationCurve2.txt")
        T, M3, dM3, E3, dE3 = load_magnetization_data("MagnetizationCurve3.txt")

        #C, _, _, _ = curie_fit(T, M1, M2, M3)

        kwargs1 = {'color':'k', 'linestyle':'-', 'label':r'$B \parallel ab$'}
        kwargs2 = {'color':'r', 'linestyle':'-', 'label':r'$B \parallel a$'}
        kwargs3 = {'color':'b', 'linestyle':'-', 'label':r'$B \parallel c$'}
        ax = plt.gca()
        plot_magnetization_curve(T, M, dM, ax, kwargs1)
        plot_magnetization_curve(T, M, dM, ax, kwargs2)
        plot_magnetization_curve(T, M, dM, ax, kwargs3)
        plt.show()




