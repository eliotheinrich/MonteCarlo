import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
from matplotlib import cm
import sys
import os

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
            T[i] = float(data[0])
            M[i] = float(data[1])
            dM[i] = float(data[2])
            E[i] = float(data[3])
            dE[i] = float(data[4])


    return T[1:], M[1:], dM[1:], E[1:], dE[1:]


def plot_magnetization_curve(T, M, dM, ax = None):
    ax_none = ax is None
    if ax_none:
        ax = plt.gca()

    ax.errorbar(T, M, yerr = dM, color = 'k', linestyle = '-')
    if ax_none:
        ax.set_xlabel(r'$T$ (K)', fontsize=20)
        ax.set_ylabel(r'$\chi$', fontsize=20)

def plot_energy_curve(T, E, dE, ax = None):
    ax_none = ax is None
    if ax_none:
        ax = plt.gca()

    ax.errorbar(T, E, yerr = dE, color = 'k', linestyle = '-')
    if ax_none:
        ax.set_xlabel(r'$T$ (K)', fontsize=20)
        ax.set_ylabel(r'$U(T)$', fontsize=20)

def plot_heat_capacity(T, dE, ax = None, nmin = 0):
    ax_none = ax is None
    if ax_none:
        ax = plt.gca()

    ax.plot(T[nmin:], ((dE/T)**2)[nmin:], 'k-')
    if ax_none:
        ax.set_xlabel(r'$T$ (K)', fontsize=20)
        ax.set_ylabel(r'$E$', fontsize=20)

if __name__ == "__main__":
    os.chdir(sys.argv[1])

    T, M, dM, E, dE = load_magnetization_data("MagnetizationCurve.txt")


    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
    plot_magnetization_curve(T, M, dM, axs[0])
    axs[0].set_ylabel('M(T)', fontsize=16)
    plot_energy_curve(T, E, dE, axs[1])
    axs[1].set_ylabel('u(T)/J', fontsize=16)
#    plot_heat_capacity(T, dE, axs[2])
#    axs[2].set_ylabel('c(T)', fontsize=16)
    axs[1].set_xlabel('T/J', fontsize=16)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()





