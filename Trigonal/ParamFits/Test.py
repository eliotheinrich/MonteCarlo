import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
from matplotlib import cm
from AnalyzeSusceptibilityCurve import *
import sys
import os

if __name__ == "__main__":
    os.chdir('data')

    folders = ['mX_1', 'mX_2']
    keyfunc = lambda x: float(x[(x.index('_')+1):])
    folders = sorted(folders, key=keyfunc)
    vals = [keyfunc(x) for x in folders]

    fig, axs = plt.subplots(nrows=3, ncols=3, sharex='col')

    X3 = []
    dX3 = []
    E3 = []
    dE3 = []

    for n,f in enumerate(folders):
        os.chdir(f)
        T, X, dX, E, dE = load_susceptibility_data("SusceptibilityCurve.txt")
        X3.append(X)
        dX3.append(dX)
        E3.append(E)
        dE3.append(dE)

        ax0 = axs[0,n]
        ax1 = axs[1,n]
        ax2 = axs[2,n]

        plot_susceptibility_curve(T, X, dX, ax0, color='k', linestyle='-')
        ax0.set_ylabel(r'$\chi$(T)', fontsize=16)
        plot_energy_curve(T, E, dE, ax1, color='k', linestyle='-')
        ax1.set_ylabel('u(T)/J', fontsize=16)
        plot_heat_capacity(T, dE, 0, ax2, 'k-')
        ax2.set_ylabel('c(T)/$k_B$', fontsize=16)
        ax2.set_xlabel('T/J', fontsize=16)

        os.chdir('..')

    ax0 = axs[0,2]
    ax1 = axs[1,2]
    ax2 = axs[2,2]

    plot_susceptibility_curve(T, X3[0], dX3[0], ax0, color='k', linestyle='-')
    plot_susceptibility_curve(T, X3[1], dX3[1], ax0, color='r', linestyle='-')
    ax0.set_ylabel(r'$\chi$(T)', fontsize=16)
    plot_energy_curve(T, E3[0], dE3[0], ax1, color='k', linestyle='-')
    plot_energy_curve(T, E3[1], dE3[1], ax1, color='r', linestyle='-')
    ax1.set_ylabel('u(T)/J', fontsize=16)
    plot_heat_capacity(T, dE3[0], 0, ax2, 'k-')
    plot_heat_capacity(T, dE3[1], 0, ax2, 'r-')
    ax2.set_ylabel('c(T)/$k_B$', fontsize=16)
    ax2.set_xlabel('T/J', fontsize=16)

    ax2.legend()
    ax0.set_ylabel(r'$\chi(T)$', fontsize=20)
    ax1.set_ylabel(r'$c(T)$', fontsize=20)
    ax2.set_ylabel(r'$u(T)$', fontsize=20)

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.show()

        




