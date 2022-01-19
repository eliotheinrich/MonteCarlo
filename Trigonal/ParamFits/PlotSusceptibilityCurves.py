import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
from matplotlib import cm
from AnalyzeSusceptibilityCurve import *
import sys
import os

if __name__ == "__main__":
    label = sys.argv[2]
    os.chdir('data')

    folders = [x for x in os.listdir() if x[:len(sys.argv[1])] == sys.argv[1]]
    keyfunc = lambda x: float(x[(x.index('_')+1):])
    folders = sorted(folders, key=keyfunc)
    vals = [keyfunc(x) for x in folders]

    L = len(folders)

    fig, axs = plt.subplots(nrows=3, ncols=L, sharex='col', sharey='row')
    
    for n,folder in enumerate(folders):
        os.chdir(folder)
        if len(os.listdir()) == 2:
            T, X, dX, E, dE = load_susceptibility_data("SusceptibilityCurve.txt")

            if L > 1:
                ax0 = axs[0,n]
                ax1 = axs[1,n]
                ax2 = axs[2,n]
            else:
                ax0 = axs[0]
                ax1 = axs[1]
                ax2 = axs[2]

            plot_susceptibility_curve(T, X, dX, ax0, color='k', linestyle='-')
            ax0.set_ylabel(r'$\chi$(T)', fontsize=16)
            plot_energy_curve(T, E, dE, ax1, color='k', linestyle='-')
            ax1.set_ylabel('u(T)/J', fontsize=16)
            plot_heat_capacity(T, dE, 0, ax2, 'k-')
            ax2.set_ylabel('c(T)/$k_B$', fontsize=16)
            ax2.set_xlabel('T/J', fontsize=16)

        else:
            T, X1, dX1, E1, dE1 = load_susceptibility_data("SusceptibilityCurve1.txt")
            T, X2, dX2, E2, dE2 = load_susceptibility_data("SusceptibilityCurve2.txt")
            T, X3, dX3, E3, dE3 = load_susceptibility_data("SusceptibilityCurve3.txt")

            if L > 1:
                ax0 = axs[0,n]
                ax1 = axs[1,n]
                ax2 = axs[2,n]
            else:
                ax0 = axs[0]
                ax1 = axs[1]
                ax2 = axs[2]

            plot_susceptibility_curve(T, X1, dX1, ax0, color='k')
            plot_susceptibility_curve(T, X2, dX2, ax0, color='r')
            plot_susceptibility_curve(T, X3, dX3, ax0, color='b')
            plot_heat_capacity(T, dE1, 2, ax1, color='k')
            plot_heat_capacity(T, dE2, 2, ax1, color='r')
            plot_heat_capacity(T, dE3, 2, ax1, color='b')
            plot_energy_curve(T, E1, dE1, ax2, color='k', label=r'$B \parallel a$')
            plot_energy_curve(T, E2, dE2, ax2, color='r', label=r'$B \parallel ab$')
            plot_energy_curve(T, E3, dE3, ax2, color='b', label=r'$B \parallel c$')

            ax0.set_title(label + " = " + str(vals[n]))
            ax2.set_xlabel(r"$T/J$", fontsize=20)

        os.chdir('..')

    if L > 1:
        ax0 = axs[0,0]
        ax1 = axs[1,0]
        ax2 = axs[2,0]
    else:
        ax0 = axs[0]
        ax1 = axs[1]
        ax2 = axs[2]

    ax2.legend()
    ax0.set_ylabel(r'$\chi(T)$', fontsize=20)
    ax1.set_ylabel(r'$c(T)$', fontsize=20)
    ax2.set_ylabel(r'$u(T)$', fontsize=20)

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.show()

        




