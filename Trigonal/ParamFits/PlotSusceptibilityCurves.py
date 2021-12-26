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

        plot_susceptibility_curve(T, X1, X2, X3, dX1, dX2, dX3, ax0)
        plot_heat_capacity(T, dE1, dE2, dE3, ax1, nmin=2)
        plot_energy_curve(T, E1, E2, E3, dE1, dE2, dE3, ax2)

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

    ax0.legend()
    ax0.set_ylabel(r'$\chi(T)$', fontsize=20)
    ax1.set_ylabel(r'$c(T)$', fontsize=20)
    ax2.set_ylabel(r'$u(T)$', fontsize=20)

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.show()

        




