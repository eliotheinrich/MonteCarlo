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

    fig, axs = plt.subplots(nrows=1, ncols=L, sharex='col', sharey='row')
    
    for n,folder in enumerate(folders):
        os.chdir(folder)
        if len(os.listdir()) == 2:
            T, X, dX, E, dE = load_susceptibility_data("SusceptibilityCurve.txt")

            if L > 1:
                ax0 = axs[n]
            else:
                ax0 = axs[n]

            plot_susceptibility_curve(T, X, dX, ax0, color='k', linestyle='-')
            ax0.set_ylabel(r'$\chi$(T) ($\mu_B/T$)', fontsize=16)
            ax0.set_xlabel('T (meV)', fontsize=16)

        else:
            T, X1, dX1, E1, dE1 = load_susceptibility_data("SusceptibilityCurve1.txt")
            T, X2, dX2, E2, dE2 = load_susceptibility_data("SusceptibilityCurve2.txt")
            T, X3, dX3, E3, dE3 = load_susceptibility_data("SusceptibilityCurve3.txt")

            if L > 1:
                ax0 = axs[n]
            else:
                ax0 = axs[n]

            plot_susceptibility_curve(T, X1, dX1, ax0, color='k', label=r'$B \parallel a$')
            plot_susceptibility_curve(T, X2, dX2, ax0, color='r', label=r'$B \parallel ab$')
            plot_susceptibility_curve(T, X3, dX3, ax0, color='b', label=r'$B \parallel c$')

            ax0.set_title(label + " = " + str(vals[n]))
            ax0.set_xlabel(r"$T$ (K)", fontsize=20)

        os.chdir('..')

    if L > 1:
        ax0 = axs[0]
    else:
        ax0 = axs[0]

    ax0.legend()
    ax0.set_ylabel(r'$\chi(T)$ ($\mu_B/T$)', fontsize=20)

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.show()

        




