import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
from matplotlib import cm
from AnalyzeSusceptibilityCurve import *
import sys
import os

if __name__ == "__main__":
    os.chdir('data')

    folders = [x for x in os.listdir() if x[:len(sys.argv[1])] == sys.argv[1]]
    folders = sorted(folders, key = lambda x: float(x[(x.index('_')+1):]))
    print(folders)

    L = len(folders)

    fig, axs = plt.subplots(nrows=2, ncols=L, sharex=True)
    print(len(axs))
    print(len(axs[0]))
    
    for n,folder in enumerate(folders):
        os.chdir(folder)
        T, X1, dX1, E1, dE1 = load_susceptibility_data("SusceptibilityCurve1.txt")
        T, X2, dX2, E2, dE2 = load_susceptibility_data("SusceptibilityCurve2.txt")
        T, X3, dX3, E3, dE3 = load_susceptibility_data("SusceptibilityCurve3.txt")

        plot_susceptibility_curve(T, X1, X2, X3, dX1, dX2, dX3, axs[0,n])
        plot_heat_capacity(T, dE1, dE2, dE3, axs[1,n])
        #plot_energy_curve(T, E1, E2, E3, dE1, dE2, dE3, axs[1,n])

        axs[1,n].set_xlabel(r"$T/J$", fontsize=20)
        os.chdir('..')

    axs[0,-1].legend()
    axs[0,0].set_ylabel(r'$\chi(T)$', fontsize=20)
    axs[1,0].set_ylabel(r'$c(T)$', fontsize=20)

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.show()

        




