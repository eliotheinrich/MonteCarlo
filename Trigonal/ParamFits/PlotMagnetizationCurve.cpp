import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
from matplotlib import cm
from AnalyzeMagnetizationCurve import *
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
        T, M, dM, E, dE = load_magnetization_data("MagnetizationCurve.txt")

        if L > 1:
            ax0 = axs[0,n]
            ax1 = axs[1,n]
            ax2 = axs[2,n]
        else:
            ax0 = axs[0]
            ax1 = axs[1]
            ax2 = axs[2]
        
        plot_magnetization_curve(T, M, dM, ax0, color='k', linestyle='-', linewidth=2.)
        plot_energy_curve(T, E, dE, ax1, color='k', linestyle='-', linewidth=2.)
        plot_heat_capacity(T, dE, 2, ax2, color='k', linestyle='-', linewidth=2.)

        ax0.set_title(label + " = " + str(vals[n]))
        ax2.set_xlabel(r"$T$ (K)", fontsize=20)
        os.chdir('..')

    if L > 1:
        ax0 = axs[0,0]
        ax1 = axs[1,0]
        ax2 = axs[2,0]
    else:
        ax0 = axs[0]
        ax1 = axs[1]
        ax2 = axs[2]

    ax0.set_ylabel(r'$M(T)$', fontsize=20)
    ax1.set_ylabel(r'$u(T)$ (meV)', fontsize=20)
    ax2.set_ylabel(r'$c(T)/k_B$', fontsize=20)

    plt.subplots_adjust(hspace=0, wspace=0)
    plt.show()

        




