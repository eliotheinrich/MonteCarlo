import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
from matplotlib import cm
import sys
import os


kB = 0.0816

# Expected / measured
susc_factor = 7.88/11.74
#susc_factor = 99./11.74

def load_susceptibility_data(filename):
    with open(filename) as f:
        s = f.readline().split('\t')
        res = int(s[0])
    
    T = np.zeros(res)
    X = np.zeros(res)
    dX = np.zeros(res)
    E = np.zeros(res)
    dE = np.zeros(res)

    with open(filename) as f:
        line = f.readline()
        for i in range(res):
            line = f.readline()

            data = np.array([x.strip('\n') for x in line.split('\t')])
            T[i] = float(data[0])/kB
            (X[i], dX[i]) = (float(x) for x in data[1].split(','))
            (E[i], dE[i]) = (float(x) for x in data[2].split(','))

    return T, X*susc_factor, dX*susc_factor, E, dE


def curie_fit(T, X):
    def curie_weiss(T, C, Tc):
        return (T - Tc)/C

    N = len(X)//2

    popt, _ = spo.curve_fit(curie_weiss, T[N:], 1/X[N:])
    (C, Tc) = popt

    return C, Tc

def plot_susceptibility_curve(T, X, dX, ax, **kwargs):
    ax.errorbar(T, X, yerr = dX, **kwargs)

def plot_energy_curve(T, E, dE, ax, *args, **kwargs):
    ax.errorbar(T, E, yerr = dE, *args, **kwargs)

def plot_heat_capacity(T, dE, nmin, ax, *args, **kwargs):
    ax.plot(T[nmin:], ((dE/T)**2)[nmin:]/kB, *args, **kwargs)

if __name__ == "__main__":
    os.chdir(sys.argv[1])

    if len(os.listdir()) == 2:
        T, X, dX, E, dE = load_susceptibility_data("SusceptibilityCurve.txt")

        fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True)
        plot_susceptibility_curve(T, X, dX, axs[0], color='k', linestyle='-')
        axs[0].set_ylabel(r'$\chi$(T)', fontsize=16)
        plot_energy_curve(T, E, dE, axs[1], color='k', linestyle='-')
        axs[1].set_ylabel('u(T)/J', fontsize=16)
        plot_heat_capacity(T, dE, 0, axs[2], 'k-')
        axs[2].set_ylabel('c(T)/$k_B$', fontsize=16)
        axs[2].set_xlabel('T/J', fontsize=16)

        plt.subplots_adjust(wspace=0, hspace=0)
        plt.show()

    else:
        T, X1, dX1, E1, dE1 = load_susceptibility_data("SusceptibilityCurve1.txt")
        T, X2, dX2, E2, dE2 = load_susceptibility_data("SusceptibilityCurve2.txt")
        T, X3, dX3, E3, dE3 = load_susceptibility_data("SusceptibilityCurve3.txt")

        #C, _, _, _ = curie_fit(T, X1, X2, X3)

        kwargs1 = {'color':'k', 'linestyle':'-', 'label':r'$B \parallel ab$'}
        kwargs2 = {'color':'r', 'linestyle':'-', 'label':r'$B \parallel a$'}
        kwargs3 = {'color':'b', 'linestyle':'-', 'label':r'$B \parallel c$'}
        ax = plt.gca()
        plot_susceptibility_curve(T, X1, dX1, ax, **kwargs1)
        plot_susceptibility_curve(T, X2, dX2, ax, **kwargs2)
        plot_susceptibility_curve(T, X3, dX3, ax, **kwargs3)
        ax.set_xlabel(r'T (meV)')
        ax.set_ylabel(r'$\chi$')
        plt.show()




