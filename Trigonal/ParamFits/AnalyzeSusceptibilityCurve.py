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
            X[i] = float(data[1])
            dX[i] = float(data[2])
            E[i] = float(data[3])
            dE[i] = float(data[4])


    return T, X*susc_factor, dX*susc_factor, E, dE


def curie_fit(Ts, X1, X2, X3, plot=False):
    def curie_weiss(T, C, Tc):
        return (T - Tc)/C

    N = len(X1)//2

    popt, _ = spo.curve_fit(curie_weiss, Ts[N:], 1/X1[N:])
    (C1, T1c) = popt
    popt, _ = spo.curve_fit(curie_weiss, Ts[N:], 1/X2[N:])
    (C2, T2c) = popt
    popt, _ = spo.curve_fit(curie_weiss, Ts[N:], 1/X3[N:])
    (C3, T3c) = popt

    C_avg = np.mean([C1, C2, C3])

    if plot:
        plt.plot(Ts, 1/X1, 'k-')
        plt.plot(Ts, (Ts - T1c)/C1, 'k--')
        plt.plot(Ts, 1/X2, 'r-')
        plt.plot(Ts, (Ts - T2c)/C2, 'r--')
        plt.plot(Ts, 1/X3, 'b-')
        plt.plot(Ts, (Ts - T3c)/C3, 'b--')
        plt.xlabel(r'$T$ (K)', fontsize=20)
        plt.ylabel(r'$\chi^{-1}$', fontsize=20)
        plt.show()

    return C_avg, T1c, T2c, T3c



def plot_susceptibility_curve(T, X1, X2, X3, dX1, dX2, dX3, ax = None):
    ax_none = ax is None
    if ax_none:
        ax = plt.gca()

    ax.errorbar(T, X1, yerr = dX1, color = 'k', linestyle = '-', label = r'$B \parallel ab$')
    ax.errorbar(T, X2, yerr = dX2, color = 'r', linestyle = '-', label = r'$B \parallel a$')
    ax.errorbar(T, X3, yerr = dX3, color = 'b', linestyle = '-', label = r'$B \parallel c$')
    if ax_none:
        ax.set_xlabel(r'$T$ (K)', fontsize=20)
        ax.set_ylabel(r'$\chi$', fontsize=20)
        ax.legend()

def plot_energy_curve(T, E1, E2, E3, dE1, dE2, dE3, ax = None):
    ax_none = ax is None
    if ax_none:
        ax = plt.gca()

    ax.errorbar(T, E1, yerr = dE1, color = 'k', linestyle = '-', label = r'$B \parallel ab$')
    ax.errorbar(T, E2, yerr = dE2, color = 'r', linestyle = '-', label = r'$B \parallel a$')
    ax.errorbar(T, E3, yerr = dE3, color = 'b', linestyle = '-', label = r'$B \parallel c$')
    if ax_none:
        ax.set_xlabel(r'$T$ (K)', fontsize=20)
        ax.set_ylabel(r'$U(T)$', fontsize=20)
        ax.legend()

def plot_heat_capacity(T, dE1, dE2, dE3, ax = None, nmin = 0):
    ax_none = ax is None
    if ax_none:
        ax = plt.gca()

    ax.plot(T[nmin:], ((dE1/T)**2)[nmin:]/kB, 'k-', label=r'$B \parallel ab$')
    ax.plot(T[nmin:], ((dE2/T)**2)[nmin:]/kB, 'r-', label=r'$B \parallel a$')
    ax.plot(T[nmin:], ((dE3/T)**2)[nmin:]/kB, 'b-', label=r'$B \parallel c$')
    if ax_none:
        ax.set_xlabel(r'$T$ (K)', fontsize=20)
        ax.set_ylabel(r'$E$', fontsize=20)
        ax.legend()

if __name__ == "__main__":
    os.chdir(sys.argv[1])

    T, X1, dX1, E1, dE1 = load_susceptibility_data("SusceptibilityCurve1.txt")
    T, X2, dX2, E2, dE2 = load_susceptibility_data("SusceptibilityCurve2.txt")
    T, X3, dX3, E3, dE3 = load_susceptibility_data("SusceptibilityCurve3.txt")

    #C, _, _, _ = curie_fit(T, X1, X2, X3)

    plot_susceptibility_curve(T, X1, X2, X3, dX1, dX2, dX3)
    plt.show()
    plt.show()




