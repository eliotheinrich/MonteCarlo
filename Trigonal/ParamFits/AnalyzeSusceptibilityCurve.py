import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
from matplotlib import cm
import sys
import os

def load_susceptibility_data(filename):
    with open(filename) as f:
        s = f.readline().split('\t')
        res = int(s[0])
        num_samples = int(s[1])
    
    Ts = np.zeros(res)
    Xs = np.zeros((res, num_samples))

    with open(filename) as f:
        line = f.readline()
        for i in range(res):
            line = f.readline()

            data = np.array([x.strip('\n') for x in line.split('\t')])
            Ti = float(data[0])
            data = data[1:]

            for n,x in enumerate(data):
                Xs[i,n] = x[1:-1]

            Ts[i] = Ti

    return Ts, Xs

def curie_fit(Ts, M1s, M2s, M3s):
    def curie_weiss(T, C, Tc):
        return (T - Tc)/C

    avg_X1 = np.mean(M1s, axis=1)
    avg_X2 = np.mean(M2s, axis=1)
    avg_X3 = np.mean(M3s, axis=1)
    N = len(avg_X1)//2

    popt, _ = spo.curve_fit(curie_weiss, Ts[N:], 1/avg_X1[N:])
    (C1, T1c) = popt
    popt, _ = spo.curve_fit(curie_weiss, Ts[N:], 1/avg_X2[N:])
    (C2, T2c) = popt
    popt, _ = spo.curve_fit(curie_weiss, Ts[N:], 1/avg_X3[N:])
    (C3, T3c) = popt

    print(C1, T1c)
    print(C2, T2c)
    print(C3, T3c)
    print(max(avg_X3))
    print(max(avg_X1))

    plt.plot(Ts, 1/avg_X1, 'k-')
    plt.plot(Ts, (Ts - T1c)/C1, 'k--')
    plt.plot(Ts, 1/avg_X2, 'r-')
    plt.plot(Ts, (Ts - T2c)/C2, 'r--')
    plt.plot(Ts, 1/avg_X3, 'b-')
    plt.plot(Ts, (Ts - T3c)/C3, 'b--')
    plt.xlabel(r'$T$', fontsize=20)
    plt.ylabel(r'$\chi^{-1}$', fontsize=20)
    plt.show()




def plot_susceptibility_curve(Ts, M1s, M2s, M3s, ax = None):
    if ax is None:
        ax = plt.gca()
    
    avg_X1 = np.mean(X1s, axis=1)
    avg_X2 = np.mean(X2s, axis=1)
    avg_X3 = np.mean(X3s, axis=1)

    err_X1 = np.std(X1s, axis=1)
    err_X2 = np.std(X2s, axis=1)
    err_X3 = np.std(X3s, axis=1)

    ax.errorbar(Ts, avg_X1, yerr = err_X1, color = 'k', linestyle = '-', label = r'$B \parallel ab$')
    ax.errorbar(Ts, avg_X2, yerr = err_X2, color = 'r', linestyle = '-', label = r'$B \parallel a$')
    ax.errorbar(Ts, avg_X3, yerr = err_X3, color = 'b', linestyle = '-', label = r'$B \parallel c$')
    ax.set_xlabel(r'$T$ (K)', fontsize=20)
    ax.set_ylabel(r'$\chi$', fontsize=20)
    ax.legend()


if __name__ == "__main__":
    os.chdir(sys.argv[1])

    _, X1s = load_susceptibility_data("SusceptibilityCurve1.txt")
    _, X2s = load_susceptibility_data("SusceptibilityCurve2.txt")
    Ts, X3s = load_susceptibility_data("SusceptibilityCurve3.txt")

    curie_fit(Ts, X1s, X2s, X3s)

    plot_susceptibility_curve(Ts, X1s, X2s, X3s)
    plt.show()




