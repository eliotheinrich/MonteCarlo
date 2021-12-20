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

def curie_fit(Ts, X1, X2, X3):
    def curie_weiss(T, C, Tc):
        return (T - Tc)/C

    N = len(X1)//2

    popt, _ = spo.curve_fit(curie_weiss, Ts[N:], 1/X1[N:])
    (C1, T1c) = popt
    popt, _ = spo.curve_fit(curie_weiss, Ts[N:], 1/X2[N:])
    (C2, T2c) = popt
    popt, _ = spo.curve_fit(curie_weiss, Ts[N:], 1/X3[N:])
    (C3, T3c) = popt

    print(C1, T1c)
    print(C2, T2c)
    print(C3, T3c)
    print(max(X3))
    print(max(X1))

    plt.plot(Ts, 1/X1, 'k-')
    plt.plot(Ts, (Ts - T1c)/C1, 'k--')
    plt.plot(Ts, 1/X2, 'r-')
    plt.plot(Ts, (Ts - T2c)/C2, 'r--')
    plt.plot(Ts, 1/X3, 'b-')
    plt.plot(Ts, (Ts - T3c)/C3, 'b--')
    plt.xlabel(r'$T$', fontsize=20)
    plt.ylabel(r'$\chi^{-1}$', fontsize=20)
    plt.show()




def plot_susceptibility_curve(Ts, X1, X2, X3, dX1, dX2, dX3, ax = None):
    if ax is None:
        ax = plt.gca()

    ax.errorbar(Ts, X1, yerr = dX1, color = 'k', linestyle = '-', label = r'$B \parallel ab$')
    ax.errorbar(Ts, X2, yerr = dX2, color = 'r', linestyle = '-', label = r'$B \parallel a$')
    ax.errorbar(Ts, X3, yerr = dX3, color = 'b', linestyle = '-', label = r'$B \parallel c$')
    ax.set_xlabel(r'$T$ (K)', fontsize=20)
    ax.set_ylabel(r'$\chi$', fontsize=20)
    ax.legend()

def get_run(X1, X2, X3, nrun):
    dX1 = np.std(X1[nrun], axis=1)
    dX2 = np.std(X2[nrun], axis=1)
    dX3 = np.std(X3[nrun], axis=1)

    avg_X1 = np.mean(X1[nrun], axis=1)
    avg_X2 = np.mean(X2[nrun], axis=1)
    avg_X3 = np.mean(X3[nrun], axis=1)

    return (avg_X1, avg_X2, avg_X3, dX1, dX2, dX3)

if __name__ == "__main__":
    os.chdir(sys.argv[1])

    filenames = [x for x in os.listdir() if x != "params.txt"]
    filenames1 = [x for x in filenames if x[19] == "1"]
    filenames2 = [x for x in filenames if x[19] == "2"]
    filenames3 = [x for x in filenames if x[19] == "3"]

    X1 = []
    X2 = []
    X3 = []
    for (f1, f2, f3) in zip(filenames1, filenames2, filenames3):
        Ts, X1t = load_susceptibility_data(f1)
        Ts, X2t = load_susceptibility_data(f2)
        Ts, X3t = load_susceptibility_data(f3)

        X1.append(X1t)
        X2.append(X2t)
        X3.append(X3t)

    X1 = np.array(X1)
    X2 = np.array(X2)
    X3 = np.array(X3)

    avg = input("Average all runs? y/n: ") == "y"
    if avg:
        dX1 = np.std(X1, axis=(0,2))
        dX2 = np.std(X2, axis=(0,2))
        dX3 = np.std(X3, axis=(0,2))

        avg_X1 = np.mean(X1, axis=(0,2))
        avg_X2 = np.mean(X2, axis=(0,2))
        avg_X3 = np.mean(X3, axis=(0,2))

        curie_fit(Ts, avg_X1, avg_X2, avg_X3)
        plot_susceptibility_curve(Ts, avg_X1, avg_X2, avg_X3, dX1, dX2, dX3)

    else:

        plot_susceptibility_curve(Ts, *get_run(X1, X2, X3, 0))
        plot_susceptibility_curve(Ts, *get_run(X1, X2, X3, 1))
        plot_susceptibility_curve(Ts, *get_run(X1, X2, X3, 2))

    plt.show()




