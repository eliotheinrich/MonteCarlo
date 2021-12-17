import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os

def load_magnetization_data(filename):
    i = filename.index('.')
    L = int(filename[13:i])

    with open(filename) as f:
        s = f.readline().split('\t')
        res = int(s[0])
        num_samples = int(s[1])
    
    Ts = np.zeros(res)
    Ms = np.zeros((res, num_samples))

    with open(filename) as f:
        line = f.readline()
        for i in range(res):
            line = f.readline()

            data = np.array([x.strip('\n') for x in line.split('\t')])
            Ti = float(data[0])
            data = data[1:]

            for n,x in enumerate(data):
                Ms[i,n] = float(x.strip()[1:-1])

            Ts[i] = Ti

    return L, Ts, Ms


def plot_magnetization_curve(Ts, Ls, Ms, T_KT = None):
    for n, L in enumerate(Ls):
        plt.plot(Ts, Ms[n], marker='o', label=f'L = {L}')

    if T_KT is not None:
        plt.axvline(T_KT, linestyle='-', color='k', alpha=0.5, label=r'$T_{KT}$')

    plt.legend()
    plt.xlim(0., max(Ts))
    plt.ylim(-0.01, 1.1)
    plt.xlabel(r'$T/J$', fontsize=15)
    plt.ylabel(r'$M$', fontsize=15)
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        os.chdir('data')
    else:
        os.chdir(sys.argv[1])

    fs = [f for f in os.listdir() if f[:13] == "magnetization"]
    _, _, Ms = load_magnetization_data(fs[0])

    N = len(fs)
    res, num_samples = Ms.shape

    Ts = np.zeros(res)
    Ms = np.zeros((N, res))
    Ls = np.zeros(N, dtype=int)

    for n, f in enumerate(fs):
        if f[:13] == 'magnetization':
            i = f.index('.')
            L, Ts, M = load_magnetization_data(f)
            Ls[n] = L
            Ms[n,:] = np.mean(M, axis=1)

    inds = np.argsort(Ls)
    Ls = Ls[inds]
    Ms = Ms[inds]


    plot_magnetization_curve(Ts, Ls, Ms, 1.38)






