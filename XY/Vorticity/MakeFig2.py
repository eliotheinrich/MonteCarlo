import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def load_vorticity_data(filename):
    i = filename.index('.')

    with open(filename) as f:
        s = f.readline().split('\t')
        res = int(s[0])
        num_samples = int(s[1])
    
    Ts = np.zeros(res)
    v1s = np.zeros((res, num_samples))
    v2s = np.zeros((res, num_samples))

    with open(filename) as f:
        f.readline()
        for i in range(res):
            line = f.readline()
            data = np.array(line.split('\t'))
            Ti = float(data[0])
            data = data[1:-1]

            for n,x in enumerate(data):
                v1, v2 = x.split(',')
                v1s[i, n] = float(v1[1:])
                v2s[i, n] = float(v2[:-1])



            Ts[i] = Ti

    return Ts, v1s, v2s

def get_mu(Ts, v1s, v2s, minn=None, maxn=None):
    vs = np.abs(v1s) + np.abs(v2s)
    betas = 1/Ts

    inds = np.argsort(betas)
    betas = betas[inds]

    logv = np.log(vs)
    logv = logv[inds]

    if minn is None:
        minn = 0
    if maxn is None:
        maxn = len(Ts)

    p = np.polyfit(betas[minn:maxn], logv[minn:maxn], 1) 
    mu = p[0]

    plt.plot(betas, logv)
    plt.plot(betas, p[0]*betas + p[1], 'r-')
    plt.xlabel(r'$(T/J)^{-1}$', fontsize=15)
    plt.ylabel(r'$\log \rho$', fontsize=15)
    plt.show()

    print(mu)

def plot_vorticity_curve(Ts, v1s, v2s, label=None, T_KT=0.88, ax=None):
    if ax is None:
        ax = plt.gca()
    if label is None:
        line = ax.plot(Ts, np.abs(v1s) + np.abs(v2s), marker='o')
    else:
        line = ax.plot(Ts, np.abs(v1s) + np.abs(v2s), marker='o', label=label)

    if T_KT is not None:
        ax.axvline(T_KT, color=line[0].get_color(), linestyle='--', alpha=0.5)
    if label is not None:
        ax.legend()
    ax.set_xlabel(r'$T/J$', fontsize=15)
    ax.set_ylabel(r'$\rho(T)$', fontsize=15)



if __name__ == "__main__":
    Ts, v1s, v2s = load_vorticity_data('data/square/vorticity64.txt')
    Tt, v1t, v2t = load_vorticity_data('data/trigonal/vorticity64.txt')

    res, num_samples = v1s.shape

    avg_v1s = np.mean(v1s, axis=1)
    avg_v2s = np.mean(v2s, axis=1)
    avg_v1t = np.mean(v1t, axis=1)
    avg_v2t = np.mean(v2t, axis=1)
    avg_v1d = np.mean(v1d, axis=1)
    avg_v2d = np.mean(v2d, axis=1)

    ax = plt.gca()
    plot_vorticity_curve(Ts, avg_v1s, avg_v2s, "L = 64, square lattice", 0.883, ax)
    plot_vorticity_curve(Ts, avg_v1t, avg_v2t, "L = 64, hexagonal lattice", 1.38, ax)
    plt.show()

    #get_mu(Ts, avg_v1s, avg_v2s, 43)
    #plt.show()





