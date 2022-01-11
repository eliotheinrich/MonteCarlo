import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os

def load_stiffness_data(filename):
    i = filename.index('.')
    L = int(filename[9:i])

    with open(filename) as f:
        s = f.readline().split('\t')
        res = int(s[0])
    
    T = np.zeros(res)
    dE = np.zeros(res)
    err_dE = np.zeros(res)
    ddE = np.zeros(res)
    err_ddE = np.zeros(res)

    with open(filename) as f:
        line = f.readline()
        for i in range(res):
            line = f.readline()

            data = line.split('\t')
            T[i] = float(data[0])
            (dE[i], err_dE[i]) = (float(x) for x in data[1].split(','))
            (ddE[i], err_ddE[i]) = (float(x) for x in data[2].split(','))

    ρ = ddE - err_dE**2/T

    return L, T, ρ

def get_L0(L, T, ρ, Lmin = 0.5, Lmax = 3., num_Ls=100):
    ρL = np.zeros_like(ρ)
    intersect_points = np.zeros_like(L, dtype=np.float32)
    inds = np.zeros_like(L)
    res = len(ρ)

    T_KT = 0.
    min_L0 = 0.
    min_dev = np.inf
    for L0 in np.linspace(Lmin, Lmax, num_Ls):
        for n,Ln in enumerate(L):
            ρL[n] = ρ[n]/(1. + 1/(2.*np.log(Ln/L0)))

            ind = np.argwhere(np.diff(np.sign(ρL[n] - 2*T/np.pi))).flatten()[-1]

            inds[n] = ind
            x1 = T[ind]
            x2 = T[ind+1]
            y1 = ρL[n,ind]
            y2 = ρL[n,ind+1]

            intersect_points[n] = -np.pi*(x1*y2 - y1*x2)/(2*(x2 - x1) + np.pi*(y1 - y2))

        std_dev = np.std(intersect_points)
        if std_dev < min_dev:
            T_KT = intersect_points[0]
            min_L0 = L0
            min_dev = std_dev

    return min_L0, T_KT

def plot_stiffness_curve(L, T, ρ, L0, T_KT):
    ρL = np.array([ρ[n]/(1. + 1./(2.*np.log(L[n]/L0))) for n in range(len(ρ))])
    for n, L in enumerate(L):
        plt.plot(T, ρL[n], marker='o', label=f'L = {L}')

    plt.plot(T, 2./np.pi*T, 'k--', label=r'$2T/\pi$')

    plt.legend()
    plt.xlim(0., max(T))
    plt.ylim(-0.05, np.max(ρL) + 0.2)
    plt.xlabel(r'$T/J$', fontsize=15)
    plt.ylabel(r'$\Upsilon(L,T)/(1+(2\log(L/L_0))^{-1})$', fontsize=15)
#    plt.ylabel(r'$\Upsilon(L,T)$', fontsize=15)
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        os.chdir('data')
    else:
        os.chdir(sys.argv[1])

    fs = [f for f in os.listdir() if f[:9] == "stiffness"]
    _, T, ρ = load_stiffness_data(fs[0])

    N = len(fs)
    res = len(T)

    T = np.zeros(res)
    ρ = np.zeros((N, res))
    L = np.zeros(N, dtype=int)

    for n, f in enumerate(fs):
        if f[:9] == 'stiffness':
            i = f.index('.')
            L[n], T, ρ[n,:] = load_stiffness_data(f)

    inds = np.argsort(L)
    ρ = ρ[inds]
    L = L[inds]
    T = T


    L0, T_KT = get_L0(L, T, ρ, 0.1, 2.0, 10000)
#    L0 = 1/1.4
#    T_KT = 0.88
    print(L0, T_KT)
    plot_stiffness_curve(L, T, ρ, L0, T_KT)






