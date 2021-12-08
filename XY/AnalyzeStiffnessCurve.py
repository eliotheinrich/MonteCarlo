import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def load_stiffness_data(filename):
    with open(filename) as f:
        n = 0
        while f.readline():
            n += 1
    
    T = np.zeros(n)
    ρ = np.zeros(n)
    with open(filename) as f:
        for i in range(n):
            line = f.readline()
            Ti, ρi = [float(x) for x in line.split('\t')]
            T[i] = Ti
            ρ[i] = ρi

    return T, ρ

def get_L0(Ls, Ts, ρs, Lmin = 0.5, Lmax = 3., num_Ls=100):
    ρsL = np.zeros_like(ρs)
    intersect_points = np.zeros_like(Ls,dtype=np.float32)
    inds = np.zeros_like(Ls)
    res = len(ρs)

    T_KT = 0.
    min_L0 = 0.
    min_dev = np.inf
    for L0 in np.linspace(Lmin, Lmax, num_Ls):
        for n,L in enumerate(Ls):
            ρsL[n] = ρs[n]/(1. + 1/(2.*np.log(L/L0)))

            ind = np.argwhere(np.diff(np.sign(ρsL[n] - 2*Ts/np.pi))).flatten()

            inds[n] = ind
            x1 = Ts[ind]
            x2 = Ts[ind+1]
            y1 = ρsL[n,ind]
            y2 = ρsL[n,ind+1]

            intersect_points[n] = -np.pi*(x1*y2 - y1*x2)/(2*(x2 - x1) + np.pi*(y1 - y2))

        std_dev = np.std(intersect_points)
        if std_dev < min_dev:
            T_KT = intersect_points[0]
            min_L0 = L0
            min_dev = std_dev

    return min_L0, T_KT


if __name__ == "__main__":
    os.chdir(sys.argv[1])

    fs = [f for f in os.listdir() if f[:9] == "stiffness"]
    T, _ = load_stiffness_data(fs[0])

    N = len(fs)
    res = len(T)

    Ts = np.zeros((N, res))
    ρs = np.zeros((N, res))
    Ls = np.zeros(N, dtype=int)

    for n, f in enumerate(fs):
        if f[:9] == 'stiffness':
            i = f.index('.')
            T, ρ = load_stiffness_data(f)
            Ts[n,:] = T
            ρs[n,:] = ρ
            Ls[n] = int(f[9:i])

    inds = np.argsort(Ls)
    Ts = Ts[inds]
    ρs = ρs[inds]
    Ls = Ls[inds]


#    L0, T_KT = get_L0(Ls, Ts[0,:], ρs, 0.1, 0.75, 10000)
    L0 = 1/1.4
    T_KT = 0.88
    print(L0, T_KT)

    for n, L in enumerate(Ls):
        plt.plot(Ts[n], ρs[n]/(1. + 1./(2.*np.log(Ls[n]/L0))), marker='o', label=f'L = {L}')

    plt.plot(Ts[0], 2./np.pi*Ts[0], 'k--', label=r'$2T/\pi$')
#    plt.axvline(T_KT, linestyle='-', color='k', alpha=0.5, label=r'$T_{KT}$')

    plt.legend()
    plt.xlim(0., max(Ts[0]))
    plt.ylim(-0.05, 1.1)
    plt.xlabel(r'$T$', fontsize=15)
    plt.ylabel(r'$\Upsilon(L)/(1+(2\log(L/L_0))^{-1})$', fontsize=15)
#    plt.title(r'$T_{KT} = $' + f'{T_KT:.3f}')
    plt.show()





