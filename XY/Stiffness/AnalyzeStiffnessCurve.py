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
        num_samples = int(s[1])
    
    Ts = np.zeros(res)
    rhos = np.zeros(res)
    dEs = np.zeros((res, num_samples))
    ddEs = np.zeros((res, num_samples))

    with open(filename) as f:
        line = f.readline()
        for i in range(res):
            line = f.readline()

            data = np.array([x.strip('\n') for x in line.split('\t')])
            Ti = float(data[0])
            data = data[1:]

            for n,x in enumerate(data):
                dE, ddE = x.split(',')
                dEs[i, n] = float(dE[1:])
                ddEs[i, n] = float(ddE[:-1])



            Ts[i] = Ti

    return L, Ts, dEs, ddEs

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

            ind = np.argwhere(np.diff(np.sign(ρsL[n] - 2*Ts/np.pi))).flatten()[-1]

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

def plot_stiffness_curve(Ts, ρs, L0, T_KT):
    ρsL = np.array([ρs[n]/(1. + 1./(2.*np.log(Ls[n]/L0))) for n in range(len(ρs))])
    for n, L in enumerate(Ls):
        plt.plot(Ts, ρsL[n], marker='o', label=f'L = {L}')

    plt.plot(Ts, 2./np.pi*Ts, 'k--', label=r'$2T/\pi$')
#    plt.axvline(T_KT, linestyle='-', color='k', alpha=0.5, label=r'$T_{KT}$')

    plt.legend()
    plt.xlim(0., max(Ts))
    plt.ylim(-0.05, np.max(ρsL) + 0.2)
    plt.xlabel(r'$T$', fontsize=15)
    plt.ylabel(r'$\Upsilon(L,T)/(1+(2\log(L/L_0))^{-1})$', fontsize=15)
#    plt.ylabel(r'$\Upsilon(L,T)$', fontsize=15)
    plt.show()

def plot_dEs(Ts, dEs, ddEs, l=0):
    fig, axs = plt.subplots(nrows = 1, ncols = 2)
    colors = cm.seismic
    res = len(dEs[0])
    for n in range(res):
        axs[0].plot(ddEs[l,n,:], color = colors(n/res))
        axs[1].plot(dEs[l,n,:], color = colors(n/res))
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        os.chdir('data')
    else:
        os.chdir(sys.argv[1])

    fs = [f for f in os.listdir() if f[:9] == "stiffness"]
    _, _, dEs, _ = load_stiffness_data(fs[0])

    N = len(fs)
    res, num_samples = dEs.shape

    Ts = np.zeros(res)
    ρs = np.zeros((N, res))
    dEs = np.zeros((N, res, num_samples))
    ddEs = np.zeros((N, res, num_samples))
    Ls = np.zeros(N, dtype=int)

    for n, f in enumerate(fs):
        if f[:9] == 'stiffness':
            i = f.index('.')
            L, Ts, dE, ddE = load_stiffness_data(f)
            Ls[n] = L
            dEs[n,:,:] = dE
            ddEs[n,:,:] = ddE
            ρs[n,:] = (np.mean(ddE, axis=1) - np.std(dE, axis=1)**2/Ts)/(L**2)

    inds = np.argsort(Ls)
    dEs = dEs[inds]
    ddEs = ddEs[inds]
    ρs = ρs[inds]
    Ls = Ls[inds]


    L0, T_KT = get_L0(Ls, Ts, ρs, 0.0001, 0.001, 10000)
#    L0 = 1/1.4
#    T_KT = 0.88
    print(L0, T_KT)
    plot_stiffness_curve(Ts, ρs, L0, T_KT)






