import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os

def avg(A):
    return np.mean(A, axis=1)

def err(A):
    return np.std(A, axis=1)

def load_stiffness_data(filename):
    i = filename.index('.')

    with open(filename) as f:
        s = f.readline().split('\t')
        res = int(s[0])
        num_samples = int(s[1])
        N = int(s[2])
        L = int(s[3])
    
    T = np.zeros(res)
    dEs = np.zeros((res, num_samples, 7))

    with open(filename) as f:
        line = f.readline()
        for i in range(res):
            line = f.readline()

            data = line.split('\t')[:-1]
            T[i] = float(data[0])
            dEs[i] = np.array([float(x) for x in data[1:]]).reshape((num_samples, 7))

    V = N*N*L
    d1E = dEs[:,:,0]
    d2E = dEs[:,:,1]
    d3E = dEs[:,:,2]
    d4E = dEs[:,:,3]
    U2s = dEs[:,:,4]
    e = dEs[:,:,5]
    s4 = dEs[:,:,6]

    U2 = (avg(d2E) - (avg(d1E**2) - avg(d1E)**2)/T)/V

    U4 = (6/T**3*avg(d1E)**4 + 12/T**2*avg(d1E)**2*(avg(d2E) - avg(d1E**2)/T) \
       + 3/T*avg(d2E)**2 - 6/T**2*avg(d1E**2)*avg(d2E) + 3/T**3*avg(d1E**2)**2 \
       + 4/T*avg(d3E)*avg(d1E) - 12/T**2*avg(d1E*d2E)*avg(d1E) + 4/T**3*avg(d1E**3)*avg(d1E) \
       - 1/T**3*avg(d1E**4) + 6/T**2*avg(d1E**2*d2E) - 3/T*avg(d2E**2) - 4/T*avg(d1E*d3E) \
       + avg(d4E))/V/V

    U22 = avg(U2s)
    U42 = (-4*avg(U2s) + 3*(avg(e) - V/T*err(U2s)**2) + 2*(V/T)**3*avg(s4))/V


    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
    axs[0].plot(T, U2, 'r-', label='General equation')
    axs[0].plot(T, U22, 'b-', label='Square XY equation')
    axs[0].set_ylabel(r'$\Upsilon_2$', fontsize=15)
    axs[0].set_title(f'L = {N}', fontsize=15)
    axs[0].legend()

    axs[1].plot(T, U4, 'r-')
    axs[1].plot(T, U42, 'b-')
    axs[1].set_ylabel(r'$\Upsilon_4$', fontsize=15)
    axs[1].set_xlabel('T', fontsize=15)
    plt.subplots_adjust(hspace=0, wspace=0)
    plt.show()

    return N, T, U2, U4

def get_L0(L, T, U2, Lmin = 0.5, Lmax = 3., num_Ls=100):
    U2L = np.zeros_like(U2)
    intersect_points = np.zeros_like(L, dtype=np.float32)
    inds = np.zeros_like(L)
    res = len(U2)

    T_KT = 0.
    min_L0 = 0.
    min_dev = np.inf
    for L0 in np.linspace(Lmin, Lmax, num_Ls):
        for n,Ln in enumerate(L):
            U2L[n] = U2[n]/(1. + 1/(2.*np.log(Ln/L0)))

            ind = np.argwhere(np.diff(np.sign(U2L[n] - 2*T/np.pi))).flatten()[-1]

            inds[n] = ind
            x1 = T[ind]
            x2 = T[ind+1]
            y1 = U2L[n,ind]
            y2 = U2L[n,ind+1]

            intersect_points[n] = -np.pi*(x1*y2 - y1*x2)/(2*(x2 - x1) + np.pi*(y1 - y2))

        std_dev = np.std(intersect_points)
        if std_dev < min_dev:
            T_KT = intersect_points[0]
            min_L0 = L0
            min_dev = std_dev

    return min_L0, T_KT

def plot_stiffness_curve(L, T, U2, L0, T_KT):
    U2L = np.array([U2[n]/(1. + 1./(2.*np.log(L[n]/L0))) for n in range(len(U2))])
    for n, L in enumerate(L):
        plt.plot(T, U2L[n], marker='o', label=f'L = {L}')

    plt.plot(T, 2./np.pi*T, 'k--', label=r'$2T/\pi$')

    plt.legend()
    plt.xlim(0., max(T))
    plt.ylim(-0.05, np.max(U2L) + 0.2)
    plt.xlabel(r'$T/J$', fontsize=15)
    plt.ylabel(r'$\Upsilon(L,T)/(1+(2\log(L/L_0))^{-1})$', fontsize=15)
    plt.show()

def plot_raw_data(L, T, U2, U4):
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
    for n,Ln in enumerate(L):
        axs[0].plot(T, U2[n], marker='o', label=f'L = {Ln}')
        axs[1].plot(T, U4[n], marker='o', label=f'L = {Ln}')

    axs[0].legend()
    axs[0].set_ylabel(r'$\Upsilon_2$', fontsize=15)
    axs[1].set_ylabel(r'$\Upsilon_4$', fontsize=15)
    axs[1].set_xlabel(r'$T$', fontsize=15)
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        os.chdir('data')
    else:
        os.chdir(sys.argv[1])

    fs = [f for f in os.listdir() if f[:9] == "stiffness"]
    _, T, _, _ = load_stiffness_data(fs[0])

    N = len(fs)
    res = len(T)

    T = np.zeros(res)
    U2 = np.zeros((N, res))
    U4 = np.zeros((N, res))
    L = np.zeros(N, dtype=int)

    for n, f in enumerate(fs):
        if f[:9] == 'stiffness':
            i = f.index('.')
            L[n], T, U2[n,:], U4[n,:] = load_stiffness_data(f)

    inds = np.argsort(L)
    U2 = U2[inds]
    U4 = U4[inds]
    L = L[inds]

    plot_raw_data(L, T, U2, U4)


#    L0, T_KT = get_L0(L, T, U2, 0.1, 2.0, 10000)
#    L0 = 1/1.4
#    T_KT = 0.88
#    print(L0, T_KT)
#    plot_stiffness_curve(L, T, U2, L0, T_KT)






