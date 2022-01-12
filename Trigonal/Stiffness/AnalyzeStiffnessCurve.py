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
        N = int(s[2])
        L = int(s[3])
    
    T = np.zeros(res)
    dEs = np.zeros((res, num_samples, 4))

    with open(filename) as f:
        line = f.readline()
        for i in range(res):
            line = f.readline()

            data = line.split('\t')[:-1]
            T[i] = float(data[0])
            dEs[i] = np.array([float(x) for x in data[1:]]).reshape((num_samples, 4))

    V = N*N*L
    d1E = dEs[:,:,0]
    d2E = dEs[:,:,1]
    d3E = dEs[:,:,2]
    d4E = dEs[:,:,3]

    U2 = (avg(d2E) - (avg(d1E**2) - avg(d1E)**2)/T)/V

    U4 = (6/T**3*avg(d1E)**4 + 12/T**2*avg(d1E)**2*(avg(d2E) - avg(d1E**2)/T) \
       + 3/T*avg(d2E)**2 - 6/T**3*avg(d1E**2)*avg(d2E) + 3/T**3*avg(d1E**2)**2 \
       + 4/T*avg(d3E)*avg(d1E) - 12/T**2*avg(d1E*d2E)*avg(d1E) + 4/T**3*avg(d1E**3)*avg(d1E) \
       + 1/T**3*avg(d1E**4) + 6/T**2*avg(d1E**2*d2E) - 3/T*avg(d2E**2) - 4/T*avg(d1E*d3E) \
       + avg(d4E))/V

    plt.plot(T, U2/max(np.abs(U2)))
    plt.plot(T, U4/max(np.abs(U4)))
    plt.show()

    plt.plot(T, np.std(d1E,axis=1), label='d1')
    #plt.plot(T, avg(d2E), label='d2')
    plt.plot(T, np.std(d3E,axis=1), label='d3')
    #plt.plot(T, avg(d4E), label='d4')
    plt.legend()
    plt.show()

    return L, T, U2, U4

def avg(A):
    return np.mean(A, axis=1)

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






