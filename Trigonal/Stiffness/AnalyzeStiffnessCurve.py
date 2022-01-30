import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os

kB = 0.0816

def load_stiffness_data(filename):
    i = filename.index('.')

    with open(filename) as f:
        s = f.readline().split('\t')
        res = int(s[0])
        dtype_size = int(s[1])
        N = int(s[2])
        L = int(s[3])


    
    T = np.zeros(res)
    ps = np.zeros((res, dtype_size))
    dps = np.zeros((res, dtype_size))

    with open(filename) as f:
        line = f.readline()
        for i in range(res):
            line = f.readline()

            data = line.split('\t')
            T[i] = float(data[0])
            for j in range(dtype_size):
                (ps[i,j], dps[i,j]) = (float(x) for x in data[j+1].split(','))

    V = N*N*L


    (d4E, err_d4E) = (ps[:,0], dps[:,0])
    (d3E, err_d3E) = (ps[:,1], dps[:,1])
    (d1E, err_d1E) = (ps[:,2], dps[:,2])
    (d2E2, err_d2E2) = (ps[:,3], dps[:,3])
    (d2E, err_d2E) = (ps[:,4], dps[:,4])
    (d1E_d3E, err_d1E_d3E) = (ps[:,5], dps[:,5])
    (d1E2_d2E, err_d1E2_d2E) = (ps[:,6], dps[:,6])
    (d1E_d2E, err_d1E_d2E) = (ps[:,7], dps[:,7])
    (d1E2, err_d1E2) = (ps[:,8], dps[:,8])
    (d1E4, err_d1E4) = (ps[:,9], dps[:,9])
    (d1E3, err_d1E3) = (ps[:,10], dps[:,10])

    U2 = (d2E - (d1E2 - d1E**2)/T)/V
    dU2 = (err_d2E - (err_d1E2 - 2*err_d1E*d1E)/T)/V

    U4 = (d4E + 1/T*(4*d3E*d1E - 3*d2E2 + 3*d2E**2 - 4*d1E_d3E) \
             + 1/T**2*(6*d1E2_d2E - 12*d1E_d2E*d1E - 6*d1E2*d2E + 12*d1E**2*d2E) \
             + 1/T**3*(-d1E4 + 4*d1E3*d1E + 3*d1E2**2 - 12*d1E2*d1E**2 + 6*d1E**4))/V**2
    dU4 = (err_d4E + 1/T*(4*(err_d3E*d1E + err_d1E*d3E) - 3*2*err_d2E*d2E - 4*err_d1E_d3E) \
                + 1/T**2*(6*err_d1E2_d2E - 12*(err_d1E*d1E_d2E + err_d1E_d2E*d1E) \
                    - 6*(err_d1E2*d2E + err_d2E*d1E2) \
                    + 12*(2*err_d1E*d1E*d2E + err_d2E*d1E**2))
                + 1/T**3*(-err_d1E4 + 4*d1E3*d1E*(err_d1E3*d1E + err_d1E*d1E3) + 3*2*err_d1E2*d1E2 \
                         - 12*(err_d1E2*d1E**2 + 2*err_d1E*d1E*d1E2) + 6*4*err_d1E*d1E**3))/V**2
    plt.plot(T, err_d1E)
    plt.show()
    return N, T, U2, dU2, U4, dU4



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

def plot_stiffness_curve(L, T, U2, dU2, L0, T_KT):
    ax = plt.gca()
    U2L = np.array([U2[n]/(1. + 1./(2.*np.log(L[n]/L0))) for n in range(len(U2))])
    dU2L = np.array([dU2[n]/(1. + 1./(2.*np.log(L[n]/L0))) for n in range(len(U2))])

    for Ln, U2n, dU2n, U4n in zip(L, U2L, dU2L, U4):
        ax.errorbar(T/kB, U2n, dU2n, marker='o', label=f'L = {Ln}')
        #axs[1].plot(T, U4n, marker='o')

    ax.plot(T/kB, 2./np.pi*T, 'k--', label=r'$2T/\pi$')

    ax.legend()
    ax.set_xlim(0., max(T/kB))
    ax.set_ylim(-0.05, np.max(U2L) + 0.2)
    ax.set_xlabel(r'$T$ (K)', fontsize=15)
    ax.set_ylabel(r'$\Upsilon(L,T)/(1+(2\log(L/L_0))^{-1})$', fontsize=15)
    plt.subplots_adjust(hspace=0, wspace=0)
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        os.chdir('data')
    else:
        os.chdir(sys.argv[1])

    fs = [f for f in os.listdir() if f[:9] == "stiffness"]
    _, T, _, _, _, _= load_stiffness_data(fs[0])

    N = len(fs)
    res = len(T)

    T = np.zeros(res)
    U2 = np.zeros((N, res))
    dU2 = np.zeros((N, res))
    U4 = np.zeros((N, res))
    dU4 = np.zeros((N, res))
    L = np.zeros(N, dtype=int)

    for n, f in enumerate(fs):
        if f[:9] == 'stiffness':
            i = f.index('.')
            L[n], T, U2[n], dU2[n], U4[n], dU4[n] = load_stiffness_data(f)

    inds = np.argsort(L)
    U2 = U2[inds]
    dU2 = dU2[inds]
    U4 = U4[inds]
    dU4 = dU4[inds]
    L = L[inds]
    T = T


#    L0, T_KT = get_L0(L, T, U2, 0.01, 0.5, 10000)
    L0 = 1
    T_KT = 0.88
    print(L0, T_KT)
    plot_stiffness_curve(L, T, U2, dU2, L0, T_KT)






