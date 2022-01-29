import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os

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

    U2 = (ps[:,4] - (ps[:,8] - ps[:,2]**2)/T)/V
    dU2 = (dps[:,4] - 2*dps[:,2]*ps[:,2]/T)/V

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

    U4 = (d4E + 1/T*(4*d3E*d1E - 3*d2E2 + 3*d2E**2 - 4*d1E_d3E) \
             + 1/T**2*(6*d1E2_d2E - 12*d1E_d2E*d1E - 6*d1E2*d2E + 12*d1E**2*d2E) \
             + 1/T**3*(-d1E4 + 4*d1E3*d1E + 3*d1E2**2 - 12*d1E2*d1E**2 + 6*d1E**4))/V**2
    dU4 = (err_d4E + 1/T*(4*d3E*d1E*(err_d3E/np.abs(d3E) + err_d1E/np.abs(d1E)) - 3*d2E**2*(2*err_d2E/np.abs(d2E)) - 4*err_d1E_d3E) \
                + 1/T**2*(6*err_d1E2_d2E - 12*d1E_d2E*d1E*(err_d1E/np.abs(d1E) + err_d1E_d2E/np.abs(d1E_d2E)) \
                    - 6*d1E2*d2E*(err_d1E2/np.abs(d1E2) + err_d2E/np.abs(d2E)) \
                    + 12*d1E**2*d2E*(2*err_d1E/np.abs(d1E) + err_d2E/np.abs(d2E)))
                + 1/T**3*(-err_d1E4 + 4*d1E3*d1E*(err_d1E3/np.abs(d1E3) + err_d1E/np.abs(d1E)) + 3*d1E2**2*(2*err_d1E2/np.abs(d1E2)) \
                         - 12*d1E2*d1E**2*(err_d1E2/np.abs(d1E2) + 2*err_d1E/np.abs(d1E)) + 6*d1E**4*(4*err_d1E/np.abs(d1E))))/V**2
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

def plot_raw_data(L, T, U2, dU2, U4, dU4):
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
    for n,Ln in enumerate(L):
        axs[0].errorbar(T, U2[n], dU2[n], marker='o', label=f'L = {Ln}')
        axs[1].errorbar(T, U4[n], dU4[n], marker='o', label=f'L = {Ln}')

    axs[0].legend()
    axs[0].set_ylabel(r'$\Upsilon_2$', fontsize=15)
    #axs[0].set_ylim(0.0, 1.2)
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
    _, T, _, _, _, _ = load_stiffness_data(fs[0])

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
            L[n], T, U2[n,:], dU2[n,:], U4[n,:], dU4[n,:] = load_stiffness_data(f)

    inds = np.argsort(L)
    U2 = U2[inds]
    dU2 = dU2[inds]
    U4 = U4[inds]
    dU4 = dU4[inds]
    L = L[inds]

    plot_raw_data(L, T, U2, dU2, U4, dU4)


    L0, T_KT = get_L0(L, T, U2, 0.1, 2., 10000)
#    L0 = 1/1.4
#    T_KT = 0.88
    print(L0, T_KT)
#    plot_stiffness_curve(L, T, U2, L0, T_KT)






