import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import sys
import os

α1 = np.array([1.0, 0.])
α2 = np.array([0.5, np.sqrt(3)/2.])

β1 = 4*np.pi/np.sqrt(3)*np.array([np.sqrt(3)/2, -1/2])
β2 = 4*np.pi/np.sqrt(3)*np.array([0,1])

def load_correlation_data(filename):
    i = filename.index('.')

    with open(filename) as f:
        s = f.readline().split('\t')
        res = int(s[0])
        N = int(s[1])
        L = int(s[3])
    
    T = np.zeros(res)
    C = np.zeros((res, N, N, L))
    dC = np.zeros((res, N, N, L))


    with open(filename) as f:
        f.readline()

        for i in range(res):
            line = f.readline()

            data = line.split('\t')
            T[i] = float(data[0])
            nums = np.array([[float(x) for x in y.split(',')] for y in data[1:]]).T
            C[i] = nums[0,:].reshape((N,N,L), order='F')
            dC[i] = nums[1,:].reshape((N,N,L), order='F')

    return T, recenter(C, axes=(1,2,3)), recenter(dC, axes=(1,2,3))

def affine_meshgrid(a1, a2, N1, N2, r0=None):
    X = np.zeros((N1, N2))
    Y = np.zeros((N1, N2))

    if r0 is None:
        r0 = np.array([0,0])

    for i in range(N1):
        for j in range(N2):
            r = a1*i + a2*j
            X[i,j] = r[0] - r0[0]
            Y[i,j] = r[1] - r0[1]

    return X, Y

def correlation_length(C):
    def exp_func(r, ξ):
        return np.exp(-r/ξ)

    N = C.shape[1]

    X, Y = affine_meshgrid(α1, α2, N, N, N*α1/2 + N*α2/2)

    r = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            r[i,j] = np.sqrt(X[i,j]**2 + Y[i,j]**2)

    C = C.reshape(N*N)
    r = r.reshape(N*N)

    #(ξ,), _ = spo.curve_fit(exp_func, r, C/C[0])
    ξ = -r[-1]/np.log(C[-1])
    #print(f'ξ = {ξ}') 
    #plt.plot(r, [exp_func(x, 2.) for x in r], label='Fit')

    #plt.plot(r, C, label='Data')
    #plt.legend()
    #plt.show()
    return ξ



def fourier_transform(C):
    N = C.shape[1]

    KX, KY = affine_meshgrid(β1/N, β2/N, N, N, β1/2 + β2/2)

    Ck = np.fft.fft2(C) 

    return KX, KY, Ck

def recenter(F, axes = (0)):
    axes = np.array(axes)
    lens = np.array(F.shape)[axes]

    for n,L in enumerate(lens):
        F = np.roll(F, L//2, axis=axes[n])

    return F

def plot_structure_factor(C, ax):
    N = C.shape[1]
    L = C.shape[-1]

    KX, KY, Ck = fourier_transform(C)

    Ck = recenter(Ck, axes=(0,1))

    ax.pcolor(KX, KY, np.abs(Ck), shading='auto')
    ax.set_aspect('equal')
    ax.axis('off')


def plot_correlation_function(C, ax):
    N = C.shape[1]
    L = C.shape[-1]

    X, Y = affine_meshgrid(α1, α2, N, N, N*α1/2 + N*α2/2)

    ax.pcolor(X, Y, C, vmin = -1, vmax = 1, shading='auto')
    ax.set_aspect('equal')
    ax.axis('off')

if __name__ == "__main__":
    os.chdir(sys.argv[1])
    T, C, dC = load_correlation_data("correlation.txt")
    N = C.shape[1]
    L = C.shape[-1]

    Ci = C[20,:,:,L//2]
    
    ξ = [correlation_length(x) for x in C[1:]]
    plt.plot(T[1:], ξ)
    plt.show()


#    ax = plt.gca()
#    plot_correlation_function(Ci, ax)
#    plt.show()

#    ax = plt.gca()
#    plot_structure_factor(Ci, ax)
#    plt.show()


