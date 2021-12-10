import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def load_vorticity_data(filename):
    i = filename.index('.')
    L = int(filename[9:i])

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

    return L, Ts, v1s, v2s

def plot_vorticity_curve(Ts, Ls, v1s, v2s):
    for n, L in enumerate(Ls):
        plt.plot(Ts, v1s[n,:], marker='o', label=f'L = {L}, ' + r'$\nu = +$')
        plt.plot(Ts, v2s[n,:], marker='o', label=f'L = {L}, ' + r'$\nu = -$')

    plt.legend()
    plt.xlabel(r'$T$', fontsize=15)
    plt.ylabel(r'$\nu_{\pm}$', fontsize=15)
    plt.show()



if __name__ == "__main__":
    os.chdir('data')

    fs = [f for f in os.listdir() if f[:9] == "vorticity"]
    _, _, v1s, _ = load_vorticity_data(fs[0])

    N = len(fs)
    res, num_samples = v1s.shape

    Ts = np.zeros(res)
    v1s = np.zeros((N, res, num_samples))
    v2s = np.zeros((N, res, num_samples))
    Ls = np.zeros(N, dtype=int)

    for n, f in enumerate(fs):
        if f[:9] == 'vorticity':
            i = f.index('.')
            L, Ts, v1, v2 = load_vorticity_data(f)
            Ls[n] = L
            v1s[n,:,:] = v1
            v2s[n,:,:] = v2

    print(Ts)
    inds = np.argsort(Ls)
    v1s = v1s[inds]
    v2s = v2s[inds]
    Ls = Ls[inds]

    avg_v1s = np.mean(v1s, axis=2)
    avg_v2s = np.mean(v2s, axis=2)

    print(avg_v1s)


    plot_vorticity_curve(Ts, Ls, avg_v1s, avg_v2s)






