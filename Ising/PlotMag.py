import numpy as np
import matplotlib.pyplot as plt
import sys

def load_mag_data(filename):
    with open(filename) as f:
        n = 0
        while f.readline():
            n += 1
    
    T = np.zeros(n)
    M = np.zeros(n)
    with open(filename) as f:
        for i in range(n):
            line = f.readline()
            Ti, Mi = [float(x) for x in line.split('\t')]
            T[i] = Ti
            M[i] = np.abs(Mi)

    return T, M

T, M = load_mag_data(sys.argv[1])

plt.plot(T, M)

plt.xlabel(r'$T$', fontsize=15)
plt.ylabel(r'$M$', fontsize=15)
plt.axvline(2./(np.log(1.+np.sqrt(2.))), linestyle='--', color='k', alpha=0.5)
plt.show()












