import sys
import numpy as np
import matplotlib.pyplot as plt

Ts = []
Es = []
Cs = []

filename = sys.argv[1]
with open(filename) as f:
    for line in f.readlines():
        T, C, E = line.split('\t')
        Ts.append(float(T))
        Cs.append(float(C))
        Es.append(float(E))

Cs = np.array(Cs)
Es = np.array(Es)


fig, axs = plt.subplots(nrows=2,ncols=1, sharex=True)
axs[0].plot(Ts, Cs)
axs[0].set_ylabel("Connectivity")
axs[1].plot(Ts, Es)
axs[1].set_ylabel("Energy")
plt.show()


