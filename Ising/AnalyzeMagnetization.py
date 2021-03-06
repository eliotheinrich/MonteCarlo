import sys
import numpy as np
import matplotlib.pyplot as plt

Es = []
Ms = []

filename = sys.argv[1]
with open(filename) as f:
    for line in f.readlines():
        M, E = line.split('\t')
        Ms.append(float(M))
        Es.append(float(E))

Ms = np.array(Ms)
Es = np.array(Es)


fig, axs = plt.subplots(nrows=2,ncols=1, sharex=True)
axs[0].plot(Ms)
axs[0].set_ylabel("Magnetization")
axs[1].plot(Es)
axs[1].set_ylabel("Energy")
plt.show()


