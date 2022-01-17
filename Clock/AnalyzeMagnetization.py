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


fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
#axs[0].plot(Mxs, label=r'$M_x$')
#axs[0].plot(Mys, label=r'$M_y$')
axs[0].plot(Ms, 'k', linewidth=2., label=r'$M$')
axs[0].set_ylabel("Magnetization")
#axs[0].legend()
axs[1].plot(Es, 'r', linewidth=2.)
axs[1].set_ylabel("Energy")
axs[1].set_xticks([])
axs[1].set_xlabel("Number of Monte-Carlo steps")
plt.show()


