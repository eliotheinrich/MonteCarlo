import sys
import numpy as np
import matplotlib.pyplot as plt

Es = []
Mxs = []
Mys = []
Mzs = []

filename = sys.argv[1]
with open(filename) as f:
    for line in f.readlines():
        Mx, My, Mz, E = line.split('\t')
        Mxs.append(float(Mx))
        Mys.append(float(My))
        Mzs.append(float(Mz))
        Es.append(float(E))

Mxs = np.array(Mxs)
Mys = np.array(Mys)
Mzs = np.array(Mzs)
Ms = np.array([np.linalg.norm(M) for M in zip(Mxs, Mys, Mzs)])
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


