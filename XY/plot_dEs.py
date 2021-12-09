import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmaps

def load_dEs(filename):
    dEs = []
    with open(filename) as f:
        for line in f.readlines():
            dEs.append(np.array(line.split('\t'), dtype=np.float32))

    return np.array(dEs)


dEs = load_dEs('dEs.txt')
ddEs = load_dEs('ddEs.txt')

avg_dEs = np.mean(dEs, axis=1)
std_dEs = np.std(dEs, axis=1)
avg_ddEs = np.mean(ddEs, axis=1)
std_ddEs = np.std(dEs, axis=1)

Ts = np.linspace(0.1, 3., len(std_dEs))

fig, axs = plt.subplots(nrows=1, ncols=2)
axs[0].plot(std_dEs**2)
axs[0].set_title(r'$\left\langle \frac{dE}{d\alpha} \right\rangle$', fontsize=20)
axs[1].plot(avg_ddEs - std_dEs**2/Ts)
axs[1].set_title(r'$\left\langle \frac{d^2 E}{d\alpha^2} \right\rangle$', fontsize=20)

plt.show()


fig, axs = plt.subplots(nrows=1, ncols=2)
axs[0].set_title(r'$\left\langle \frac{dE}{d\alpha} \right\rangle$', fontsize=20)
axs[1].set_title(r'$\left\langle \frac{d^2 E}{d\alpha^2} \right\rangle$', fontsize=20)

colors = cmaps.seismic
num_terms = 20
for i in range(num_terms):
    axs[0].plot(dEs[i,:], color=colors(i/num_terms))
    axs[1].plot(ddEs[i,:], color=colors(i/num_terms))
plt.show()
