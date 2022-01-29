import numpy as np
import matplotlib.pyplot as plt

def load_data(filename):
    with open(filename) as f:
        N = int(f.readline())
        print(N)
        data = f.readline()[:-1].strip()
        data = np.array([float(x) for x in data.split('\t')])
        data = data.reshape((N, N))

        return N, data


N, E = load_data('density.txt')

x = np.linspace(0, N, N)
X, Y = np.meshgrid(x, x)
plt.pcolor(X, Y, E)
plt.show()


