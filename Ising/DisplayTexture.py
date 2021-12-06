import numpy as np
from plots import line
from numpy.linalg import norm
import re
import matplotlib.pyplot as plt
import numpy.linalg as la
import sys


# Load spins from texture file
def load_texture(filename):
    with open(filename, 'r') as f:
        firstline = f.readline()
        N1, N2, N3 = np.array([int(i.strip()) for i in firstline.split('\t')])
        spins = np.zeros((N1,N2,N3))

        for i,line in enumerate(f.readlines()):
            data1 = line.split(';')
            data2 = [x.split(',')[:-1] for x in data1[:-1]]
            for j in range(N2):
                for k in range(N3):
                    spins[i,j,k] = float(data2[j][k])

        return spins



# Plot square Ising spin texture
def display_square_spins(ax, spins, layer, color=None):
    N = len(spins)

    # NN vectors
    α1 = np.array([1,0])
    α2 = np.array([0,1])

    for i in range(N):
        for j in range(N):
            (x,y) = i*α1 + j*α2

            line(ax, [x,y], np.array([x,y]) + α1, 'k--', alpha=0.5, linewidth=0.5)
            line(ax, [x,y], np.array([x,y]) + α2, 'k--', alpha=0.5, linewidth=0.5)
            line(ax, [x,y], np.array([x,y]) - α1, 'k--', alpha=0.5, linewidth=0.5)
            line(ax, [x,y], np.array([x,y]) - α2, 'k--', alpha=0.5, linewidth=0.5)

    X,Y = np.meshgrid(range(N), range(N))
    plt.pcolor(X, Y, spins[:,:,layer])

    plt.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect('equal')


if __name__ == "__main__":
    filename = sys.argv[1]
    spins = load_texture(filename)
    fig = plt.figure()
    ax = fig.gca()
    display_square_spins(ax, spins, int(sys.argv[2]))
    fig.set_size_inches(6, 6)
    plt.show()






