import numpy as np
from numpy.linalg import norm
import re
import matplotlib.pyplot as plt
import numpy.linalg as la
import sys


# Load spins from texture file
def load_texture(filename):
    with open(filename, 'r') as f:
        firstline = f.readline()
        (num_sublattices, N) = [int(i.strip()) for i in firstline.split('\t')]
        spins = np.zeros((N,N,num_sublattices,3))

        for i,line in enumerate(f.readlines()):
            data1 = re.findall("\((.*?)\)", line)
            data2 = [re.findall("\[(.*?)\]", unitcell) for unitcell in data1]
            for j in range(N):
                for s in range(num_sublattices):
                    spins[i,j,s,:] = np.array(data2[j][s].split(' '), dtype=float)

        return spins

# Plot hexagonal spin texture
def display_hexagonal_spins(ax, spins, color=True):
    N = len(spins)

    # NN vectors
    α1 = np.array([-np.sqrt(3)/2., -1/2.])
    α2 = np.array([np.sqrt(3)/2., -1/2.]) 
    α3 = np.array([0., 1.])

    # NNN vectors
    β1 = np.array([np.sqrt(3),0])
    β2 = np.array([np.sqrt(3)/2,3/2])
    β3 = np.array([-np.sqrt(3)/2,3/2])



    xs0 = np.zeros(N**2)
    ys0 = np.zeros(N**2)
    spinxs0 = np.zeros(N**2)
    spinys0 = np.zeros(N**2)
    spinzs0 = np.zeros(N**2)

    xs1 = np.zeros(N**2)
    ys1 = np.zeros(N**2)
    spinxs1 = np.zeros(N**2)
    spinys1 = np.zeros(N**2)
    spinzs1 = np.zeros(N**2)
    for i in range(N):
        for j in range(N):
            (x1,y1) = i*β1 + j*β2
            (x2,y2) = i*β1 + j*β2 + α3

            xs0[i*N + j] = x1
            ys0[i*N + j] = y1
            xs1[i*N + j] = x2
            ys1[i*N + j] = y2

            spin0 = spins[i,j,0,:]
            spin1 = spins[i,j,1,:]
            spinxs0[i*N + j] = spin0[0]
            spinys0[i*N + j] = spin0[1]
            spinzs0[i*N + j] = spin0[2]
            spinxs1[i*N + j] = spin1[0]
            spinys1[i*N + j] = spin1[1]
            spinzs1[i*N + j] = spin1[2]



    scale = 0.75
    L = len(xs0)
    width = .003

    if color:
        red = np.array([1,0,0])
        blue = np.array([0,0,1])
        for i in range(L):
            s0 = (spinzs0[i] + 1)/2
            s1 = (spinzs1[i] + 1)/2
            color0 = s0*red + (1 - s0)*blue
            color1 = s1*red + (1 - s1)*blue

            ax.quiver(xs0[i], ys0[i], spinxs0[i], spinys0[i], scale_units='x', scale=scale, width=width, color=color0)
            ax.quiver(xs1[i], ys1[i], spinxs1[i], spinys1[i], scale_units='x', scale=scale, width=width, color=color1)
    else:
        ax.quiver(xs0, ys0, spinxs0, spinys0, scale_units='x', scale=scale, width=width)
        ax.quiver(xs1, ys1, spinxs1, spinys1, scale_units='x', scale=scale, width=width)

    plt.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect('equal')


if __name__ == "__main__":
    filename = sys.argv[1]
    spins = load_texture(filename)
    fig = plt.figure()
    ax = fig.gca()
    display_hexagonal_spins(ax, spins)
    fig.patch.set_alpha(0.)
    ax.patch.set_alpha(0.)
    fig.set_size_inches(6, 6)
    plt.show()






