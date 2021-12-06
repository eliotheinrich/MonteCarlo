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
def display_trigonal_spins(ax, spins, layer, color=None):
    N = len(spins)

    # NN vectors
    α1 = np.array([-np.sqrt(3)/2., -1/2.])
    α2 = np.array([np.sqrt(3)/2., -1/2.]) 
    α3 = np.array([0., 1.])

    # NNN vectors
    β1 = np.array([np.sqrt(3),0])
    β2 = np.array([np.sqrt(3)/2,3/2])
    β3 = np.array([-np.sqrt(3)/2,3/2])



    xs = np.zeros(N**2)
    ys = np.zeros(N**2)
    spinxs = np.zeros(N**2)
    spinys = np.zeros(N**2)
    spinzs = np.zeros(N**2)

    for i in range(N):
        for j in range(N):
            (x,y) = i*β1 + j*β2
            xs[i*N + j] = x
            ys[i*N + j] = y


            spin = spins[i,j,layer,:]
            spinxs[i*N + j] = spin[0]
            spinys[i*N + j] = spin[1]
            spinzs[i*N + j] = spin[2]

            if i < N-1:
                line(ax, [x,y], np.array([x,y]) + β1, 'k--', alpha=0.5, linewidth=0.5)
            if j < N-1:
                line(ax, [x,y], np.array([x,y]) + β2, 'k--', alpha=0.5, linewidth=0.5)
            if i < N-1 and j > 0:    
                line(ax, [x,y], np.array([x,y]) + β1 - β2, 'k--', alpha=0.5, linewidth=0.5)




    scale = 0.75
    width = .003

    if color is None:
        color = 'k'
    ax.quiver(xs, ys, spinxs, spinys, scale_units='x', scale=scale, width=width, color=color)

    plt.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect('equal')


if __name__ == "__main__":
    filename = sys.argv[1]
    spins = load_texture(filename)
    fig = plt.figure()
    ax = fig.gca()
    display_trigonal_spins(ax, spins, 0, 'red')
#    display_trigonal_spins(ax, spins, 1, 'blue')
#    fig.patch.set_alpha(0.)
#    ax.patch.set_alpha(0.)
    fig.set_size_inches(6, 6)
    plt.show()






