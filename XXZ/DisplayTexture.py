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
        (N1, N2, N3, sl) = np.array([int(i.strip()) for i in firstline.split('\t')])

        spins = np.array([float(x) for x in f.readline().split('\t')])
        spins = spins.reshape((N1*N2*N3*sl, 3))

        return spins.reshape((N1, N2, N3, sl, 3), order='F')

# Detect spin vortices
def find_vortices(spins, layer):
    N = len(spins)

    α1 = np.array([1,0])
    α2 = np.array([0,1])

    vortex_x = []
    vortex_y = []
    
    φ = np.zeros((N,N))

    for i in range(N):
        for j in range(N):
            φ[i,j] = np.arctan2(spins[i,j,layer,1], spins[i,j,layer,0])
            #print(φ[i,j])

    #print(np.max(φ))
    #print(np.min(φ))

    for i in range(N):
        for j in range(N):
            (x, y) = (i+0.5)*α1 + (j+0.5)*α2
            winding =((φ[i,j] - φ[(i+1)%N,j]) + (φ[(i+1)%N,j] - φ[(i+1)%N,(j+1)%N]) \
             + (φ[(i+1)%N,(j+1)%N] - φ[i,(j+1)%N]) + (φ[i,(j+1)%N] - φ[i,j]))
            if ((φ[i,j] - φ[(i+1)%N,j]) + (φ[(i+1)%N,j] - φ[(i+1)%N,(j+1)%N]) \
             + (φ[(i+1)%N,(j+1)%N] - φ[i,(j+1)%N]) + (φ[i,(j+1)%N] - φ[i,j]))%(2*np.pi) != 0:
                vortex_x.append(x)
                vortex_y.append(y)

    return np.array(vortex_x,dtype=float), np.array(vortex_y, dtype=float)


# Plot square XY spin texture
def display_trigonal_spins(ax, spins, layer, color=None):
    N = len(spins)

    # NN vectors
    α1 = np.array([1,0])
    α2 = np.array([0,1])



    xs = np.zeros(N**2)
    ys = np.zeros(N**2)
    spinxs = np.zeros(N**2)
    spinys = np.zeros(N**2)


    for i in range(N):
        for j in range(N):
            (x,y) = i*α1 + j*α2
            xs[i*N + j] = x
            ys[i*N + j] = y


            spin = spins[i,j,layer,:]
            spinxs[i*N + j] = spin[0]
            spinys[i*N + j] = spin[1]

            line(ax, [x,y], np.array([x,y]) + α1, 'k--', alpha=0.5, linewidth=0.5)
            line(ax, [x,y], np.array([x,y]) + α2, 'k--', alpha=0.5, linewidth=0.5)
            line(ax, [x,y], np.array([x,y]) - α1, 'k--', alpha=0.5, linewidth=0.5)
            line(ax, [x,y], np.array([x,y]) - α2, 'k--', alpha=0.5, linewidth=0.5)




    scale = 0.75
    width = .003

    if color is None:
        color = 'k'
    ax.quiver(xs, ys, spinxs, spinys, scale_units='x', scale=scale, width=width, color=color)
    vortex_x, vortex_y = find_vortices(spins, layer)
#    plt.plot(vortex_x, vortex_y, 'b.', markersize=3)

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
    fig.patch.set_alpha(0.)
    ax.patch.set_alpha(0.)
    fig.set_size_inches(6, 6)
    plt.show()






