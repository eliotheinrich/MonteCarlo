import os
import sys
sys.path.append("../../")

import numpy as np
import matplotlib.pyplot as plt
from DisplayTexture import *
import re
from matplotlib import animation

def load_spin_series(fileprefix, num_items):
    s0 = load_texture(fileprefix + "0.txt")
    N = len(s0)
    layers = len(s0[0,0])
    spins = np.zeros((num_items, N, N, layers, 3))
    for i in range(num_items):
        filename = fileprefix + str(i) + ".txt"
        spins[i,:,:,:,:] = load_texture(filename)
                
    return spins

def make_animation(spins, fps, filename, layers=None):
    if layers is None:
        layers = [0]

    frn = len(spins)

    fig = plt.figure()
    ax = plt.gca()

    def animate(i):
        ax.cla()
        if len(layers) > 1:
            display_trigonal_spins(ax, spins[i,:,:,:,:], layers[0], color='r')
            display_trigonal_spins(ax, spins[i,:,:,:,:], layers[1], color='b')
        else:
            display_trigonal_spins(ax, spins[i,:,:,:,:], layers[0], color='k')


    ani = animation.FuncAnimation(
        fig,
        animate,
        frn,
        interval=1000/fps
    )

    my_writer = animation.PillowWriter(fps=fps, codec='libx264', bitrate=2)
    ani.save(filename, fps=fps)

    ax.cla()
    if len(layers) > 1:
        display_trigonal_spins(ax, spins[-1,:,:,:,:], layers[0], color='r')
        display_trigonal_spins(ax, spins[-1,:,:,:,:], layers[1], color='b')
    else:
        display_trigonal_spins(ax, spins[-1,:,:,:,:], layers[0], color='k')

    plt.savefig("final_texture.jpg")


spins = load_spin_series("texture/texture1-", 30)
make_animation(spins, 5, "animation1.gif", [0,1]) 
spins = load_spin_series("texture/texture2-", 30)
make_animation(spins, 5, "animation2.gif", [0,1]) 




