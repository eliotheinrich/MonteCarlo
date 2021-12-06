import os
import sys
sys.path.append("../")

import numpy as np
import matplotlib.pyplot as plt
from DisplayTexture import display_hexagonal_spins
import re
from matplotlib import animation

def load_spin_series(filename):
    num_items = len(open(filename).readlines(  ))
    with open(filename, 'r') as f:
        line = f.readline()
        rows = re.findall("\[(.*?)\]", line)
        N = len(rows)

    spins = np.zeros((num_items, N, N, 2, 3))
    with open(filename, 'r') as f:
        for line_num,line in enumerate(f.readlines()):
            rows = re.findall("\[(.*?)\]", line)
            for i,row in enumerate(rows):
                unitcells = row.split('\t')
                for j, cell in enumerate(unitcells):
                    atoms = re.findall("\((.*?)\)", cell)
                    spins[line_num, i, j, 0, :] = np.array(atoms[0].split(' '),dtype=float)
                    spins[line_num, i, j, 1, :] = np.array(atoms[1].split(' '),dtype=float)
                

    return spins

def make_animation(spins, fps, filename):
    frn = len(spins)

    fig = plt.figure()
    ax = plt.gca()

    def animate(i):
        ax.cla()
        display_hexagonal_spins(ax, spins[i,:,:,:,:], color=False)


    ani = animation.FuncAnimation(
        fig,
        animate,
        frn,
        interval=1000/fps
    )

    my_writer = animation.PillowWriter(fps=fps, codec='libx264', bitrate=2)
    ani.save(filename, fps=fps)

    ax.cla()
    display_hexagonal_spins(ax, spins[-1,:,:,:,:])
    plt.savefig("final_texture.jpg")


spins = load_spin_series("Log.txt")
make_animation(spins, 20, "animation.gif") 




