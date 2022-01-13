import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import os

sys.path.append('..')
from DisplayTexture import *

def make_animation(spins, layer, fps, filename):
    frn = len(spins)

    fig = plt.figure()
    ax = plt.gca()

    def animate(i):
        print(i)
        ax.cla()
        display_trigonal_spins(ax, spins[i], layer, color='r')
        display_trigonal_spins(ax, spins[i], layer+1, color='b')



    ani = animation.FuncAnimation(
        fig,
        animate,
        frn,
        interval=1000/fps
    )

    my_writer = animation.PillowWriter(fps=fps, codec='libx264', bitrate=2)
    ani.save(filename, fps=fps)

def get_idx(filename):
    i = filename.index('.')
    return int(filename[5:i])

if __name__ == "__main__":
    cwd = os.getcwd()
    os.chdir('data')
    spins = []
    filenames = os.listdir()
    filenames = sorted(filenames, key = lambda x: get_idx(x))

    for f in filenames:
        spins.append(load_texture(f))

    spins = np.array(spins)
    os.chdir(cwd)

    make_animation(spins, 0, 25, "animation2.gif")

