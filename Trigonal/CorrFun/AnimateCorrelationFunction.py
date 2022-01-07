import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
import os
from AnalyzeCorrelationFunction import *

def make_animation(T, C, layer, fps, filename, kspace=False):
    frn = len(T)

    fig = plt.figure()
    ax = plt.gca()

    def animate(i):
        print(i)
        ax.cla()

        if kspace:
            plot_structure_factor(C[i,:,:,layer], ax)
        else:
            plot_correlation_function(C[i,:,:,layer], ax)
        ax.set_title(f'T = {T[i]}')



    ani = animation.FuncAnimation(
        fig,
        animate,
        frn,
        interval=1000/fps
    )

    my_writer = animation.PillowWriter(fps=fps, codec='libx264', bitrate=2)
    ani.save(filename, fps=fps)

if __name__ == "__main__":
    cwd = os.getcwd()
    os.chdir(sys.argv[1])
    T, C, dC = load_correlation_data("correlation.txt")
    os.chdir(cwd)

    N = C.shape[1]
    L = C.shape[-1]

    make_animation(T, C, L//2, 10, "animation_real.gif", kspace=False)

