import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from AnalyzeStiffnessCurve import *
from matplotlib import animation

def make_animation(T, U2, L, L0, fps, filename):
    frn = len(L0)

    fig = plt.figure()
    ax = plt.gca()

    def animate(i):
        print(i)
        ax.cla()
        for Ln, U2n in zip(L, U2):
            plt.plot(T, U2n/(1. + 1./(2.*np.log(Ln/L0[i]))), marker='o', label=f'L = {L[n]}')

        ax.plot(T, 2./np.pi*T, 'k--', label=r'$2T/\pi$')
 #       ax.axvline(T_KT, linestyle='-', color='k', alpha=0.5, label=r'$T_{KT}$')

        ax.legend()
        ax.set_xlim(0., max(T))
        ax.set_ylim(-0.05, 1.1)
        ax.set_xlabel(r'$T$', fontsize=15)
        ax.set_ylabel(r'$\Upsilon(L)/(1 + (2\log(L/L_0))^{-1})$', fontsize=15)
        ax.set_title(r'$L_0 = $' + f'{L0[i]:.3f}')



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

    filenames = [f for f in os.listdir() if f[:9] == "stiffness"]
    _, T, _, _ = load_stiffness_data(filenames[0])

    N = len(filenames)
    res = len(T)

    U2 = np.zeros((N, res))
    U4 = np.zeros((N, res))
    L = np.zeros(N, dtype=int)

    for n, f in enumerate(filenames):
        if f[:9] == 'stiffness':
            i = f.index('.')
            L[n], T, U2[n], U4[n] = load_stiffness_data(f)

    inds = np.argsort(L)
    U2 = U2[inds]
    U4 = U4[inds]
    L = L[inds]
    L0 = np.sin(np.linspace(0, 2*np.pi, 50)) + 1.1

    os.chdir(cwd)

    make_animation(T, U2, L, L0, 15, "animation.gif") 






