import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from matplotlib import animation
from AnalyzeStiffnessCurve import load_stiffness_data

def make_animation(Ts, ρs, Ls, L0s, fps, filename):
    frn = len(L0s)

    fig = plt.figure()
    ax = plt.gca()
    T_KT = 0.88

    def animate(i):
        print(i)
        ax.cla()
        for n, L in enumerate(Ls):
            plt.plot(Ts, ρs[n]/(1. + 1./(2.*np.log(Ls[n]/L0s[i]))), marker='o', label=f'L = {L}')

        ax.plot(Ts, 2./np.pi*Ts, 'k--', label=r'$2T/\pi$')
 #       ax.axvline(T_KT, linestyle='-', color='k', alpha=0.5, label=r'$T_{KT}$')

        ax.legend()
        ax.set_xlim(0., max(Ts))
        ax.set_ylim(-0.05, 1.1)
        ax.set_xlabel(r'$T$', fontsize=15)
        ax.set_ylabel(r'$\Upsilon(L)/(1 + (2\log(L/L_0))^{-1})$', fontsize=15)
        ax.set_title(r'$L_0 = $' + f'{L0s[i]:.3f}')



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
    os.chdir('data')

    fs = [f for f in os.listdir() if f[:9] == "stiffness"]
    _, Ts, _, _ = load_stiffness_data(fs[0])

    N = len(fs)
    res = len(Ts)

    Ts = np.zeros(res)
    ρs = np.zeros((N, res))
    dEs = np.zeros((N, res))
    ddEs = np.zeros((N, res))
    Ls = np.zeros(N, dtype=int)

    for n, f in enumerate(fs):
        if f[:9] == 'stiffness':
            i = f.index('.')
            L, Ts, ρs[n,:], _ = load_stiffness_data(f)
            Ls[n] = L

    inds = np.argsort(Ls)
    ρs = ρs[inds]
    Ls = Ls[inds]
    L0s = np.sin(np.linspace(0, 2*np.pi, 50)) + 1.1

    os.chdir(cwd)

    make_animation(Ts, ρs, Ls, L0s, 15, "animation.gif") 






