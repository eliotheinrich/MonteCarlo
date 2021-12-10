import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from matplotlib import animation

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



def load_stiffness_data(filename):
    i = filename.index('.')
    L = int(filename[9:i])

    with open(filename) as f:
        s = f.readline().split('\t')
        res = int(s[0])
        num_samples = int(s[1])
    
    Ts = np.zeros(res)
    rhos = np.zeros(res)
    dEs = np.zeros((res, num_samples))
    ddEs = np.zeros((res, num_samples))

    with open(filename) as f:
        f.readline()
        for i in range(res):
            line = f.readline()
            data = np.array(line.split('\t'))
            Ti = float(data[0])
            data = data[1:-1]

            for n,x in enumerate(data):
                dE, ddE = x.split(',')
                dEs[i, n] = float(dE[1:])
                ddEs[i, n] = float(ddE[:-1])



            Ts[i] = Ti

    return L, Ts, dEs, ddEs

if __name__ == "__main__":
    cwd = os.getcwd()
    os.chdir(sys.argv[1])

    fs = [f for f in os.listdir() if f[:9] == "stiffness"]
    _, _, dEs, _ = load_stiffness_data(fs[0])

    N = len(fs)
    res, num_samples = dEs.shape

    Ts = np.zeros(res)
    ρs = np.zeros((N, res))
    dEs = np.zeros((N, res, num_samples))
    ddEs = np.zeros((N, res, num_samples))
    Ls = np.zeros(N, dtype=int)

    for n, f in enumerate(fs):
        if f[:9] == 'stiffness':
            i = f.index('.')
            L, Ts, dE, ddE = load_stiffness_data(f)
            Ls[n] = L
            dEs[n,:,:] = dE
            ddEs[n,:,:] = ddE
            ρs[n,:] = (np.mean(ddE, axis=1) - np.std(dE, axis=1)**2/Ts)/(L**2)

    inds = np.argsort(Ls)
    ρs = ρs[inds]
    Ls = Ls[inds]
    L0s = np.sin(np.linspace(0, 2*np.pi, 50)) + 1.1

    os.chdir(cwd)

    make_animation(Ts, ρs, Ls, L0s, 15, "animation.gif") 






