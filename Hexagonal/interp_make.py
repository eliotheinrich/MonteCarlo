import os
import sys

sys.path.append(os.environ['RESEARCH_PATH'] + '/ChromiumTrihalides/Code')

from ibcc import *

p1 = params_phys.copy()
p2 = params_phys.copy()

set_NNN_DMI(p1, 3, True)
set_NNN_DMI(p2, 6, True)

p1 = np.array(list(p1.values()))
p2 = np.array(list(p2.values()))

def interp(p1, p2, N):
    p1 = np.array(p1)
    p2 = np.array(p2)
    ps = np.zeros((N, len(p1)))

    for n,i in enumerate(ps):
        ps[n,:] = p1 + (p2 - p1)*n/N

    return ps


ps = interp(p1, p2, 5)





N = 20
nsteps = 50000000



os.system("g++ -std=c++11 -O3 GenerateData.cpp -o GenerateData")

for p in ps:
    params_s = ''

    for val in p:
        params_s += str(val) + ' '

    params_s += str(N) + ' '
    params_s += str(nsteps)

    os.system("./GenerateData " + params_s)
    os.system("python DisplayTexture.py SpinTexture.txt")

os.system("rm GenerateData")

