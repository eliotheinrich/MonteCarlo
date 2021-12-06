import os
import sys

sys.path.append(os.environ['RESEARCH_PATH'] + '/ChromiumTrihalides/Code')

from ibcc import *

params = params_skyrmionb


params_s = ''

for p in params.values():
    params_s += str(p) + ' '

N = 3
nsteps = 2000000

params_s += str(N) + ' '
params_s += str(nsteps)


os.system("g++ -std=c++11 -O3 GenerateData.cpp -o GenerateData")

num_runs = 1
for i in range(num_runs):
    os.system("./GenerateData " + params_s)
    os.system("python DisplayTexture.py SpinTexture.txt")

#os.system("rm GenerateData")
