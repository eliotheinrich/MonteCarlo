import numpy as np
import matplotlib.pyplot as plt

Ms = []
As = []
Es = []

with open('Log.txt') as f:
    for n,line in enumerate(f.readlines()):
        data = line.split('\t')
        Ms.append(np.array(data[0].split(' '), dtype=float))
        As.append(float(data[1]))
        Es.append(float(data[2]))

Ms = np.array(Ms)
As = np.array(As)
Es = np.array(Es)

ax = plt.gca()
plt.plot(As)
plt.title("Acceptance rate")
plt.xlabel("Iterations")
plt.ylabel("%")
plt.show()
'''
plt.plot(Ms[:,0], label='x')
plt.plot(Ms[:,1], label='y')
plt.plot(Ms[:,2], label='z')
plt.title("Magnetization")
plt.show()
'''

plt.plot(Es[len(Es)//10:])
plt.title("Energy")
plt.show()


