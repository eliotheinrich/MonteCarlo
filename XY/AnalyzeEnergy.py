import sys
import numpy as np
import matplotlib.pyplot as plt

Es = []

filename = sys.argv[1]
with open(filename) as f:
    for line in f.readlines():
        Es.append(float(line.strip()))

Es = np.array(Es)

plt.plot(Es)
plt.title("Energy")
plt.show()


