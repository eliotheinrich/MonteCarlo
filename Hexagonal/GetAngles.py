from DisplayTexture import *


# Load spins from texture file
spins = load_texture("SpinTexture.txt")
spinA = np.array([0,0,0], dtype=float)
spinB = np.array([0,0,0], dtype=float)


N = len(spins)

# Get average magnetization by sublattice
for i in range(N):
    for j in range(N):
        spinA += spins[i,j,0,:]
        spinB += spins[i,j,1,:]

spinA = spinA/(N**2)
spinB = spinB/(N**2)


# Compute equilibrium angles
thetaA = np.arccos(spinA[2])
thetaB = np.arccos(spinB[2])
PA = np.arctan(spinA[1]/spinA[0])
PB = np.arctan(spinB[1]/spinB[0])

print(thetaA, thetaB, PA, PB)






