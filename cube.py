import numpy as np

springs = np.zeros(28)
masses = np.zeros(8)

each_mass = [0.1] * len(masses)


def position(num_mass):
    pos = np.zeros([num_mass, 3])
    for i in range(num_mass):
        if i > 3:
            pos[i][0] = 0.1
        if i % 4 > 1:
            pos[i][1] = 0.1
        if i % 2 == 1:
            pos[i][2] = 0.1
    return pos


print(position(8))
