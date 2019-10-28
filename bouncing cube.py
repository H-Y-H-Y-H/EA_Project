import numpy as np
from numba import jit, njit
import time

# Parameters
springs = {}
masses = {}


class Mass:
    def __init__(self, mass, x, y, z):
        self.m = mass
        self.p = [x, y, z]
        self.v = [0] * 3
        self.a = [0] * 3

    # def acceleration_2_velocity(self, v0):


class Spring:
    def __init__(self, constraint, length, mass1, mass2):
        self.k = constraint
        self.lgh = length
        self.cnn1 = mass1
        self.cnn2 = mass2


def init_cube():
    for i in range(8):
        x = y = z = 0
        if i > 3:
            x = 0.1
        if i % 4 > 1:
            y = 0.1
        if i % 2 == 1:
            z = 0.1
        masses[i] = Mass(0.1, x, y, z)
        # print(masses[i].p)

    k = 0
    for i in range(0, 8):
        for j in range(i + 1, 8):
            springs[k] = Spring(10000, 0.1, i, j)
            # print(k, springs[k].cnn1, springs[k].cnn2)
            k += 1





