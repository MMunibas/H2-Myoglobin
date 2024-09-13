import os
import numpy as np
from glob import glob

# MDAnalysis
import MDAnalysis

import matplotlib.pyplot as plt

from ase import Atoms
from ase import io
from ase.io import read


coordinates = np.load("coordfile.npz")

h2_positions = coordinates['h2']
iron_pos = coordinates['iron']
iron_dist = coordinates['FE_dist']
h2_dist = coordinates['h2_dist']
fe_upper = 1.1
fe_lower = 0.9
h2_upper = 2.9
h2_lower = 2.0


for i in range(len(h2_positions)):
    if iron_dist[i] >= fe_lower and iron_dist[i] < fe_upper:
        if h2_dist[i] < h2_upper and h2_dist[i] > h2_lower:
            print(h2_positions[i], "h2pos")
            print(iron_pos, "ironpos")
            print(iron_dist[i], "irondist")
            print(h2_dist[i], "h2dist")
            print(i)


