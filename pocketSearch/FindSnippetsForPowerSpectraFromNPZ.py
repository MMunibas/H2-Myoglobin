import os
import numpy as np
from glob import glob

# MDAnalysis
import MDAnalysis

import matplotlib.pyplot as plt

from ase import Atoms
from ase import io
prefix = "pocketArray/v14.H2.eval/1/"
distfile = "pocket_arrays_v14.H2_1.DCD_files.1-37.H2_num1.npz"
h2_read = 1
min_intervall_length = 20000
angstroemcap = 5


distfile = prefix + distfile
# 2 4 143000 165000
#
# 2 4 170000 180000

print("reading data")

data = np.load(distfile)
pocket_arr = data['pocket_arr']
pocket_XE1_dist_arr = data['pocket_XE1']
pocket_XE2_dist_arr = data['pocket_XE2']
pocket_XE3_dist_arr = data['pocket_XE3']
pocket_XE4_dist_arr = data['pocket_XE4']
pocket_XE5_dist_arr = data['pocket_XE5']
pocket_XE6_dist_arr = data['pocket_XE6']
pocket_XE7_dist_arr = data['pocket_XE7']
pocket_XE8_dist_arr = data['pocket_XE8']
pocket_XE9_dist_arr = data['pocket_XE9']

pockets = [pocket_XE1_dist_arr,pocket_XE2_dist_arr,pocket_XE3_dist_arr,pocket_XE4_dist_arr,pocket_XE5_dist_arr,pocket_XE6_dist_arr,pocket_XE7_dist_arr,pocket_XE8_dist_arr,pocket_XE9_dist_arr]
unknown_pocket = np.zeros(len(pocket_XE8_dist_arr))
chunksize = 1000
for x in range(len(pocket_XE2_dist_arr)//chunksize):
    in_known_pocket = False
    for pocket in pockets:
        if np.sum(pocket[x*chunksize:(x+1)*chunksize])/chunksize <= 3:
            in_known_pocket = True

    if not in_known_pocket:
        unknown_pocket[x*chunksize:(x+1)*chunksize] = 1.5
        pocket_arr[x*chunksize:(x+1)*chunksize] = 10
    else:
        unknown_pocket[x*chunksize:(x+1)*chunksize] = 15

oldmindex = 0
cnt= 0
for x in range(len(pocket_XE8_dist_arr)):
        mindex = np.argmin([pocket[x] for pocket in pockets])
        if mindex == oldmindex and pockets[mindex][x] < angstroemcap:
            cnt += 1
        else:
            if cnt > min_intervall_length:
                print(str(oldmindex+1) + " " + str(h2_read) +  " " + str(x-cnt) + " " + str(x))
            oldmindex = mindex
            cnt = 0

print(unknown_pocket)

num_timesteps = len(pocket_arr)

max_val = np.max([
    np.max(pocket_XE1_dist_arr), np.max(pocket_XE2_dist_arr), np.max(pocket_XE3_dist_arr), np.max(pocket_XE4_dist_arr),
    np.max(pocket_XE5_dist_arr), np.max(pocket_XE6_dist_arr), np.max(pocket_XE7_dist_arr)])

max_val = min(max_val, 20)

x = np.arange(num_timesteps)

fig, (ax1, ax2) = plt.subplots(2, 1)

# plot in which pocket the H2 is
nplot = 100
ax1.plot(x[::nplot], pocket_arr[::nplot], linewidth=2.0, label='pocket state')

# plot all distances to H2 for debugging
ax2.plot(x[::nplot], pocket_XE1_dist_arr[::nplot], linewidth=1.0, label='dist to XE1')
ax2.plot(x[::nplot], pocket_XE2_dist_arr[::nplot], linewidth=1.0, label='dist to XE2')
ax2.plot(x[::nplot], pocket_XE3_dist_arr[::nplot], linewidth=1.0, label='dist to XE3')
ax2.plot(x[::nplot], pocket_XE4_dist_arr[::nplot], linewidth=1.0, label='dist to XE4')
ax2.plot(x[::nplot], pocket_XE5_dist_arr[::nplot], linewidth=1.0, label='dist to XE5')
ax2.plot(x[::nplot], pocket_XE6_dist_arr[::nplot], color='black', linewidth=2.0, label='dist to XE6')
ax2.plot(x[::nplot], pocket_XE7_dist_arr[::nplot], color='gray', linewidth=1.0, label='dist to XE7')
ax2.plot(x[::nplot], pocket_XE8_dist_arr[::nplot], color='gold', linewidth=1.0, label='dist to XE8')
ax2.plot(x[::nplot], pocket_XE9_dist_arr[::nplot], color='yellow', linewidth=1.0, label='dist to XE9')
ax2.plot(x[::nplot], unknown_pocket[::nplot], color='green', linewidth=1.0, label='UnknownPocket')
ax2.legend()

print("terminated correctly")

plt.show()


