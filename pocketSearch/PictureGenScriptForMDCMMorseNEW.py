import numpy as np
import matplotlib.pyplot as plt

distfile = ("pocketArray/v14.5_H2.eval/" + 'pocket_arrays_v14.5_H2.DCD_files.1-37.H2_num2.npz')

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


num_timesteps = len(pocket_arr)

max_val = np.max([
    np.max(pocket_XE1_dist_arr), np.max(pocket_XE2_dist_arr), np.max(pocket_XE3_dist_arr), np.max(pocket_XE4_dist_arr),
    np.max(pocket_XE5_dist_arr), np.max(pocket_XE6_dist_arr), np.max(pocket_XE7_dist_arr)])

max_val = min(max_val, 20)

x = np.arange(num_timesteps)

# plot in which pocket the H2 is
nplot = 100
# ax1.plot(x[::nplot], pocket_arr[::nplot], linewidth=2.0, label='pocket state')



x = np.arange(num_timesteps)


SMALL_SIZE = 15
MEDIUM_SIZE = 22
BIGGER_SIZE = 28

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

# plot in which pocket the H2 is
# plt.plot(x, pocket_arr, linewidth=2.0, label='pocket state')
fig = plt.figure(figsize=(16, 9))
x = x/1000
# plot all distances to H2 for debugging
plt.plot(x[::nplot], pocket_XE1_dist_arr[::nplot], linewidth=1.0, label='dist to XE1')
plt.plot(x[::nplot], pocket_XE2_dist_arr[::nplot], linewidth=1.0, label='dist to XE2')
plt.plot(x[::nplot], pocket_XE3_dist_arr[::nplot], linewidth=1.0, label='dist to XE3')
plt.plot(x[::nplot], pocket_XE4_dist_arr[::nplot], linewidth=1.0, label='dist to XE4')
plt.plot(x[::nplot], pocket_XE5_dist_arr[::nplot], linewidth=1.0, label='dist to XE5')
plt.plot(x[::nplot], pocket_XE6_dist_arr[::nplot], color='black', linewidth=2.0, label='dist to XE6')
plt.plot(x[::nplot], pocket_XE7_dist_arr[::nplot], color='gray', linewidth=1.0, label='dist to XE7')
plt.plot(x[::nplot], pocket_XE8_dist_arr[::nplot], color='gold', linewidth=1.0, label='dist to XE8')
plt.plot(x[::nplot], pocket_XE9_dist_arr[::nplot], color='yellow', linewidth=1.0, label='dist to XE9')
# plt.xlim(0,1000)
plt.legend(fancybox=True, loc='upper right', title='Distance H$_2$\nto Pocket', framealpha=1)
plt.xlabel("Time (ps)")
plt.ylabel("Distance ($\mathrm{\AA}$)")
plt.savefig(dpi=400, fname="CGenFF", )
plt.show()

