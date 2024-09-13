import subprocess
import re
from glob import glob

import numpy as np
from matplotlib import pyplot as plt


output_file_path_prefix = "outputmyo2/output_orca_scan_"

input_file_list = glob(output_file_path_prefix + "*")
print(input_file_list)
print(len(input_file_list))
num_files = len(input_file_list)


energy_line_prefix = "FINAL SINGLE POINT ENERGY (QM/MM)"

energy_list = []

for i in range(num_files):
    output_file_path = output_file_path_prefix + f"{i:03d}"

    found = False

    with open(output_file_path, "r") as file:
        lines = file.readlines()
    file.close()

    for line in lines:
        if energy_line_prefix in line:
            line = float(line.split()[-1])
            energy_list.append((line))
            found = True

    if not found:
        energy_list.append(100000)



min_elem = min(energy_list)
minindex = np.argmin(energy_list)

energy_limit = 20000
for i in range(len(energy_list)):
    if energy_list[i] > energy_limit:
        energy_list[i] = energy_limit
    else:
        energy_list[i] -= min_elem

setter = 0.7

while max(energy_list) > setter:
    energy_list[energy_list.index(max(energy_list))] = setter

print(max(energy_list))


#h2_pos=h2_pos_list, fe_pos=iron_atom_coords, dist_h2_FE=dist_h2_FE, dist_h2=dist_h2
coordfile = np.load("coord-dist_save_file.npz")
fe_pos = coordfile['fe_pos']

distH2 = coordfile['dist_h2']
# distH2Fe = coordfile['dist_h2_FE'][::-1]
h2pos = coordfile['h2_pos']
distH2Fe = []
for i, hpos in enumerate(h2pos):
    c = np.linalg.norm(fe_pos-(hpos[0]+hpos[1])/2)
    distH2Fe.append(c)




print(len(distH2), len(distH2Fe), len(energy_list))
print(coordfile['h2_pos'][minindex])
print(distH2[minindex], distH2Fe[minindex])

x = 2
x_accuracy = 0.1
y = 2
y_accuracy = 0.07



print(sorted(set(np.round(distH2, decimals=6))))
print(sorted(set(np.round(distH2Fe, decimals=4))))

print("new pos (not min)")

for i in range(len(h2pos)):
    if x-x_accuracy < distH2[i] < x+x_accuracy and y - y_accuracy < distH2Fe[i] < y + y_accuracy:
        print(h2pos[i])
        print(distH2[i], distH2Fe[i])

MEDIUM_SIZE = 18
plt.rc('font', size=MEDIUM_SIZE)  # controls default text sizes
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labe
plt.tricontourf(distH2, distH2Fe, np.array(energy_list)*627.509, 42, cmap="inferno_r")
plt.plot([0.8], [1.7], 'x', ms=7, markeredgewidth=3, color='k')
plt.plot([2.6], [1], 'x', ms=7, markeredgewidth=3, color='blue')
plt.plot([2], [2], 'x', ms=7, markeredgewidth=3, color='green')
plt.colorbar(label='\nEnergy (kcal/mol)')
plt.ylabel("Distance COM(H$_2$) to Fe", fontsize=20)
plt.xlabel(r"Bond Distance H$_2$ ($\mathrm{\AA}$)", fontsize=20)

plt.show()