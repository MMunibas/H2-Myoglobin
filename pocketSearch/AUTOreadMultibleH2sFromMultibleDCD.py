import os
import numpy as np
from glob import glob

# MDAnalysis
import MDAnalysis

import matplotlib.pyplot as plt

from ase import Atoms
from ase import io


def calc_pocket_pos(pocket_indices):
    positions_pocket = frame._pos[pocket_indices]
    mean_position_pocket = np.mean(positions_pocket, axis=0)
    dist_H2_pocket = np.linalg.norm(mean_position_H2 - mean_position_pocket)
    return dist_H2_pocket


def read_pocket_list(res_list):
    pocket_indices = []
    for list_index in range(len(dcd_list)):
        for res in res_list:
            pocket_indices.append(dcd_list[list_index].select_atoms(res).ix[0])
    return np.array(pocket_indices)


path_psf = '../Data/v14/step3_pbcsetup.psf'

file_dcd_prefix = "../Data/v14/dyna11"
file_dcd_suffix = ".dcd"
start_dcd = 1
end_dcd = 1
num_files = end_dcd - start_dcd + 1
paths = []

for i in range(start_dcd, end_dcd + 1):
    paths.append(file_dcd_prefix + str(i) + file_dcd_suffix)

num_timesteps = 25000 * num_files

H2_res_list = [156]

for h2_index in range(1):

    H1_string = "resname H2 and resid " + str(H2_res_list[h2_index]) + " and name H1"
    H2_string = "resname H2 and resid " + str(H2_res_list[h2_index]) + " and name H2"

    distfile = "pocket_arrays_V10.DCD_files." + str(start_dcd) + "-" + str(end_dcd) + ".H2_num" + str(
        h2_index + 1) + ".npz"

    # ATOM   2533  H1  H2    156      41.957  35.323  35.289  1.00  0.00      H2
    # ATOM   2534  H2  H2    156      41.227  35.422  35.348  1.00  0.00      H2
    # ATOM   2535  H1  H2    157      41.293  39.779  34.422  1.00  0.00      H2
    # ATOM   2536  H2  H2    157      40.592  39.964  34.468  1.00  0.00      H2
    # ATOM   2537  H1  H2    158      33.230  37.348  37.951  1.00  0.00      H2
    # ATOM   2538  H2  H2    158      32.820  37.903  37.693  1.00  0.00      H2
    # ATOM   2539  H1  H2    159      25.680  38.170  34.168  1.00  0.00      H2
    # ATOM   2540  H2  H2    159      25.035  38.435  34.446  1.00  0.00      H2
    # ATOM   2541  H1  H2    160      50.607  36.386  37.612  1.00  0.00      H2
    # ATOM   2542  H2  H2    160      50.612  37.031  37.650  1.00  0.00      H2
    # ["resname H2 and resid 156 and name H1", "resname H2 and resid 156 and name H2"]

    H2_atoms = [H1_string, H2_string]

    pocket_XE1 = ["resname LEU and resid 89 and name C", "resname HSD and resid 93 and name C",
                  "resname LEU and resid 104 and name C", "resname PHE and resid 138 and name C",
                  "resname ILE and resid 142 and name C", "resname TYR and resid 146 and name C"]

    pocket_XE2 = ["resname LEU and resid 72 and name C", "resname ILE and resid 107 and name C",
                  "resname SER and resid 108 and name C", "resname LEU and resid 135 and name C",
                  "resname PHE and resid 138 and name C", "resname ARG and resid 139 and name C"]

    pocket_XE3 = ["resname TRP and resid 7 and name C", "resname LEU and resid 76 and name C",
                  "resname GLY and resid 80 and name C", "resname ALA and resid 134 and name C",
                  "resname LEU and resid 137 and name C", "resname PHE and resid 138 and name C"]

    pocket_XE4 = ["resname GLY and resid 25 and name C", "resname ILE and resid 28 and name C",
                  "resname LEU and resid 29 and name C", "resname GLY and resid 65 and name C",
                  "resname VAL and resid 68 and name C", "resname LEU and resid 72 and name C"]

    # original
    # pocket_XE5 = ["resname LEU and resid 29 and name C", "resname HEME and resid 1154 and name FE"]

    # more accurate replacement for XE5
    pocket_XE5 = ['resname VAL and resid 68 and name CA', 'resname LEU and resid 29 and name CA',
                  'resname GLY and resid 65 and name CA', 'resname ILE and resid 107 and name CA',
                  'resname HSD and resid 64 and name CA', 'resname LEU and resid 69 and name CA',
                  'resname THR and resid 67 and name CA', 'resname ILE and resid 28 and name CA',
                  'resname LEU and resid 104 and name CA', 'resname PHE and resid 43 and name CA',
                  'resname GLY and resid 25 and name CA', 'resname LEU and resid 72 and name CA',
                  'resname LEU and resid 32 and name CA', 'resname THR and resid 39 and name CA',
                  'resname ALA and resid 71 and name CA', 'resname ILE and resid 99 and name CA',
                  'resname HSD and resid 93 and name CA', 'resname LEU and resid 61 and name CA',
                  'resname LYS and resid 42 and name CA', 'resname LEU and resid 89 and name CA',
                  'resname TYR and resid 103 and name CA']

    # got replaced by XE9
    # pocket_XE6 = ['resname VAL and resid 10 and name CA', 'resname MET and resid 131 and name CA', 'resname ALA and resid 127 and name CA', 'resname GLN and resid 128 and name CA', 'resname PHE and resid 123 and name CA', 'resname VAL and resid 13 and name CA', 'resname LEU and resid 9 and name CA', 'resname ALA and resid 130 and name CA', 'resname TRP and resid 14 and name CA', 'resname LEU and resid 11 and name CA', 'resname HSD and resid 12 and name CA', 'resname LEU and resid 115 and name CA', 'resname GLY and resid 129 and name CA', 'resname ASP and resid 122 and name CA', 'resname GLY and resid 124 and name CA', 'resname ILE and resid 112 and name CA', 'resname ASN and resid 132 and name CA', 'resname ASP and resid 126 and name CA', 'resname TRP and resid 7 and name CA', 'resname GLN and resid 8 and name CA', 'resname GLU and resid 6 and name CA', 'resname HSD and resid 116 and name CA', 'resname ALA and resid 125 and name CA', 'resname ALA and resid 134 and name CA', 'resname ILE and resid 111 and name CA', 'resname PRO and resid 120 and name CA', 'resname LYS and resid 133 and name CA', 'resname HSD and resid 119 and name CA', 'resname LEU and resid 135 and name CA', 'resname LEU and resid 76 and name CA', 'resname VAL and resid 17 and name CA', 'resname ALA and resid 15 and name CA']
    pocket_XE6 = ['resname TRP and resid 14 and name CA', 'resname VAL and resid 13 and name CA',
                  'resname VAL and resid 10 and name CA', 'resname LEU and resid 11 and name CA',
                  'resname VAL and resid 17 and name CA', 'resname GLY and resid 73 and name CA',
                  'resname MET and resid 131 and name CA', 'resname HSD and resid 12 and name CA',
                  'resname ALA and resid 15 and name CA', 'resname LEU and resid 69 and name CA',
                  'resname LEU and resid 115 and name CA', 'resname LEU and resid 76 and name CA',
                  'resname LYS and resid 16 and name CA', 'resname ILE and resid 111 and name CA',
                  'resname LEU and resid 72 and name CA', 'resname GLU and resid 18 and name CA',
                  'resname LEU and resid 9 and name CA', 'resname ILE and resid 112 and name CA',
                  'resname PHE and resid 123 and name CA', 'resname TRP and resid 7 and name CA',
                  'resname ALA and resid 134 and name CA', 'resname LEU and resid 135 and name CA',
                  'resname LYS and resid 77 and name CA', 'resname HSD and resid 24 and name CA']

    pocket_XE7 = ['resname ILE and resid 111 and name CA', 'resname VAL and resid 17 and name CA',
                  'resname LEU and resid 69 and name CA', 'resname LEU and resid 115 and name CA',
                  'resname TRP and resid 14 and name CA', 'resname HSD and resid 24 and name CA',
                  'resname ILE and resid 112 and name CA', 'resname VAL and resid 114 and name CA',
                  'resname GLY and resid 25 and name CA', 'resname LEU and resid 72 and name CA',
                  'resname ILE and resid 28 and name CA', 'resname GLY and resid 73 and name CA',
                  'resname VAL and resid 21 and name CA', 'resname ALA and resid 110 and name CA',
                  'resname SER and resid 108 and name CA', 'resname VAL and resid 13 and name CA',
                  'resname GLU and resid 18 and name CA', 'resname LYS and resid 16 and name CA',
                  'resname MET and resid 131 and name CA', 'resname HSD and resid 113 and name CA',
                  'resname VAL and resid 68 and name CA', 'resname ILE and resid 107 and name CA',
                  'resname THR and resid 70 and name CA', 'resname ALA and resid 15 and name CA',
                  'resname VAL and resid 10 and name CA', 'resname LEU and resid 135 and name CA',
                  'resname ASP and resid 27 and name CA']

    pocket_XE8 = ['resname GLY and resid 73 and name CA', 'resname LEU and resid 11 and name CA',
                  'resname LEU and resid 76 and name CA', 'resname TRP and resid 14 and name CA',
                  'resname LYS and resid 77 and name CA', 'resname VAL and resid 10 and name CA',
                  'resname ALA and resid 74 and name CA', 'resname LEU and resid 72 and name CA',
                  'resname VAL and resid 13 and name CA', 'resname ILE and resid 75 and name CA',
                  'resname ALA and resid 15 and name CA', 'resname TRP and resid 7 and name CA',
                  'resname HSD and resid 12 and name CA', 'resname VAL and resid 17 and name CA',
                  'resname MET and resid 131 and name CA', 'resname ALA and resid 134 and name CA',
                  'resname GLN and resid 8 and name CA', 'resname LEU and resid 69 and name CA',
                  'resname LYS and resid 78 and name CA', 'resname LEU and resid 9 and name CA',
                  'resname LYS and resid 79 and name CA', 'resname GLU and resid 18 and name CA',
                  'resname THR and resid 70 and name CA', 'resname ALA and resid 71 and name CA']

    pocket_XE9 = ['resname TRP and resid 14 and name CA', 'resname VAL and resid 13 and name CA',
                  'resname VAL and resid 10 and name CA', 'resname MET and resid 131 and name CA',
                  'resname LEU and resid 115 and name CA', 'resname PHE and resid 123 and name CA',
                  'resname LEU and resid 11 and name CA', 'resname VAL and resid 17 and name CA',
                  'resname GLN and resid 128 and name CA', 'resname ILE and resid 112 and name CA',
                  'resname ILE and resid 111 and name CA', 'resname HSD and resid 12 and name CA',
                  'resname ALA and resid 127 and name CA', 'resname ALA and resid 15 and name CA',
                  'resname LYS and resid 16 and name CA', 'resname LEU and resid 9 and name CA',
                  'resname ASN and resid 132 and name CA', 'resname GLY and resid 73 and name CA',
                  'resname ASP and resid 122 and name CA', 'resname ALA and resid 130 and name CA',
                  'resname ALA and resid 134 and name CA', 'resname LEU and resid 76 and name CA',
                  'resname LEU and resid 135 and name CA', 'resname HSD and resid 119 and name CA',
                  'resname VAL and resid 114 and name CA', 'resname HSD and resid 116 and name CA',
                  'resname LEU and resid 69 and name CA', 'resname LEU and resid 72 and name CA',
                  'resname GLU and resid 18 and name CA']

    # initialize arrays which will hold the distance from H2 to the corresponding pocket
    pocket_arr = np.zeros(num_timesteps)
    pocket_XE1_dist_arr = np.zeros(num_timesteps)
    pocket_XE2_dist_arr = np.zeros(num_timesteps)
    pocket_XE3_dist_arr = np.zeros(num_timesteps)
    pocket_XE4_dist_arr = np.zeros(num_timesteps)
    pocket_XE5_dist_arr = np.zeros(num_timesteps)
    pocket_XE6_dist_arr = np.zeros(num_timesteps)
    pocket_XE7_dist_arr = np.zeros(num_timesteps)
    pocket_XE8_dist_arr = np.zeros(num_timesteps)
    pocket_XE9_dist_arr = np.zeros(num_timesteps)
    pocket_XE10_dist_arr = np.zeros(num_timesteps)

    count = -1

    # Open dcd file
    dcd_list = []
    for i in range(num_files):
        dcd_list.append(MDAnalysis.Universe(path_psf, paths[i]))

    # Get atom indices for pocket
    XE1_pocket_indices = read_pocket_list(pocket_XE1)
    XE2_pocket_indices = read_pocket_list(pocket_XE2)
    XE3_pocket_indices = read_pocket_list(pocket_XE3)
    XE4_pocket_indices = read_pocket_list(pocket_XE4)
    XE5_pocket_indices = read_pocket_list(pocket_XE5)
    XE6_pocket_indices = read_pocket_list(pocket_XE6)
    XE7_pocket_indices = read_pocket_list(pocket_XE7)
    XE8_pocket_indices = read_pocket_list(pocket_XE8)
    XE9_pocket_indices = read_pocket_list(pocket_XE9)

    # Get atom indices for H2
    H2_arr = read_pocket_list(H2_atoms)

    for i in range(num_files):
        print("evaluating file " + str(i+1) + " for H2 " + str(h2_index+1))
        for jj, frame in enumerate(dcd_list[i].trajectory[:num_timesteps]):

            if not count%1000:
                print("frame: "+str(count))

            count += 1

            # Get updated H2 position
            positions_H2 = frame._pos[H2_arr]
            mean_position_H2 = np.mean(positions_H2, axis=0)

            # Get updated position of Pocket
            pocket_XE1_dist_arr[count] = calc_pocket_pos(XE1_pocket_indices)
            pocket_XE2_dist_arr[count] = calc_pocket_pos(XE2_pocket_indices)
            pocket_XE3_dist_arr[count] = calc_pocket_pos(XE3_pocket_indices)
            pocket_XE4_dist_arr[count] = calc_pocket_pos(XE4_pocket_indices)
            pocket_XE5_dist_arr[count] = calc_pocket_pos(XE5_pocket_indices)
            pocket_XE6_dist_arr[count] = calc_pocket_pos(XE6_pocket_indices)
            pocket_XE7_dist_arr[count] = calc_pocket_pos(XE7_pocket_indices)
            pocket_XE8_dist_arr[count] = calc_pocket_pos(XE8_pocket_indices)
            pocket_XE9_dist_arr[count] = calc_pocket_pos(XE9_pocket_indices)

            # set current pocket to the smallest element = the pocket closest to the H2
            dist_list = [pocket_XE1_dist_arr[count], pocket_XE2_dist_arr[count], pocket_XE3_dist_arr[count],
                         pocket_XE4_dist_arr[count], pocket_XE5_dist_arr[count], pocket_XE6_dist_arr[count],
                         pocket_XE7_dist_arr[count], pocket_XE8_dist_arr[count], pocket_XE9_dist_arr[count]]

            dist_list = np.array(dist_list)
            pocket_arr[count] = np.argmin(dist_list) + 1

        np.savez(distfile, pocket_arr=pocket_arr, pocket_XE1=pocket_XE1_dist_arr,
                 pocket_XE2=pocket_XE2_dist_arr, pocket_XE3=pocket_XE3_dist_arr, pocket_XE4=pocket_XE4_dist_arr,
                 pocket_XE5=pocket_XE5_dist_arr, pocket_XE6=pocket_XE6_dist_arr, pocket_XE7=pocket_XE7_dist_arr,
                 pocket_XE8=pocket_XE8_dist_arr, pocket_XE9=pocket_XE9_dist_arr)
