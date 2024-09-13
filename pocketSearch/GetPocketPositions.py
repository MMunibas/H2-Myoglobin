import os
import numpy as np
from glob import glob

# MDAnalysis
import MDAnalysis

import matplotlib.pyplot as plt

from ase import Atoms
from ase import io


def calc_pocket_pos(pocket_indices):
    positions_pocket = pdb.trajectory[0]._pos[pocket_indices]
    mean_position_pocket = np.mean(positions_pocket, axis=0)
    return mean_position_pocket


def read_pocket_list(res_list):
    pocket_indices = []
    for res in res_list:
        pocket_indices.append(pdb.select_atoms(res).ix[0])
    return np.array(pocket_indices)



path_psf = '../Data/MyoglobinStructure/step3_pbcsetup.psf'

myo_pdb = "../Data/MyoglobinStructure/step4_equilibration.pdb"


num_timesteps = 10000

H2_atoms = ["resname H2 and resid 156 and name H1", "resname H2 and resid 156 and name H2"]

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


# Open pdb file
pdb = MDAnalysis.Universe(path_psf, myo_pdb)

# Get atom indices for pocket

pocket_list = []

XE1_pocket_indices = read_pocket_list(pocket_XE1)
XE2_pocket_indices = read_pocket_list(pocket_XE2)
XE3_pocket_indices = read_pocket_list(pocket_XE3)
XE4_pocket_indices = read_pocket_list(pocket_XE4)
XE5_pocket_indices = read_pocket_list(pocket_XE5)
XE6_pocket_indices = read_pocket_list(pocket_XE6)
XE7_pocket_indices = read_pocket_list(pocket_XE7)
XE8_pocket_indices = read_pocket_list(pocket_XE8)
XE9_pocket_indices = read_pocket_list(pocket_XE9)

pocket_list.append(XE1_pocket_indices)
pocket_list.append(XE2_pocket_indices)
pocket_list.append(XE3_pocket_indices)
pocket_list.append(XE4_pocket_indices)
pocket_list.append(XE5_pocket_indices)
pocket_list.append(XE6_pocket_indices)
pocket_list.append(XE7_pocket_indices)
pocket_list.append(XE8_pocket_indices)
pocket_list.append(XE9_pocket_indices)


for i in range(9):
    print(calc_pocket_pos(pocket_list[i]))




























