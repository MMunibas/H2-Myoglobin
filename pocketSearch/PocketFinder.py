import os
from typing import List, Any

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


split_dcd = ['.', 1]  # split a each '.' and take element [1] in filename_dcd

path_psf = '../Data/v14/step3_pbcsetup.psf'

datadir = '/Data/v14/dyna11.dcd'

file_list_dcd = [datadir]

# for i in range(1, 11):
#     file_list_dcd.append('../Data/dcdFiles/dyna.' + str(i) + '.dcd')

start_frame = 0
end_frame = 250

idcds = np.array([
    int(file_dcd.split('/')[-1].split(split_dcd[0])[split_dcd[1]])
    for file_dcd in file_list_dcd])

sort_dcd = np.argsort(idcds)
file_list_dcd = np.array(file_list_dcd)[sort_dcd]
idcds = idcds[sort_dcd]

# initialize arrays which will hold the distance from H2 to the corresponding pocket

pocket_storage_pseudodict = []

count = -1

# Loop over dcd files
for ii, idcd in enumerate(file_list_dcd):

    # Open dcd file
    file_dcd = idcd
    dcd = MDAnalysis.Universe(path_psf, file_dcd)

    # XE4c6 = dcd.select_atoms("resname LEU and resid 72 and name C").ix[0]

    # Get atom indices for H2
    H1 = dcd.select_atoms("resname H2 and resid 156 and name H1").ix[0]
    H2 = dcd.select_atoms("resname H2 and resid 156 and name H2").ix[0]
    H2_arr = np.array([H1, H2])

    # Get trajectory parameter
    # Nframes: Number of frames in dcd trajectory file
    # Nskip: equivalent to 'nsavc' in the CHARMM input file
    # dt: equivalent to 'timestep' in the CHARMM input file
    Nframes = len(dcd.trajectory)
    Nskip = int(dcd.trajectory.skip_timestep)
    dt = np.round(
        float(dcd.trajectory._ts_kwargs['dt']), decimals=8) / Nskip

    # Get atom types
    atoms = np.array([ai for ai in dcd._topology.names.values])
    Natoms = len(atoms)

    # Get atom masses
    masses = dcd._topology.masses.values

    # Get atom charges
    charges = dcd._topology.charges.values

    # Iterate over frames
    for jj, frame in enumerate(dcd.trajectory):
        count += 1

        if count > start_frame and count < end_frame + 2:

            # Get updated H2 position
            positions_H2 = frame._pos[H2_arr]
            mean_position_H2 = np.mean(positions_H2, axis=0)

            # Calculate distances between residues of molecule "PROA"
            distances = []
            for residue in dcd.residues:
                # Check if the residue belongs to molecule "PROA"
                if residue.segid == "PROA":
                    # Get the alpha carbon atom of the residue
                    alpha_carbon = residue.atoms.select_atoms('name CA')
                    if len(alpha_carbon) == 0:
                        continue  # Skip if alpha carbon is not found
                    alpha_carbon = alpha_carbon.positions[0]  # Get the position of the first alpha carbon
                    dist = np.linalg.norm(alpha_carbon - mean_position_H2)
                    if dist < 10:
                        distances.append((residue, dist))
            sorted_residues = sorted(distances, key=lambda x: x[1])

            res_list = []
            for selection in sorted_residues:
                for atom in selection[0].atoms:
                    if atom.name == "CA":
                        res = f"resname {atom.resname} and resid {atom.resid} and name {atom.name}"
                        res_list.append(res)

            aminolist = []
            for aminoacid in res_list:
                aminolist.append(dcd.select_atoms(aminoacid).ix[0])
            aminotupel = (calc_pocket_pos(aminolist), res_list)

            pocket_storage_pseudodict.append(aminotupel)


            if count == end_frame:
                pocket_storage_pseudodict = sorted(pocket_storage_pseudodict, key=lambda x: x[0])
                for x in range(0,10):
                    print(pocket_storage_pseudodict[x])
                break

        if count % 50 == 0 and count < end_frame:
            print(count)