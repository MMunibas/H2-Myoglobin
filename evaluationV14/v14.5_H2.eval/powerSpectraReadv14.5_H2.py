# Basics
import os
import sys

import MDAnalysis
import numpy as np
from glob import glob

# Statistics
from statsmodels.tsa.stattools import acovf

# Matplotlib
import matplotlib
import matplotlib.pyplot as plt

# Miscellaneous
import ase.units as units


def read_pocket_list(res_list, dcd_file):
    pocket_indices = []
    for res in res_list:
        pocket_indices.append(dcd_file.select_atoms(res).ix[0])
    return np.array(pocket_indices)


def moving_average(data_set, periods=9):
    weights = np.ones(periods) / periods
    return np.convolve(data_set, weights, 'same')


# -----------------------------
# Parameters
# -----------------------------

# Temperatures (K)
T = 300.0

# Time for speed of light in vacuum to travel 1 cm in ps
jiffy = 0.01 / units._c * 1e12
# 33.3564 [cm/ps]

# Boltzman constant in Hartree/Kelvin
kB_Ha_K = units.kB / units.Hartree
# 3.16681e-06 [Hartree/K]

# Conversion factor from Hartree (au) to inverse centimeter cm**-1
# 1 Hartree = [Hartree to Joule) / (planck constant [Js] * speed of light [m/s]) * (1/m to 1/cm)] cm**-1
au2cminv = 1.0 * units.Hartree / units.J / (units._hplanck * units._c) * 1.e-2
# 219474.63 [1/Hartree/cm]

# -----------------------------
# Get data
# -----------------------------
version = "14.5_H2"
preprefix = "/cluster/data/toepfer/v" + version + ".nonBonded.morse.mdcm/"


path_psf = preprefix + 'step3_pbcsetup.psf'

file_dcd_prefix = preprefix + 'dyna'
file_dcd_suffix = ".dcd"

timedt = 0.001  # Time / frame in ps

# ... Get H2 bond distances in Angstrom

h2_distances = []

start_dcd = 1
end_dcd = 37
num_files = end_dcd - start_dcd + 1
paths = []

H2_res_list = [156, 157, 158, 159, 160]

for i in range(start_dcd, end_dcd + 1):
    paths.append(file_dcd_prefix + str(i) + file_dcd_suffix)

# init dcd files
dcd_list = []
for i in range(num_files):
    if os.path.exists(paths[i]):
        dcd_list.append(MDAnalysis.Universe(path_psf, paths[i]))
# dcd = MDAnalysis.Universe(path_psf, file_dcd, file_dcd2)
num_files = len(dcd_list)

for h2_index in range(len(H2_res_list)):
    print("starting H" + str(h2_index+1))

    distfile = "h2distsV" + version + ".H2-" + str(h2_index + 1) + ".dcd." + str(start_dcd) + "-" + str(end_dcd) + ".npy"
    # !!!!!!!!!!!!!!! overwrite dist file to read data !!!!!!!!!!!!!!!!!!!!!
    # distfile = ""

    # Get atom indices for H2
    H2_atoms = ["resname H2 and resid " + str(H2_res_list[h2_index]) + " and name H1",
                "resname H2 and resid " + str(H2_res_list[h2_index]) + " and name H2"]
    H2_arr = read_pocket_list(H2_atoms, dcd_list[0])

    if not os.path.exists(distfile):

        num_timesteps = 0
        h2_distances = []

        for i in range(num_files):
            for jj, frame in enumerate(dcd_list[i].trajectory):

                if not jj % 1000:
                    print("frame: " + str(jj))
                num_timesteps += 1
                # Get updated H2 position
                positions_H2 = frame._pos[H2_arr]
                h2_dist = np.linalg.norm((positions_H2[0] - positions_H2[1]))
                h2_distances.append(h2_dist)

        # ... Get time sequence [0.00025, 0.0005, 0.00075, ...] in ps
        time = np.arange(0, num_timesteps) * timedt
        h2_distances = np.array(h2_distances)

        np.save(distfile, h2_distances)

    else:
        h2_distances = np.load(distfile)
        num_timesteps = len(h2_distances)
        time = np.arange(0, num_timesteps) * timedt

    # -----------------------------
    # Predict IR spectra
    # -----------------------------

    # Frequency array
    nf = len(time)
    dt = time[1] - time[0]
    Nfrq = int(nf / 2) + 1
    freq = np.arange(Nfrq) / float(nf) / dt * jiffy

    # Compute IR spectra
    acv = acovf(h2_distances, fft=True)

    acv = acv * np.blackman(nf)
    spec = np.abs(np.fft.rfftn(acv))

    beta = 1.0 / (kB_Ha_K * au2cminv) / T
    spec = spec * freq * (1 - np.exp(-freq * beta))

    # -----------------------------
    # Plot IR spectra
    # -----------------------------

    # Apply moving average

    favg = 1.0
    Nave = int(favg / (freq[1] - freq[0]))
    #spec = moving_average(spec, Nave)

    # Get CN vibration amplitude
    ampl_h2 = np.max(spec[np.logical_and(freq > 3500., freq < 4500.)])

    plt.plot(freq, spec / ampl_h2)

    plt.xlim(500., 4500)

    plt.title("h2spec" + version + "-" + str(h2_index + 1))

    plt.savefig("h2spec" + version + "-" + str(h2_index + 1) + ".dcd." + str(start_dcd) + "-" + str(end_dcd) + ".png",
                format='png')

    specfile = "h2spec" + version + "-" + str(h2_index + 1) + ".dcd." + str(start_dcd) + "-" + str(end_dcd) + ".npz"

    np.savez(specfile, freq=freq, spec=spec)
