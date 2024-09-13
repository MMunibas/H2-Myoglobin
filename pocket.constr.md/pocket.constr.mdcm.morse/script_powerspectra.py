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
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea

# Miscellaneous
import ase.units as units

# -----------------------------
# Parameters
# -----------------------------

# Data files
datadir = 'data_pocket'
datafile_tag = 'data_h2.{:d}.npz'
datasubfile_tag = 'data_h2.{:d}.{:d}.npz'

# Working directory
workdir_tag = 'pocket_*'
split_workdir = ['_', -1]

# DCD and PSF file tags
dcdfile_tag = 'dyna*.dcd'
split_dcdfile = [['.', 0], ['dyna', 1]]
psffile = 'step3_pbcsetup.psf'

# H2 segid
segid_h2 = 'H2'

# Temperatures (K)
T = 300.0

# Single H2 frequency
gas_freq = 4465.5

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
# Read Data
# -----------------------------

# Check data directory
if not os.path.exists(datadir):
    os.makedirs(datadir)

# Detect pocket working directories
workdirs = glob(workdir_tag)

# Get the pocket directory index numbers
ipockets = np.array([
    int(workdir.split('/')[-1].split(split_workdir[0])[split_workdir[1]])
    for workdir in workdirs])

# Sort pocket directories
sort_pockets = np.argsort(ipockets)
workdirs = np.array(workdirs)[sort_pockets]
ipockets = ipockets[sort_pockets]

# Iterate over pockets
for iw, workdir in enumerate(workdirs):

    # Pocket index
    ipocket = ipockets[iw]

    # Check if results already exists
    if os.path.exists(
        os.path.join(datadir, datafile_tag.format(ipocket))
    ):
        continue
    
    # Detect all dcd files
    dcdfiles = glob(os.path.join(workdir, dcdfile_tag))

    # Get the dcd index numbers
    idcds = []
    for dcdfile in dcdfiles:
        idcd = dcdfile.split('/')[-1]
        for (split_tag, split_id) in split_dcdfile:
            idcd = idcd.split(split_tag)[split_id]
        idcds.append(int(idcd))
    idcds = np.array(idcds)
    
    # Sort dcd files
    sort_dcds = np.argsort(idcds)
    dcdfiles = np.array(dcdfiles)[sort_dcds]
    idcds = idcds[sort_dcds]
    
    # Prepare auxiliary variables
    Nframes_h2 = None
    dt_h2 = None

    # Prepare H2 distances array
    distances_h2 = []

    # Loop over dcd files
    for ii, idcd in enumerate(idcds):
        
        # Check if results already exists
        if os.path.exists(
            os.path.join(datadir, datasubfile_tag.format(ipocket, idcd))
        ):
            
            # Load dcd result data
            dcd_data = np.load(
                os.path.join(datadir, datasubfile_tag.format(ipocket, idcd)))
            distances = dcd_data['distances_h2']
            dt = dcd_data['dt_h2']
            
            # Check time step
            if dt_h2 is None:
                dt_h2 = dt
            elif dt_h2 != dt:
                raise SyntaxError()
        
        else:

            # Open dcd file
            dcdfile = dcdfiles[ii]
            dcd = MDAnalysis.Universe(
                os.path.join(workdir, psffile),
                dcdfile)

            # Get trajectory parameter
            Nframes = len(dcd.trajectory)
            Nskip = int(dcd.trajectory.skip_timestep)
            dt = np.round(
                float(dcd.trajectory._ts_kwargs['dt']), decimals=8)/Nskip

            # Check Nframes
            if Nframes_h2 is None:
                Nframes_h2 = Nframes
            elif Nframes_h2 != Nframes:
                continue
            
            # Check time step
            if dt_h2 is None:
                dt_h2 = dt*Nskip
            elif dt_h2 != dt*Nskip:
                raise SyntaxError()

            # Get H2 atom selection
            dcd_h2 = dcd.select_atoms(f'segid {segid_h2:s}')

            # Prepare H2 positions array
            positions_h2 = np.zeros([Nframes, 2, 3], dtype=float)
            
            # Iterate over frames
            iprint = Nframes//100
            for jj, frame in enumerate(dcd.trajectory):
                
                if not jj%iprint:
                    print(
                        f"Frame {jj:{len(str(Nframes)):d}d} of {Nframes:d} of "
                        + f"{dcdfile:s}.")

                # Get H2 atom positions in frame
                positions_h2[jj, :, :] = frame._pos[dcd_h2._ix]

            # Compute bond distances
            distances = np.sqrt(
                np.sum(
                    (positions_h2[:, 1, :] - positions_h2[:, 0, :])**2,
                    axis=-1)
                )

            np.savez(
                os.path.join(datadir, datasubfile_tag.format(ipocket, idcd)),
                distances_h2=distances,
                dt_h2=dt*Nskip)

        # Append to list
        distances_h2.append(distances)

    # Concatenate H2 distances
    distances_h2 = np.concatenate(distances_h2)
    
    # Store results
    np.savez(
        os.path.join(datadir, datafile_tag.format(ipocket)),
        distances_h2=distances_h2,
        dt_h2=dt_h2)

# -----------------------------
# Predict & Plot Powerspectra
# -----------------------------

dpi = 200

# Fontsize
TINY_SIZE = 9
XSMALL_SIZE = 11
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE, weight='bold')  # controls default text sizes
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=XSMALL_SIZE)    # legend fontsize

# Line colors
color_scheme = [
    'b', 'r', 'g', 'purple', 'orange', 'magenta', 'brown', 'darkblue',
    'darkred', 'darkgreen', 'darkgrey', 'olive']

# Spectra shift
spectra_shift = 0.5

# Plot size
figsize = (6, 8)
left = 0.15
bottom = 0.10
column = [0.75, 0.0]
row = [0.70, 0.0]

# Initialize plot
fig = plt.figure(figsize=figsize)
axs = fig.add_axes([left + 0*np.sum(column), bottom, column[0], row[0]])

def moving_average(data_set, periods=9):
    weights = np.ones(periods) / periods
    return np.convolve(data_set, weights, 'same')

pocket_label = {
    1: "Xe1",
    2: "Xe2",
    3: "Xe3",
    4: "Xe4",
    5: "B-state",
    6: "6",
    7: "7",
    8: "8",
    9: "9",
    }


# Iterate over pockets
for iw, workdir in enumerate(workdirs):

    # Pocket index
    ipocket = ipockets[iw]

    # Read data
    data_h2 = np.load(os.path.join(datadir, datafile_tag.format(ipocket)))
    distances_h2 = data_h2['distances_h2']
    dt_h2 = data_h2['dt_h2']

    # Time array
    times_h2 = np.arange(distances_h2.shape[0], dtype=float)*dt_h2

    # Frequency array
    nf = len(times_h2)
    Nfrq = int(nf / 2) + 1
    freq = np.arange(Nfrq)/float(nf)/dt_h2*jiffy

    # Compute IR spectra
    acv = acovf(distances_h2, fft=True)

    acv = acv*np.blackman(nf)
    spec = np.abs(np.fft.rfftn(acv))

    beta = 1.0/(kB_Ha_K*au2cminv)/T
    spec = spec*freq*(1 - np.exp(-freq*beta))

    # Apply moving average
    favg = 4.0
    Nave = int(favg / (freq[1] - freq[0]))
    if Nave:
        spec = moving_average(spec, Nave)

    # Get H2 max vibration amplitude
    select = np.logical_and(freq > 4200., freq < 4600.)
    ampl_h2 = np.max(spec[select])

    # Compute peak position
    peak = np.sum(freq[select]*spec[select])/np.sum(spec[select])

    # Plot powerspectra
    #label = f"Pocket {ipocket:d}"
    #label = f"{ipocket:d} ({peak:5.1f})"
    #nn = 0
    #for ii in range(3*(iw//3), (3*(iw//3 + 1))):
        #nnii = len(pocket_label[ii + 1])
        #if nnii > nn:
            #nn = nnii
    label = f"{pocket_label[ipocket]:s} ({peak:>5.1f})"
    #label = f"{pocket_label[ipocket]:s}"
    yshift = (len(workdirs) - 1 - iw)*spectra_shift
    axs.plot(
        freq, spec/ampl_h2 + yshift, color=color_scheme[iw], label=label)

    label = f"{peak:5.1f}"
    axs.scatter(
        peak, 0.8 + yshift, color=color_scheme[iw])

# Plot gas frequency
ymax = len(workdirs)*spectra_shift + spectra_shift
axs.plot(
    [gas_freq, gas_freq], [0., ymax],
    '--k')

# Plot options
axs.set_xlim(4425., 4550)
axs.set_ylim(0., ymax)

axs.set_yticks(np.arange(0, ymax + spectra_shift/2., spectra_shift))
axs.set_yticklabels([])

axs.set_xlabel(r'Frequency (cm$^{-1}$)', fontweight='bold')
axs.get_xaxis().set_label_coords(0.50, -0.08)

axs.set_ylabel(r'Intensity (arb. units)', fontweight='bold')
axs.get_yaxis().set_label_coords(-0.10, 0.50)

axs.legend(
    loc=(0.1, 1.02), framealpha=1.0, ncol=2,
    title=r'H$_2$ in Pocket (Mean Frequency in cm$^{-1}$)')

tbox = TextArea(
    'C', 
    textprops=dict(
        color='k', fontsize=35, ha='center', va='center')
    )
anchored_tbox = AnchoredOffsetbox(
    loc="upper right", child=tbox, pad=0., frameon=False,
    bbox_to_anchor=(0.12, 0.97),
    bbox_transform=axs.transAxes, borderpad=0.)
axs.add_artist(anchored_tbox)


plt.savefig(
    os.path.join(datadir, 'h2_spectra_pockets_mdcm_morse.png'), 
    format='png', dpi=dpi)
plt.close()   





# TOC graphics


figsize = (4, 4)
left = 0.15
bottom = 0.10
column = [0.70, 0.00]
row = [0.80, 0.00]

fig = plt.figure(num=0, figsize=figsize)

axs = fig.add_axes([left, bottom, column[0], row[0]], projection='3d')

axs.view_init(elev=15., azim=-75.)

plt.axis('off')

pocket_label = ['Pocket Xe4', 'Pocket 7']


# Iterate over pockets
for ii, iw in enumerate([3, 6]):

    # Pocket workdir
    workdir = workdirs[iw]

    # Pocket index
    ipocket = ipockets[iw]

    # Read data
    data_h2 = np.load(os.path.join(datadir, datafile_tag.format(ipocket)))
    distances_h2 = data_h2['distances_h2']
    dt_h2 = data_h2['dt_h2']

    # Time array
    times_h2 = np.arange(distances_h2.shape[0], dtype=float)*dt_h2

    # Frequency array
    nf = len(times_h2)
    Nfrq = int(nf / 2) + 1
    freq = np.arange(Nfrq)/float(nf)/dt_h2*jiffy

    # Compute IR spectra
    acv = acovf(distances_h2, fft=True)

    acv = acv*np.blackman(nf)
    spec = np.abs(np.fft.rfftn(acv))

    beta = 1.0/(kB_Ha_K*au2cminv)/T
    spec = spec*freq*(1 - np.exp(-freq*beta))

    # Apply moving average
    favg = 4.0
    Nave = int(favg / (freq[1] - freq[0]))
    if Nave:
        spec = moving_average(spec, Nave)

    # Get H2 max vibration amplitude
    select = np.logical_and(freq > 4420., freq < 4500.)
    ampl_h2 = np.max(spec[select])

    # Compute peak position
    peak = np.sum(freq[select]*spec[select])/np.sum(spec[select])

    axs.plot(
        freq[select],
        [0.5*iw]*len(freq[select]),
        spec[select]/ampl_h2 + iw*0.05,
        '-',
        color=color_scheme[iw],
        lw=2,
        zorder=ii+1)

    axs.text(
        freq[select][-1] + 5,
        0.5*iw + (ii - 2)*0.5,
        spec[select][-1]/ampl_h2 + iw*0.05,
        pocket_label[ii])


data = np.load("md_h2_morse_spec.npz")
freq = data['freq']
spec = data['spec']
select = np.logical_and(freq > 4420., freq < 4500.)
ampl_h2 = np.max(spec[select])
iw += 4
axs.plot(
    freq[select],
    [0.5*iw]*len(freq[select]),
    spec[select]/ampl_h2 + iw*0.05,
    '--',
    color='k',
    lw=2,
    zorder=0)

axs.text(
    freq[select][-1] + 5,
    0.5*iw + (ii - 2)*0.5,
    spec[select][-1]/ampl_h2 + iw*0.05,
    'Gas')


axs.text(
    freq[select][-1],
    -10.0,
    spec[select][0]/ampl_h2,
    #r'FT($\langle r_\mathrm{H2}(t) \cdot r_\mathrm{H2}(0) \rangle$)'
    r'$\omega$')

#axs.set_xlim(4450, 4500)
axs.set_ylim(0, 10)
axs.set_zlim(0, 2)

plt.savefig(
    'toc_spec.png',
    format='png', dpi=600)
plt.close()
               
