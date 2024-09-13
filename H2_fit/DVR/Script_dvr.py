#!/usr/bin/python

# Basics
import numpy as np

# Matplotlib
import matplotlib
import matplotlib.pyplot as plt

# DVR
from ase import units
from fourier_DVR_1D import Domain_Fourier_DVR_1D

#---------------------------------
# Input Section
#---------------------------------

# Hydrogen Mass
mH = 1.008

# Harmonic potential
kh2 = 350.000*units.kcal/units.mol*units.Bohr**2/units.Hartree
rh2 = 0.7414/units.Bohr
def Vharm(x):
    return kh2*(x - rh2)**2 

# Morse potential
deh2 = 111.763853757505*units.kcal/units.mol/units.Hartree
reh2 = 7.47755368e-01/units.Bohr
beta = 1.94869413e+00*units.Bohr
def Vmorse(x):
    return deh2*(1 - np.exp(-beta*(x - reh2)))**2

V = Vmorse

#---------------------------------
# DVR Section
#---------------------------------

# Cinversion units
A2Bohr = 1./units.Bohr
kcalmol2Hartree = units.kcal/units.mol/units.Hartree
Hartree2rcm = 219474.63
Debye2eBohr = 1./2.541746473
Debye2eA = 0.2081943

# DVR parameter
m = mH*mH/(mH + mH)*units._amu/units._me
x_min = 0.2*A2Bohr
x_max = 2.2*A2Bohr
n_DVR = 201
n_g = 2001
n_states = 15
n_plot = 5
scale = 0.003

# Solve
domain = Domain_Fourier_DVR_1D(x_min, x_max, n_DVR)
E, E_four = domain.solve(m, V, n_states=n_states)

# Evaluate eigenstates on grid
x = np.linspace(x_min, x_max, n_g)
dx = x[1] - x[0]
psi_x = domain.grid(x, E_four[:,:n_plot])

# Print energies
print('Eigenenergies')
print('-------------')
print
for i, e in enumerate(E[:n_plot]):
    if i:
        print(
            '%3i: %12.6f Ha,  %4.6f cm-1,   %4.6f cm-1' % (
                i, e, e*Hartree2rcm, (e - E[i - 1])*Hartree2rcm))
    else:
        print('%3i: %12.6f Ha,  %4.6f cm-1' % (i, e, e*Hartree2rcm))

# Print Potential
plt.plot(x*units.Bohr, V(x)*units.Hartree)
plt.show()

