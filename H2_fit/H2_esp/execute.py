#!/usr/bin/python

# Basics
import os
import numpy as np

# Load FDCM modules
from pydcm import Scan, DCM, dcm_fortran

# MDCM:
#-------

file_mdcm_clcl = "source/charges6_h2_sym.dcm"

# Prepare some cube file list
scan_fesp = ["source/h2_esp.cube"]
scan_fdns = ["source/h2_dens.cube"]

weight_list = [1.]

Nfiles = len(scan_fesp)
Nchars = int(np.max([
    len(filename) for filelist in [scan_fesp, scan_fdns] 
    for filename in filelist]))

esplist = np.empty([Nfiles, Nchars], dtype='c')
dnslist = np.empty([Nfiles, Nchars], dtype='c')

for ifle in range(Nfiles):
    esplist[ifle] = "{0:{1}s}".format(scan_fesp[ifle], Nchars)
    dnslist[ifle] = "{0:{1}s}".format(scan_fdns[ifle], Nchars)

# Load cube files, read MDCM global and local files
dcm_fortran.load_cube_files(Nfiles, Nchars, esplist, dnslist)

dcm_fortran.load_clcl_file(file_mdcm_clcl)

# Write MDCM global from local and Fitted ESP cube files
dcm_fortran.write_cxyz_files()
dcm_fortran.write_mdcm_cube_files()

# Get current local MDCM array
clcl = dcm_fortran.mdcm_clcl

# Constraint charges
fix_q = []
for i in range(3, len(clcl), 4):
    fix_q.append(clcl[i])
# or not constraint, fix_q = None
fix_q = None

clcl_red = []
for i in range(2, len(clcl)//2, 4):
    clcl_red.append(clcl[i])
    clcl_red.append(clcl[i+1])

def mdcm_rmse(
    clcl_red, lmax=[-0.6, +0.6], lfor=1.0e4,
    qtot=0.0,
    qmax=[-1.0, +1.0], qfor=1.0e4,
    qfix=fix_q):

    for ii, ired in enumerate(range(2, len(clcl)//2, 4)):
        clcl[ired] = clcl_red[2*ii]
        clcl[len(clcl)//2 + ired] = -clcl_red[2*ii]
        clcl[ired+1] = clcl_red[2*ii+1]
        clcl[len(clcl)//2 + ired+1] = clcl_red[2*ii+1]

    # Constraint total charge
    if qfix is None:
        qsum = 0.0
        for i in range(3, len(clcl), 4):
            qsum += clcl[i]
        diff = qtot - qsum
        clcl[i//2] += diff/2.
        clcl[i] += diff/2.
        Nchg = float(len(clcl)//4)
        qsum = 0.0
        for i in range(3, len(clcl), 4):
            qsum += clcl[i]
    else:
        qsum = 0.0
        for iq, qi in enumerate(range(3, len(clcl), 4)):
            clcl[qi] = qfix[iq]
            qsum += qfix[iq]
    
    dcm_fortran.set_clcl(clcl)
    rmse = dcm_fortran.get_rmse()
    wrmse = dcm_fortran.get_rmse_weighted(Nfiles, list(weight_list))
    srmse = dcm_fortran.get_rmse_each(Nfiles)
    
    # Check maximum amplitude constraint 
    for i in range(0, len(clcl), 4):
        for j in range(2, 3):
            if clcl[i + j] < lmax[0]:
                rmse += lfor*(lmax[0] - clcl[i + j])**2
            elif clcl[i + j] > lmax[1]:
                rmse += lfor*(clcl[i + j] - lmax[1])**2
    
    # Check maximum charge constraint 
    for i in range(3, len(clcl), 4):
        if clcl[i] < qmax[0]:
            rmse += qfor*(qmax[0] - clcl[i])**2
        elif clcl[i] > qmax[1]:
            rmse += qfor*(clcl[i] - qmax[1])**2
    
    out = f"Total charge: {qsum:.6f}\n"
    out += f" weighted RMSE: {wrmse:.6f} , pure ESP RMSE: {rmse:.6f}\n"
    a02A = 0.52917721067
    pp = "Local charge positions and magnitude:\n"
    for i in range(0, len(clcl), 4):
        pp += (
            f"{clcl[i + 0]*a02A: 11.10f} " + 
            f"{clcl[i + 1]*a02A: 11.10f} " +
            f"{clcl[i + 2]*a02A: 11.10f} " + 
            f"{clcl[i + 3]: 11.10f}\n"
        )
    out += pp
    out += "\n"
    print(out)
    with open('out_mdcm.out', 'a+') as f:
        f.write(out)

    return rmse

# Apply simple minimization without any feasibility check (!)
# Leads to high amplitudes of MDCM charges and local positions
from scipy.optimize import minimize
res = minimize(mdcm_rmse, clcl_red, tol=1e-5)
print(res)

# Create mode cube files
dcm_fortran.write_cxyz_files()
dcm_fortran.write_mdcm_cube_files()

# Not necessary but who knows when it become important to deallocate all 
# global arrays
dcm_fortran.dealloc_all()
