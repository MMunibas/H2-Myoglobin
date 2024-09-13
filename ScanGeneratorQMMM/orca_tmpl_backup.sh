#!/bin/bash
#SBATCH --job-name=orca_scan
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=1000

#--------------------------------------
# Modules
#--------------------------------------

my_orca=/home/toepfer/programs/orca/orca
ulimit -s 10420

#--------------------------------------
# Prepare Run
#--------------------------------------

export SLURMFILE=slurm-$SLURM_JOBID.out

#--------------------------------------
# Run ORCA
#--------------------------------------

$my_orca test.inp > test.out

