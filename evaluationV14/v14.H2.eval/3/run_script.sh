#!/bin/bash
#SBATCH --job-name=script
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=2000

#--------------------------------------
# Modules
#--------------------------------------


#--------------------------------------
# Prepare Run
#--------------------------------------

export SLURMFILE=slurm-$SLURM_JOBID.out

#--------------------------------------
# Run Scripts
#--------------------------------------

#python readV14_H2_1.py
python powerSpectraReadv14.H2_1.py

