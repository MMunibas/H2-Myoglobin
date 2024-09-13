#!/bin/bash
#SBATCH --job-name=v13.p7
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=4000

#--------------------------------------
# Modules
#--------------------------------------

module load gcc/gcc-12.2.0-cmake-3.25.1-openmpi-4.1.4 
my_charmm=/cluster/home/toepfer/programs/c49a2/build/cmake/charmm
ulimit -s 10420

#--------------------------------------
# Prepare Run
#--------------------------------------

export SLURMFILE=slurm-$SLURM_JOBID.out

#--------------------------------------
# Run CHARMM
#--------------------------------------
#mpirun -np 24 $my_charmm -i step1_pdbreader.inp -o step1_pdbreader.out
#mpirun -np 24 $my_charmm -i step2.1_waterbox.inp -o step2.1_waterbox.out
#mpirun -np 24 $my_charmm -i step2.2_ions.inp -o step2.2_ions.out
#mpirun -np 24 $my_charmm -i step2_solvator.inp -o step2_solvator.out
#mpirun -np 24 $my_charmm -i step3_pbcsetup.inp -o step3_pbcsetup.out
#mpirun -np 24 $my_charmm -i step4_equilibration.inp -o step4_equilibration.out
mpirun -np 24 $my_charmm -i step4_preparation.inp -o step4_preparation.out
#mpirun -np 24 $my_charmm -i prod22.inp -o prod22.out
mpirun -np 24 $my_charmm -i prod_pocket.inp -o prod_pocket.out

