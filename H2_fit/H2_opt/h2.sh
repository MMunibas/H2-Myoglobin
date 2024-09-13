#!/bin/bash

#SBATCH --job-name=h2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1375
# 
hostname

 
#set -xv
export g16root=/opt/cluster/programs/g16-c.01
export GAUSS_SCRDIR=/scratch/toepfer/h2
mkdir -p /scratch/toepfer/h2
source $g16root/g16/bsd/g16.profile

$g16root/g16/g16 /data/toepfer/Project_Heme/H2_fit/H2_opt/h2.com /scratch/toepfer/h2/h2.out

# don't delete the result file if not able to copy to fileserver 
cp /scratch/toepfer/h2/h2.out /data/toepfer/Project_Heme/H2_fit/H2_opt/h2.out 
status=$?
if [ $status -eq 0 ] 
then 
   rm -rf /scratch/toepfer/h2
else
   host=`/bin/hostname`
   /usr/bin/Mail -v -s "Error at end of batch job" $USER <<EOF

At the end of the batch job the system could not copy the output file
	$host:/scratch/toepfer/h2/h2.out
to
	/data/toepfer/Project_Heme/H2_fit/H2_opt/h2.out
Please copy this file by hand or inform the system manager.

EOF
 
fi
