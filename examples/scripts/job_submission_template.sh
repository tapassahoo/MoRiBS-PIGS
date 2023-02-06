#!/bin/bash
#SBATCH --job-name=jobn2p0
#SBATCH --output=jobn2p0.log
#SBATCH --error=jobn2p0.err
#SBATCH --time=00-12:00
#SBATCH --account=rrg-pnroy
#SBATCH --constraint=broadwell
#SBATCH --mem-per-cpu=1024mb
#SBATCH --cpus-per-task=1
export OMP_NUM_THREADS=1
echo $SLURM_JOB_ID
python execution_moribs_driver_linear_molecule.py 
 
