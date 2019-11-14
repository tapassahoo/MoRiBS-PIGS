#!/bin/bash
#SBATCH --job-name=asym
#SBATCH --output=log.out
#SBATCH --time=0-00:30

#SBATCH --mem-per-cpu=4096mb
#SBATCH --cpus-per-task=1
export OMP_NUM_THREADS=1
#           T    P   iodevn ith0 ithend Arot Brot Crot  maxj
./asymrho.x 0.5 4 0 0 0 27.877 14.512 9.285 10
