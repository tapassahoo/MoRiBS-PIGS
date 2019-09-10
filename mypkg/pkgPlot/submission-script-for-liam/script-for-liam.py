#!/usr/bin/python
import os
from os import system
from subprocess import call

def jobstring_sbatch(T,R,S):
	'''
	This function creats jobstring for #SBATCH script
	'''
	job_name       = "MyjobT"+str(T)+"R"+str(R)+"S"+str(S)
	walltime       = "40-00:00"
	omp_thread     = 1
	logPath        = "/home/tapas/TEST/"

	exe_file       = "/home/l2mcgrat/MMTK/bin/python PolSimStart.py H2O_T"+str(T)+"t1.rho  H2O_T"+str(T)+"t1.eng  H2O_T"+str(T)+"t1.esq "+str(S)+" "+str(R)+" 1"

	job_string     = """#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --output=%s.out
#SBATCH --time=%s
#SBATCH --mem-per-cpu=1200mb
#SBATCH --cpus-per-task=%s
export OMP_NUM_THREADS=%s
%s
""" % (job_name, job_name, walltime, omp_thread, omp_thread, exe_file)

	return job_string

# Main function
RList = [3.0, 4.0]

TList = [15, 70, 95, 293]
SList = [0.05, 0.15, 0.15, 0.18, 0.1, 0.2, 0.2, 0.18]

print(TList)
print(RList)
print(SList)

it = 0
for R in RList:
	for T in TList:
		S = SList[it]
		fname="job-T"+str(T)+"R"+str(R)+"S"+str(it)
		print(fname)
		it += 1
		fileName = open(fname,'w')
		fileName.write(jobstring_sbatch(T,R,S))
		fileName.close()
		call(["sbatch", fname])
		call(["rm", fname])

