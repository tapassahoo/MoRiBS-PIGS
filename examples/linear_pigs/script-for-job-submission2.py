#!/usr/bin/python
# Example PBS cluster job submission in Python
 
import time
from subprocess import call
from os import system
import os
 
#function defined for the calculation of rotational Bconst for linear rotor
def bconstant():
    energyj0      = -36117.5942855
    energyj1      = -35999.1009407
    bconst        = 0.5*(energyj1-energyj0)     # in cm^-1
    return bconst

def replace(string_old, string_new, file1, file2):
	f1 = open(file1, 'r')
	f2 = open(file2, 'w')
	for line in f1:
		f2.write(line.replace(string_old, string_new))
	f1.close()
	f2.close()

def beads(tau,temperature):
	beta=1.0/temperature
	numbbeads1=beta/tau+1
	numbbeads2 = int(round(numbbeads1,0))
	if (numbbeads2 % 2 == 0):
		numbbeads2 = numbbeads2 + 1
	return numbbeads2

def jobstring(temperature):
	job_name = "job_Temp%3.2f" % temperature
	walltime = "100:00:00"
	processors = "nodes=1:ppn=12"
	command_pimc_run = "../../../pimc"

	job_string = """    #!/bin/bash
	#PBS -N %s
	#PBS -l walltime=%s
	#PBS -l %s
	#PBS -o %s.out
	#PBS -e %s.err
	cd $PBS_O_WORKDIR
	%s""" % (job_name, walltime, processors, job_name, job_name, command_pimc_run)
	print job_string
	return job_string

#generating linden.out
path_enter_linden  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/linear_prop/"
os.chdir(path_enter_linden)
call(["make", "clean"])
call(["make"])
path_exit_linden  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"
os.chdir(path_exit_linden)

#initial parameters for qmc.input
call(["rm", "yw*"])
molecule 	= "H2"
numbbeads 	= 129
numbblocks	= 2000
ntemp   	= 10
tempmin 	= 0.0
tempmax 	= 2.5
dtemp   	= (tempmax - tempmin)/ntemp
ntau    	= 10
dtau    	= 0.001
temperature	= 0.5
tau 		= 0.005

src_path  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"
# Loop over your jobs
for i in range(1, ntemp+1):
 
	temperature = i*dtemp	
#	tau=i*dtau
	numbbeads = beads(tau,temperature)
	beta=1.0/temperature
	fldr = "e0vsbeta"+str(beta)+"beta"+str(numbbeads)+"beads"+str(numbblocks)+"blocks"

	# create folder for run pimc 
#	folder_run = fldr+"%3.2fK" % temperature
	folder_run = fldr
	call(["mkdir", folder_run])
	call(["mkdir", "-p", folder_run+"/results"])

	# copy files to running folder
	dest_path = src_path +folder_run
	source_file  = src_path + "h2n2ogr.pot"
	call(["cp", source_file, dest_path])
	source_file  = src_path + "xyz.init"
	call(["cp", source_file, dest_path])

	# for qmc.input
	file_pigs = "pigs%3.2fK" % temperature
	replace("filename_input", str(file_pigs), "qmc1.input", "qmc2.input")
	replace("temperature_input", str(temperature), "qmc2.input", "qmc3.input")
	replace("numbbeads_input", str(numbbeads), "qmc3.input", "qmc4.input")
	replace("numbblocks_input", str(numbblocks), "qmc4.input", "qmc.input")
	call(["rm", "qmc2.input", "qmc3.input", "qmc4.input"])

	source_file  = src_path + "qmc.input"
	call(["cp", source_file, dest_path])
     
	# Write submit file for the current cycle
	path_enter = src_path+folder_run
	os.chdir(dest_path)

	#generating rotational density matrics .rot
	param2  = "%3.2f" % temperature
	command_linden_run = "../../../linear_prop/linden.x "+str(param2)+" "+str(numbbeads)+" "+str(bconstant())+" 1500 -1"
	system(command_linden_run)

	if (i == 4 or i == 8 ):
		param1 = "%d" % temperature
		file_rotdens = molecule+"_T"+str(param1)+"t"+str(numbbeads)+".rot"
	else:
		file_rotdens = molecule+"_T"+str(temperature)+"t"+str(numbbeads)+".rot"

	call(["cp", "linden.out", file_rotdens])
	#

	fname = 'submit_'+str(i)
	fwrite = open(fname, 'w')
	fwrite.write(jobstring(temperature))
	fwrite.close()
	call(["qsub", fname])
	os.chdir(src_path)
