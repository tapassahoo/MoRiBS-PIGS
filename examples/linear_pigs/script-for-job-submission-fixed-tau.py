#!/usr/bin/python
# Example PBS cluster job submission in Python
 
import time
from subprocess import call
from os import system
import os
import decimal
 

def dropzeros(number):
    mynum          = decimal.Decimal(number).normalize()
    # e.g 22000 --> Decimal('2.2E+4')
    return mynum.__trunc__() if not mynum % 1 else float(mynum)


def bconstant():
	'''
	This function calculates rotational Bconstant for linear rotor
	'''
	energyj0       = -36117.5942855
	energyj1       = -35999.1009407
	bconst         = 0.5*(energyj1-energyj0)     # in cm^-1
	return bconst


def replace(string_old, string_new, file1, file2):
	'''
	This function replaces old string by new string
	'''
	f1             = open(file1, 'r')
	f2             = open(file2, 'w')
	for line in f1:
		f2.write(line.replace(string_old, string_new))
	f1.close()
	f2.close()


def beads(tau,beta):
	'''
	This function determins number of beads
	'''
	numbbeads1     =beta/tau+1
	numbbeads2     = int(round(numbbeads1,0))
	if (numbbeads2 % 2 == 0):
		numbbeads2 = numbbeads2 + 1
	return numbbeads2

def exact_value(tau,numbbeads):
	'''
	This function determins exact value of beta
	'''
	beta_exact     = tau*(numbbeads - 1)
	return beta_exact

def jobstring(file_name,value):
	'''
	This function creats jobstring for #PBS script
	'''
	job_name       = "job_"+str(file_name)+"%5.3f" % value
	walltime       = "100:00:00"
	processors     = "nodes=1:ppn=12"
	command_pimc_run = "./pimc"

	job_string     = """    #!/bin/bash
	#PBS -N %s
	#PBS -l walltime=%s
	#PBS -l %s
	#PBS -o %s.out
	#PBS -e %s.err
	export OMP_NUM_THREADS=12
	cd $PBS_O_WORKDIR
	%s""" % (job_name, walltime, processors, job_name, job_name, command_pimc_run)
	print job_string
	return job_string

def outputstring2(numbbeads,tau,temperature,beta,temperature_exact,beta_exact):
	'''
    This function creats jobstring for #PBS script
	'''
	argu1          = "%7d" % numbbeads
	argu2          = "%5.3f" % tau
	argu3	       = "%7.5f" % temperature
	argu4	       = "%7.5f" % beta
	argu5	       = "%7.5f" % temperature_exact
	argu6	       = "%7.5f" % beta_exact
	output2 ="numbbeads = "+argu1+", tau = "+argu2+", temperature = "+argu3+", beta = "+argu4+", exact_temperature = "+argu5+", exact_beta = "+argu6+"\n"
	return output2

def modify_input(temperature,numbbeads,numbblocks,distance):
	'''
	This function modifies parameters in qmc1.input
	'''
	replace("temperature_input", str(temperature), "qmc1.input", "qmc2.input")
	replace("numbbeads_input", str(numbbeads), "qmc2.input", "qmc3.input")
	replace("numbblocks_input", str(numbblocks), "qmc3.input", "qmc4.input")
	replace("distance_input", str(distance), "qmc4.input", "qmc.input")
	call(["rm", "qmc2.input", "qmc3.input", "qmc4.input"])


def rotmat(molecule,temperature,numbbeads):
	'''
	This function generates rotational density matrics .rot
	'''
	#temperature1    = dropzeros(temperature)
	temperature1    = "%5.3f" % temperature
	command_linden_run = "../../../linear_prop/linden.x "+str(temperature)+" "+str(numbbeads)+" "+str(bconstant())+" 1500 -1"
	system(command_linden_run)
	file_rotdens    = molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	call(["cp", "linden.out", file_rotdens])


#generating linden.out
path_enter_linden = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/linear_prop/"
os.chdir(path_enter_linden)
call(["make", "clean"])
call(["make"])
path_exit_linden = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"
os.chdir(path_exit_linden)

#initial parameters for qmc.input
call(["rm", "yw*"])
molecule         = "H2"
numbblocks	     = 10000
temperature	     = 10.0
tau 		     = 0.0005

ntemp   	     = 10
tempmin 	     = 0.0
tempmax 	     = 2.5
dtemp   	     = (tempmax - tempmin)/ntemp

ntau    	     = 20
dtau    	     = 0.01

rmin             = 3.0
rmax             = 10.0
nr               = 70
dr               = (rmax-rmin)/nr

dbeta            = 0.001
nbeta            = 100

src_path         = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"
nrange           = nbeta                                                                              #change
displacement     = dbeta                                                                              #change
file1_name       = "tau"+str(tau)+"blocks"+str(numbblocks)+"e0vsbeta"                                 #change
file2_name       = "K-1"                                                                              #change
argument2        = "beta"                                                                             #change
value_min        = 0.0                                                                                #change

fw = open("Input-parameters-beta-2-molecule-tau"+str(tau)+"-blocks"+str(numbblocks)+".txt", "a")      #change
# Loop over your jobs
for i in range(1, nrange+1): 
 
	value = value_min + i*displacement          
	r_dis        = "%5.3f" % value

#	tau          = value                                    #change
	beta         = value #1.0/temperature    #value         #change
	temperature  = 1.0/beta                                 #change

	numbbeads    = beads(tau,beta)                         
	beta_exact   = exact_value(tau,numbbeads)
	temperature_exact  = 1.0/beta_exact

	fldr         = file1_name+r_dis+file2_name
	folder_run   = fldr
	call(["mkdir", folder_run])
	call(["mkdir", "-p", folder_run+"/results"])

	# copy files to running folder
	dest_path    = src_path +folder_run
	source_file  = src_path + "h2n2ogr.pot"
	call(["cp", source_file, dest_path])
	source_file  = src_path + "xyz.init"
	call(["cp", source_file, dest_path])
	source_file  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"
	call(["cp", source_file, dest_path])

	argument1    = value        
	modify_input(temperature,numbbeads,numbblocks,argument1)
	source_file = src_path + "qmc.input"
	call(["cp", source_file, dest_path])
     
	# Write submit file for the current cycle
	os.chdir(dest_path)
	rotmat(molecule,temperature,numbbeads)

    #job submission
	fname        = 'submit_'+str(i)
	fwrite       = open(fname, 'w')

	fwrite.write(jobstring(argument2, value))

	fwrite.close()
	call(["qsub", fname])
	os.chdir(src_path)

	fw.write(outputstring2(numbbeads,tau,temperature,beta,temperature_exact,beta_exact))

fw.close()
