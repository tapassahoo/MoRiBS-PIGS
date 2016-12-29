#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *

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

def jobstring(file_name,value):
	'''
	This function creats jobstring for #PBS script
	'''
	job_name       = "job_"+str(file_name)+str(value)
	walltime       = "100:00:00"
	processors     = "nodes=1:ppn=16"
	command_pimc_run = "./pimc"

	job_string     = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
#PBS -l %s
#PBS -o %s.out
#PBS -e %s.err
export OMP_NUM_THREADS=16
cd $PBS_O_WORKDIR
%s""" % (job_name, walltime, processors, job_name, job_name, command_pimc_run)
	print job_string
	return job_string

def outputstring1(numbbeads,tau,temperature):
	'''
	This function gives us the exact values of the agruments
	'''
	argu1          = "%7d"   % numbbeads
	argu2          = "%20.15f" % tau
	argu3          = "%7.5f" % temperature
	output ="numbbeads = "+argu1+", tau = "+argu2+", temperature = "+argu3+"\n"
	return output

def outputstring2(numbbeads,tau,mean_pot,mean_tot,mean_rot,error_pot,error_tot,error_rot):
	'''
	This function gives us the output 
	'''
	argu1          = "%5d"   % numbbeads
	argu2          = "%11.6f" % tau
	argu3          = "%10.5f" % mean_pot
	argu4          = "%10.5f" % mean_tot
	argu5          = "%10.5f" % mean_rot
	argu6          = "%10.5f" % error_pot
	argu7          = "%10.5f" % error_tot
	argu8          = "%10.5f" % error_rot
	output  = " "+argu1+" "+argu2+"   "+argu3+"     "+argu4+"     "+argu5
	output += "        "+argu6+"          "+argu7+"          "+argu8+"\n"
	return output

def formatting(status):
	'''
	This function gives us the output 
	'''
	if status == "analysis":
		output     = "#  Beads   Beta     Avg. Potential   Avg. Total   Avg. rotational  Error of Potential     Error of Total    Error of Rotational \n"
		output    += "#          (1/K)      Energy (K)     Energy (K)      Energy (K)        Energy (K)           Energy (K)          Energy (K) \n"
		output    += "#==============================================================================================================================\n"
		return output

def modify_input(temperature,numbbeads,numbblocks,numbmolecules,distance):
	'''
	This function modifies parameters in qmc_run.input
	'''
	replace("temperature_input", str(temperature), "qmc_run.input", "qmc2.input")
	replace("numbbeads_input", str(numbbeads), "qmc2.input", "qmc3.input")
	replace("numbblocks_input", str(numbblocks), "qmc3.input", "qmc4.input")
	replace("numbmolecules_input", str(numbmolecules), "qmc4.input", "qmc5.input")
	replace("distance_input", str(distance), "qmc5.input", "qmc.input")
	call(["rm", "qmc_run.input", "qmc2.input", "qmc3.input", "qmc4.input", "qmc5.input"])


def rotmat(molecule,temperature,numbbeads):
	'''
	This function generates rotational density matrix - linden.dat
	'''
	#temperature1    = dropzeros(temperature)
	temperature1    = "%5.3f" % temperature
	numbbeads1		= numbbeads - 1
	command_linden_run = "../../../linear_prop/linden.x "+str(temperature)+" "+str(numbbeads1)+" "+str(bconstant())+" 1500 -1"
	print command_linden_run
	system(command_linden_run)
	file_rotdens    = molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	call(["mv", "linden.out", file_rotdens])


#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
molecule            = "H2"                                                         #change param1
numbmolecules       = 2                                                            #change param2
numbblocks	        = 10000                                                        #change param3
tau	                = 0.001                                                        #change param4
Rpt                 = 3.5                                                          #change param5

dbeta               = 0.02                                                         #change param6
#status              = "submission"                                                 #change param7
status              = "analysis"                                                   #change param7

src_path            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"
nrange              = 10								                           #change param8
file1_name          = "blocks"+str(numbblocks)+"molecule"+str(numbmolecules)+str(molecule)+"e0vsbeta" 
																				   #change param9
file2_name          = "K-1"                                                      
argument2           = "beta"                                                      

#===============================================================================
#                                                                              |
#   compilation of linden.f to generate rotational density matrix - linden.out |
#                                                                              |
#===============================================================================
if status == "submission":
	call(["rm", "yw*"])
	path_enter_linden = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/linear_prop/"
	os.chdir(path_enter_linden)
	call(["make", "clean"])
	call(["make"])
	path_exit_linden  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"
	os.chdir(path_exit_linden)
	file_input        = "Input-parameters-vs-beta-"+str(numbmolecules)+"-"+str(molecule)+"-fixed-tau"+str(tau)+"-blocks"+str(numbblocks)+".txt"
	call(["rm", file_input])
	fsubmit           = open(file_input, "a")  


#===============================================================================
#                                                                              |
#   Analysis of output files 												   |
#                                                                              |
#===============================================================================
if status == "analysis":
	file_output       = "Energy-vs-beta-"+str(numbmolecules)+"-"+str(molecule)+"-fixed-tau"+str(tau)+"-blocks"+str(numbblocks)+".txt"  
	call(["rm", file_output])
	fanalyze          = open(file_output, "a")           
	fanalyze.write(formatting(status))


# Loop over jobs
for i in range(1,nrange+1):
	
	value        = i*dbeta          
	beta         = "%5.3f" % value
	temperature  = 1./value
	numbbeads    = beads(tau,value)

	folder_run   = file1_name+str(beta)+file2_name
	dest_path    = src_path +folder_run

	if status == "submission":
		call(["rm", "-rf", folder_run])
		call(["mkdir", folder_run])
		call(["mkdir", "-p", folder_run+"/results"])

		# copy files to running folder
		source_file  = src_path + "xyz.init"
		call(["cp", source_file, dest_path])
		source_file  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"
		call(["cp", source_file, dest_path])

		source_file  = src_path + "qmc_run.input"
		call(["cp", source_file, dest_path])
    
		# Write submit file for the current cycle
		os.chdir(dest_path)
		argument1    = Rpt
		modify_input(temperature,numbbeads,numbblocks,numbmolecules,argument1)
		rotmat(molecule,temperature,numbbeads)

		#job submission
		fname        = 'submit_'+str(i)
		fwrite       = open(fname, 'w')
	
		fwrite.write(jobstring(argument2,numbbeads))

		fwrite.close()
		call(["qsub", fname])
		os.chdir(src_path)

		fsubmit.write(outputstring1(numbbeads,tau,temperature))

	if status == "analysis":

		try:
			#Reading input data using numpy module
			col_block, col_pot, col_tot, col_rot = loadtxt(dest_path+"/results/pigs.eng",unpack=True, usecols=[0,1,2,3])
			mean_pot     = np.sum(col_pot)/len(col_block)
			mean_tot     = np.sum(col_tot)/len(col_block)
			mean_rot     = np.sum(col_rot)/len(col_block)
			x2			 = np.multiply(col_pot, col_pot)
			y2			 = np.multiply(col_tot, col_tot)
			z2			 = np.multiply(col_rot, col_rot)
			mean_sq_pot  = np.sum(x2)/len(col_block)
			mean_sq_tot  = np.sum(y2)/len(col_block)
			mean_sq_rot  = np.sum(z2)/len(col_block)
			error_pot    = sqrt((mean_sq_pot-mean_pot*mean_pot)/len(col_block))
			error_tot    = sqrt((mean_sq_tot-mean_tot*mean_tot)/len(col_block))
			error_rot    = sqrt((mean_sq_rot-mean_rot*mean_rot)/len(col_block))
			print i, len(col_block)
			
			fanalyze.write(outputstring2(numbbeads,value,mean_pot,mean_tot,mean_rot,error_pot,error_tot,error_rot))
		except:
			print "no file ", i
			pass

if status == "submission":
	fsubmit.close()

if status == "analysis":
	fanalyze.close()
	call(["cat",file_output])
