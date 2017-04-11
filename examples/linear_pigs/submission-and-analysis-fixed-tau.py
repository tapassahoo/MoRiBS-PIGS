#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
TypeCal             = 'PIGS'
TypeCal             = 'ENT'

#status              = "submission"                                                
status              = "analysis"                                                  
status_rhomat       = "Yes"                                                      
#RUNDIR              = "work"
RUNDIR              = "scratch"

#molecule            = "HF-C60"                                                 
molecule            = "HF"                                                     
#molecule            = "H2"                                                   
molecule_rot        = "HF"                                                   

numbblocks	        = 10000                                                 
numbmolecules       = 2                                                    
numbpass            = 50
skip                = 5

tau                 = 0.002                                               

Rpt                 = 10.0                                               
dipolemoment        = 1.86


nrange              = 101 #31  			  						        

if (TypeCal == "PIGS"):
	file1_name      = "Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

if (TypeCal == "ENT"):
	numbmolecules   *= 2
	file1_name      = "Entanglement-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

if (TypeCal == "PIMC"):
	file1_name      = "Finite-Temperature-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 


file2_name          = ""                                                           #change param13
argument2           = "beads"                                                      #change param14
value_min           = 1                                                            #change param15
var                 = "beta"                                                       #change param16

src_path            = os.getcwd()
dest_path           = "/work/tapas/linear_rotors/"                                 #change param17
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"                  #change param18

if status   == "submission":
	if RUNDIR == "scratch":
		dest_path   = "/scratch/tapas/linear_rotors/" 
		final_path  = "/work/tapas/linear_rotors/"                                 #change param17

trunc               = numbblocks
preskip             = 1000


#===============================================================================
#                                                                              |
#   compilation of linden.f to generate rotational density matrix - linden.out |
#   Yet to be generalized                                                      |
#                                                                              |
#===============================================================================
if status == "submission":
	if status_rhomat == "Yes":
		support.compile_rotmat()
#===============================================================================
#                                                                              |
#   Analysis of output files 												   |
#                                                                              |
#===============================================================================
if status == "analysis":
	if (TypeCal == "PIGS"):
		file_output      = "Energy-vs-"+str(var)+"-fixed-"
		file_output     += "tau"+str(tau)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output     += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"
		file_output_angularDOF  = "AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF += "tau"+str(tau)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"
		file_output_angularDOF1  = "AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF1 += "tau"+str(tau)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF1 += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+"-for-zdir.txt"
		call(["rm", file_output, file_output_angularDOF,file_output_angularDOF1])

	if (TypeCal == "ENT"):
		file_output      = "Entropy-vs-"+str(var)+"-fixed-"
		file_output     += "tau"+str(tau)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output     += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"
		call(["rm", file_output])

	if (TypeCal == "ENT"):
		fanalyze             = open(file_output, "a")
		fanalyze.write(support.fmt_entropy(status,var))
	else:
		fanalyze             = open(file_output, "a")           
		fanalyze.write(support.fmt_energy(status,var))
		fanalyze_angularDOF  = open(file_output_angularDOF, "a")           
		fanalyze_angularDOF.write(support.fmt_angle(status,var))
		fanalyze_angularDOF1  = open(file_output_angularDOF1, "a")           
		fanalyze_angularDOF1.write(support.fmt_angle1(status,var))

# Loop over jobs
for i in range(nrange):                                                  #change param19

	'''
	if (i>0):

		#value        = pow(2,i) + value_min
		value        = beads_list[i]
	'''
	if (i>1 and i % skip == 0 ):

		if i % 2 != 0:
			value        = i
		else:
			value        = i+value_min


		numbbeads    = value #support.dropzeros(value)
		beta         = tau*(value-1)
		temperature  = 1.0/beta

		folder_run   = file1_name+str(numbbeads)+file2_name
		dest_dir     = dest_path + folder_run 

		if status   == "submission":
			if RUNDIR != "scratch":
				os.chdir(dest_path)
				call(["rm", "-rf", folder_run])
				call(["mkdir", folder_run])
				call(["mkdir", "-p", folder_run+"/results"])
				os.chdir(src_path)


				# copy files to running folder
				src_file      = src_path + "/IhRCOMC60.xyz"
				call(["cp", src_file, dest_dir])
				call(["cp", run_file, dest_dir])

				src_file      = src_path + "/qmc_run.input"
				call(["cp", src_file, dest_dir])

				# Write submit file for the current cycle
				os.chdir(dest_dir)

			argument1     = Rpt
			level         = support.levels(numbbeads)
			step1         = 1.0;#step[jj]
			support.modify_input(temperature,numbbeads,numbblocks,numbpass,molecule_rot,numbmolecules,argument1,level,step1,dipolemoment)

			if status_rhomat == "Yes":
				support.rotmat(molecule_rot,temperature,numbbeads)
	
			#job submission
			fname         = 'jobsubmit_'+str(i)
			fwrite        = open(fname, 'w')
	

			if RUNDIR == "scratch":
				os.chdir(final_path)
				call(["rm", "-rf", folder_run])
				os.chdir(src_path)
				final_dir     = final_path + folder_run 
				call(["rm", "-rf", final_dir])
				input_file    = "qmc"+argument2+str(i)+".input"
				print input_file
				call(["mv", "qmc.input", input_file])
				fwrite.write(support.jobstring_scratch(argument2,i,numbmolecules, dest_dir, molecule_rot, temperature, numbbeads, final_dir))
			else: 
				fwrite.write(support.jobstring(argument2,numbbeads,numbmolecules))
				

			fwrite.close()
			call(["qsub", fname])
			if RUNDIR != "scratch":
				os.chdir(src_path)

		if status == "analysis":

			variable          = beta
			try:
				if (TypeCal == "ENT"):
					fanalyze.write(support.outputstr_entropy(numbbeads,variable,dest_dir,trunc,preskip))
				else:
					fanalyze.write(support.outputstr_energy(numbbeads,variable,dest_dir,trunc,preskip))
					print dest_dir
					fanalyze_angularDOF.write(support.outputstr_angle(numbbeads,variable,dest_dir,trunc,preskip))
					fanalyze_angularDOF1.write(support.outputstr_angle1(numbbeads,variable,dest_dir,trunc,preskip))
			except:
				pass

if status == "analysis":
	fanalyze.close()
	if (TypeCal != "ENT"):
		fanalyze_angularDOF.close()
		fanalyze_angularDOF1.close()
	call(["cat",file_output])
	'''
	print
	print
	call(["cat",file_output_angularDOF])
	print
	print
	call(["cat",file_output_angularDOF1])
	'''
