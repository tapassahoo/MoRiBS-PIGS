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
TypeCal             = 'PIMC'
#TypeCal             = 'PIGS'
#TypeCal             = 'ENT'

#molecule            = "HF-C60"                                                    #change param1
molecule            = "HF"                                                         #change param1
#molecule            = "H2"                                                        #change param1
molecule_rot        = "HF"

#print 5/(support.bconstant(molecule_rot)/0.695)
#print 7/(support.bconstant(molecule_rot)/0.695)
#exit()

numbblocks	        = 40000                                                        #change param2
numbmolecules       = 2                                                           #change param3
numbpass            = 50
beta     	        = 0.2                                                         #change param4

Rpt                 = 10.0                                                        #change param7
dipolemoment        = 1.86
skip                = 10

status              = "submission"                                                 #change param9
#status              = "analysis"                                                   #change param10
status_rhomat       = "Yes"                                                        #change param10 
#RUNDIR              = "work"
RUNDIR              = "scratch"

nrange              = 101 		  						                           #change param12

if (TypeCal == "PIGS"):
	file1_name      = "Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

if (TypeCal == "PIMC"):
	file1_name      = "PIMC-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 


if (TypeCal == "ENT"):
	numbmolecules  *= 2
	file1_name      = "Entanglement-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 


file2_name          = ""                                                           #change param10
argument2           = "beads"                                                      #change param11
value_min           = 1                                                            #change param12
var                 = "tau"                                                        #change param13

src_path            = os.getcwd()
dest_path           = "/work/tapas/linear_rotors/"                                 #change param13
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"                   #change param14

if status   == "submission":
	if RUNDIR == "scratch":
		dest_path       = "/scratch/tapas/linear_rotors/" 
		final_path      = "/work/tapas/linear_rotors/"                                 #change param17

temperature         = 1.0/beta   

trunc               = numbblocks
preskip             = 0

#==================================== MCStep ===================================# 
if (molecule_rot == "H2"):
	#step           = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
	#step           = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
	step            = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6

if (molecule_rot == "HF"):
	step           = [0.7,1.5,1.5,1.5,1.5,1.5,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.0,1.0,1.0,1.0,1.0]  # 2 HF beta 0.2 K-1 #change param6 for 10 Angstrom

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
		file_output             = "Energy-vs-"+str(var)+"-fixed-"
		file_output            += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output            += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"

		file_output_angularDOF  = "AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"

		file_output_angularDOF1  = "AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF1 += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF1 += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+"-for-zdir.txt"

		call(["rm", file_output, file_output_angularDOF, file_output_angularDOF1])


	if (TypeCal == "PIMC"):
		file_output             = "PIMC-Energy-vs-"+str(var)+"-fixed-"
		file_output            += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output            += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"

		file_output_angularDOF  = "PIMC-AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+".txt"

		file_output_angularDOF1  = "PIMC-AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF1 += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF1 += "-System"+str(numbmolecules)+str(molecule)+"-trunc"+str(trunc)+"-for-zdir.txt"

		call(["rm", file_output, file_output_angularDOF, file_output_angularDOF1])


	if (TypeCal == "ENT"):
		file_output      = "Entropy-vs-"+str(var)+"-fixed-"
		file_output     += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
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
for i in range(nrange):                                                  #change param13

	if (TypeCal == 'PIMC'):

		if (i>1 and i%skip == 0):
			value = i
			tau          = beta/value
			numbbeads    = value
			folder_run   = file1_name+str(numbbeads)+file2_name
			dest_dir     = dest_path + folder_run 
			print dest_dir

			if status == "submission":
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
				istep         = i/skip
				step1         = step[istep]
				support.modify_input(temperature,numbbeads,numbblocks,numbpass,molecule_rot,numbmolecules,argument1,level,step1,dipolemoment)
				if status_rhomat == "Yes":
					support.rotmat(TypeCal,molecule_rot,temperature,numbbeads)
	
				#job submission
				fname         = 'jobsubmit_'+str(i)
				fwrite        = open(fname, 'w')
				print fname
	
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
				#call(["qsub", fname])
				if RUNDIR != "scratch":
					os.chdir(src_path)

			if status == "analysis":

				variable          = tau
				try:
					fanalyze.write(support.outputstr_energy(numbbeads,variable,dest_dir,trunc,preskip))
					fanalyze_angularDOF.write(support.outputstr_angle(numbbeads,variable,dest_dir,trunc,preskip))
					fanalyze_angularDOF1.write(support.outputstr_angle1(numbbeads,variable,dest_dir,trunc,preskip))
				except:
					pass

	else:
		if (i>1 and i % skip == 0 ):

			if i % 2 != 0:
				value        = i
			else:
				value        = i+value_min
			tau          = beta/(value-1)
			numbbeads    = value
			folder_run   = file1_name+str(numbbeads)+file2_name
			dest_dir     = dest_path + folder_run 
			print dest_dir

			if status == "submission":
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
				istep         = i/skip
				step1         = step[istep]
				support.modify_input(temperature,numbbeads,numbblocks,numbpass,molecule_rot,numbmolecules,argument1,level,step1,dipolemoment)
				if status_rhomat == "Yes":
					support.rotmat(TypeCal,molecule_rot,temperature,numbbeads)
	
				#job submission
				fname         = 'jobsubmit_'+str(i)
				fwrite        = open(fname, 'w')
				print fname
	
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
				#call(["qsub", fname])
				if RUNDIR != "scratch":
					os.chdir(src_path)

			if status == "analysis":

				variable          = tau
				try:
					if (TypeCal == "ENT"):
						fanalyze.write(support.outputstr_entropy(numbbeads,variable,dest_dir,trunc,preskip))
					else:
						fanalyze.write(support.outputstr_energy(numbbeads,variable,dest_dir,trunc,preskip))
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
