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
#molecule            = "HF-C60"                                                     #change param1
molecule            = "HF"                                                         #change param1
#molecule            = "H2"                                                         #change param1
molecule_rot        = "HF"                                                         #change param2

numbblocks	        = 200                                                        #change param3
numbmolecules       = 1                                                            #change param4
tau                 = 0.001                                                        #change param5

Rpt                 = 10.0                                                         #change param6
dipolemoment        = 1.86

status              = "submission"                                                 #change param8
#status              = "analysis"                                                   #change param8
status_rhomat       = "Yes"                                                        #change param9 

nrange              = 8  			  						                       #change param10

file1_name          = "Rpt"+str(Rpt)+"Angstrom-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
file1_name         += "-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads"

file2_name          = ""                                                           #change param13
argument2           = "beads"                                                      #change param14
value_min           = 1                                                            #change param15
var                 = "beta"                                                       #change param16

src_path            = os.getcwd()
dest_path           = "/scratch/tapas/linear_rotors/"                                 #change param17
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"                  #change param18
trunc               = 5000


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
	file_output          = "Energy-vs-beta-"+str(numbmolecules)+"-"+str(molecule)
	file_output         += "-fixed-tau"+str(tau)+"-blocks"+str(numbblocks)+"Rpt"+str(Rpt)+"Angstrom-trunc"+str(trunc)+".txt"  
	file_output_angularDOF = "AngularDOF-vs-beta-"+str(numbmolecules)+"-"+str(molecule)
	file_output_angularDOF+= "-fixed-tau"+str(tau)+"-blocks"+str(numbblocks)+"Rpt"+str(Rpt)+"Angstrom-trunc"+str(trunc)+".txt"
	file_output_angularDOF1 = "AngularDOF-vs-beta-"+str(numbmolecules)+"-"+str(molecule)
	file_output_angularDOF1+= "-fixed-tau"+str(tau)+"-blocks"+str(numbblocks)+"Rpt"+str(Rpt)+"Angstrom-trunc"+str(trunc)+"for-zdir.txt"
	call(["rm", file_output, file_output_angularDOF,file_output_angularDOF1])

	fanalyze             = open(file_output, "a")           
	fanalyze.write(support.fmt_energy(status,var))

	fanalyze_angularDOF  = open(file_output_angularDOF, "a")           
	fanalyze_angularDOF.write(support.fmt_angle(status,var))
	fanalyze_angularDOF1  = open(file_output_angularDOF1, "a")           
	fanalyze_angularDOF1.write(support.fmt_angle1(status,var))

# Loop over jobs
for i in range(nrange):                                                  #change param19

	if (i>0):

		value        = pow(2,i) + value_min

		numbbeads    = support.dropzeros(value)
		beta         = tau*(value-1)
		temperature  = 1.0/beta

		folder_run   = file1_name+str(numbbeads)+file2_name
		dest_dir     = dest_path + folder_run 

		if status   == "submission":
			# Write submit file for the current cycle
			argument1     = Rpt
			level         = support.levels(numbbeads)
			step1         = 1.0;#step[jj]
			support.modify_input(temperature,numbbeads,numbblocks,molecule_rot,numbmolecules,argument1,level,step1,dipolemoment)

			if status_rhomat == "Yes":
				support.rotmat(molecule_rot,temperature,numbbeads)
	
			#job submission
			fname         = 'submit_'+str(i)
			fwrite        = open(fname, 'w')
	
			fwrite.write(support.jobstring(argument2,numbbeads,numbmolecules, dest_dir, molecule_rot, temperature, numbbeads))

			fwrite.close()
			call(["qsub", fname])
			os.chdir(src_path)

		if status == "analysis":

			variable          = beta
			fanalyze.write(support.outputstr_energy(numbbeads,variable,dest_dir,trunc))
			fanalyze_angularDOF.write(support.outputstr_angle(numbbeads,variable,dest_dir,trunc))
			fanalyze_angularDOF1.write(support.outputstr_angle1(numbbeads,variable,dest_dir,trunc))

if status == "analysis":
	fanalyze.close()
	fanalyze_angularDOF.close()
	fanalyze_angularDOF1.close()
	call(["cat",file_output])
	print
	print
	call(["cat",file_output_angularDOF])
	print
	print
	call(["cat",file_output_angularDOF1])
