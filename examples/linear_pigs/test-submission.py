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
molecule            = "HF"                                                         #change param1
molecule_rot        = "HF"                                                         #change param2

numbblocks	        = 200                                                        #change param3
numbmolecules       = 2                                                            #change param4
tau                 = 0.001                                                        #change param5

Rpt                 = 20.0                                                         #change param6
dipolemoment        = 1.86

status              = "submission"                                                 #change param8
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
final_path           = "/work/tapas/linear_rotors/"                                 #change param17


#===============================================================================
#                                                                              |
#   compilation of linden.f to generate rotational density matrix - linden.out |
#   Yet to be generalized                                                      |
#                                                                              |
#===============================================================================
if status == "submission":
	if status_rhomat == "Yes":
		support.compile_rotmat()

# Loop over jobs
for i in range(nrange):                                                  #change param19

	if (i>0):

		value        = pow(2,i) + value_min

		numbbeads    = support.dropzeros(value)
		beta         = tau*(value-1)
		temperature  = 1.0/beta

		folder_run   = file1_name+str(numbbeads)+file2_name
		dest_dir     = dest_path + folder_run 
		final_dir     = final_path + folder_run 
		call(["rm", "-rf", final_dir])

		if status   == "submission":
			# Write submit file for the current cycle
			argument1     = Rpt
			level         = support.levels(numbbeads)
			step1         = 1.0;#step[jj]
			support.modify_input(temperature,numbbeads,numbblocks,molecule_rot,numbmolecules,argument1,level,step1,dipolemoment)
			input_file    = "qmc_beads"+str(numbbeads)+".input"
			print input_file
			call(["mv", "qmc.input", input_file])

			if status_rhomat == "Yes":
				support.rotmat(molecule_rot,temperature,numbbeads)
	
			#job submission
			fname         = 'submit_'+str(i)
			fwrite        = open(fname, 'w')
	
			fwrite.write(support.jobstring_scratch(argument2,numbbeads,numbmolecules, dest_dir, molecule_rot, temperature, numbbeads, final_dir))

			fwrite.close()
			call(["qsub", fname])
			#os.chdir(src_path)
