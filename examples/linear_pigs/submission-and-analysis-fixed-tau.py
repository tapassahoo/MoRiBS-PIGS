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
molecule            = "HF"                                                     #change param1
#molecule            = "H2"                                                         #change param1
molecule_rot        = "HF"                                                         #change param2

numbblocks	        = 50000                                                        #change param3
numbmolecules       = 2                                                            #change param4
tau                 = 0.001                                                       #change param5

Rpt                 = 3.5                                                          #change param6

status              = "submission"                                                 #change param8
status              = "analysis"                                                   #change param8
status_rhomat       = "Yes"                                                        #change param9 

if (molecule_rot == "H2"):
	nrange          = 8  			  						                       #change param10
if (molecule_rot == "HF"):
	nrange          = 8  			  						                       #change param10

if (molecule_rot == "H2"):
	#step            = [1., 1., 3.0, 1., 1.5, 1., 1.0, 1., 0.7, 1., 0.5, 1., 0.5, 1., 0.4, 1., 0.4, 1., 0.3, 1., 0.3]  #tau = 0.02       #change param11
	step            = [1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.2, 1.1, 1.1, 1.1]  #tau = 0.001      #change param11
	#step            = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  #tau = 0.0005      #change param11
	file1_name      = "Rpt"+str(Rpt)+"Angstrom-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads"

if (molecule_rot == "HF"):
	#step            = [1., 1., 2.0, 1., 1.0, 1., 0.6, 1., 0.4, 1., 0.3, 1., 0.2, 1., 0.2, 1., 0.4, 1., 0.3, 1., 0.3]  #tau = 0.02      #change param12
	#step            = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8]  #tau = 0.001      #change param11
	step            = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  #tau = 0.002      #change param11
	file1_name      = "tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

file2_name          = ""                                                           #change param13
argument2           = "beads"                                                      #change param14
value_min           = 1                                                            #change param15
var                 = "Beta"                                                       #change param16

src_path            = os.getcwd()
dest_path           = "/work/tapas/linear_rotors/"                                 #change param17
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"                  #change param18


#===============================================================================
#                                                                              |
#   compilation of linden.f to generate rotational density matrix - linden.out |
#   Yet to be generalized                                                      |
#                                                                              |
#===============================================================================
if status == "submission":
	if status_rhomat == "Yes":
		path_enter_linden = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/linear_prop/"
		os.chdir(path_enter_linden)
		call(["make", "clean"])
		call(["make"])
		path_exit_linden  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"
		os.chdir(path_exit_linden)

#===============================================================================
#                                                                              |
#   Analysis of output files 												   |
#                                                                              |
#===============================================================================
if status == "analysis":
	file_output          = "Energy-vs-beta-"+str(numbmolecules)+"-"+str(molecule)+"-fixed-tau"+str(tau)+"-blocks"+str(numbblocks)+".txt"  
	file_output_angularDOF = "AngularDOF-vs-beta-"+str(numbmolecules)+"-"+str(molecule)+"-fixed-tau"+str(tau)+"-blocks"+str(numbblocks)+".txt"  
	call(["rm", file_output, file_output_angularDOF])

	fanalyze             = open(file_output, "a")           
	fanalyze.write(support.fmt_energy(status,var))

	fanalyze_angularDOF  = open(file_output_angularDOF, "a")           
	fanalyze_angularDOF.write(support.fmt_angle(status,var))

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
			step1         = 3.0;#step[jj]
			support.modify_input(temperature,numbbeads,numbblocks,molecule_rot,numbmolecules,argument1,level,step1)

			if status_rhomat == "Yes":
				support.rotmat(molecule_rot,temperature,numbbeads)
	
			#job submission
			fname         = 'submit_'+str(i)
			fwrite        = open(fname, 'w')
	
			fwrite.write(support.jobstring(argument2,numbbeads,numbmolecules))

			fwrite.close()
			call(["qsub", fname, ])
			os.chdir(src_path)

		if status == "analysis":

			variable          = beta
			try:
				#Reading input data using numpy module
				col_block, col_pot, col_tot, col_rot = loadtxt(dest_dir+"/results/pigs.eng",unpack=True, usecols=[0,1,2,3])
				mean_pot      = np.sum(col_pot)/len(col_block)
				mean_tot      = np.sum(col_tot)/len(col_block)
				mean_rot      = np.sum(col_rot)/len(col_block)
				x2			  = np.multiply(col_pot, col_pot)
				y2			  = np.multiply(col_tot, col_tot)
				z2			  = np.multiply(col_rot, col_rot)
				mean_sq_pot   = np.sum(x2)/len(col_block)
				mean_sq_tot   = np.sum(y2)/len(col_block)
				mean_sq_rot   = np.sum(z2)/len(col_block)
				error_pot     = sqrt((mean_sq_pot-mean_pot*mean_pot)/len(col_block))
				error_tot     = sqrt((mean_sq_tot-mean_tot*mean_tot)/len(col_block))
				error_rot     = sqrt((mean_sq_rot-mean_rot*mean_rot)/len(col_block))
				print i, len(col_block)
			
				fanalyze.write(support.outputstr_energy(numbbeads,variable,mean_pot,mean_tot,mean_rot,error_pot,error_tot,error_rot))


				col_block, col_costheta, col_theta = loadtxt(dest_dir+"/results/pigs.dof",unpack=True, usecols=[0,1,2])
				mean_costheta = np.sum(col_costheta)/len(col_block)
				mean_theta    = np.sum(col_theta)/len(col_block)
				x2			  = np.multiply(col_costheta, col_costheta)
				y2			  = np.multiply(col_theta, col_theta)
				mean_sq_costheta = np.sum(x2)/len(col_block)
				mean_sq_theta = np.sum(y2)/len(col_block)
				error_costheta   = sqrt((mean_sq_costheta - mean_costheta*mean_costheta)/len(col_block))
				error_theta   = sqrt((mean_sq_theta - mean_theta*mean_theta)/len(col_block))
			
				fanalyze_angularDOF.write(support.outputstr_angle(numbbeads,variable,mean_costheta,mean_theta,error_costheta,error_theta))
			except:
				print "no file ", i
				pass

if status == "analysis":
	fanalyze.close()
	fanalyze_angularDOF.close()
	call(["cat",file_output])
	print
	print
	call(["cat",file_output_angularDOF])
