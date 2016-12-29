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
molecule            = "HF-C60"                                                     #change param1
#molecule            = "H2"                                                         #change param1
molecule_rot        = "HF"
#print 5/(support.bconstant(molecule_rot)/0.695)
#print 7/(support.bconstant(molecule_rot)/0.695)
#exit()
numbblocks	        = 80000                                                        #change param2
numbmolecules       = 2                                                            #change param3
temperature	        = 5.0                                                         #change param4
Rpt                 = 3.5                                                          #change param7

skip		        = 2                                                            #change param8
#status              = "submission"                                                 #change param9
status              = "analysis"                                                   #change param10
status_rhomat       = "Yes"                                                        #change param10 

nrange             = 8 #41		  						                           #change param12

if (molecule_rot == "H2"):
	#step                = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K         #change param6
	#step                = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K         #change param6
	step                = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K         #change param6
	file1_name          = "Rpt"+str(Rpt)+"Angstrom-Temperature"+str(temperature)+"K-Blocks"+str(numbblocks)+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads"

if (molecule_rot == "HF"):
	#step                = [3.0,3.0,3.0,3.0,3.0,1.8,1.2,0.80,0.2]                            #change param6
	step                = [3.0,3.0,3.0,3.0,3.0,2.0,1.2,0.80,0.2]                            #change param6
	file1_name          = "newTemperature"+str(temperature)+"K-Blocks"+str(numbblocks)+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

file2_name          = ""                                                           #change param10
argument2           = "beads"                                                      #change param11
value_min           = 1                                                            #change param12
var                 = "tau"                                                        #change param13

src_path            = os.getcwd()
dest_path           = "/work/tapas/linear_rotors/"                                 #change param13
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"                  #change param14

beta                = 1.0/temperature   

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

	file_input = "Input-parameters-vs-beads-"+str(numbmolecules)+"-"+str(molecule)+"-fixed-beta"+str(beta)+"-blocks"+str(numbblocks)+".txt"
	call(["rm", file_input])
	fsubmit    = open(file_input, "a")  


#===============================================================================
#                                                                              |
#   Analysis of output files 												   |
#                                                                              |
#===============================================================================
if status == "analysis":
	file_output          = "newEnergy-vs-tau-"+str(numbmolecules)+"-"+str(molecule)+"-fixed-beta"+str(beta)+"-blocks"+str(numbblocks)+".txt"  
	file_output_angularDOF = "newAngularDOF-vs-tau-"+str(numbmolecules)+"-"+str(molecule)+"-fixed-beta"+str(beta)+"-blocks"+str(numbblocks)+".txt"  
	call(["rm", file_output, file_output_angularDOF])

	fanalyze             = open(file_output, "a")           
	fanalyze.write(support.formatting(status,var))

	fanalyze_angularDOF  = open(file_output_angularDOF, "a")           
	fanalyze_angularDOF.write(support.formatting1(status,var))



# Loop over jobs
for i in range(nrange):                                                  #change param13
	if (i>0):
 
		value        = pow(2,i) + value_min

		numbbeads    	 = support.dropzeros(value)
	
		tau              = beta/(value-1)

		folder_run   = file1_name+str(numbbeads)+file2_name
		dest_dir      = dest_path + folder_run 

		if status == "submission":
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
			step1         = step[i]
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

			fsubmit.write(support.outputstring1(numbbeads,tau,temperature))

		if status == "analysis":

			variable = tau
			try:
				#Reading input data using numpy module
				col_block, col_pot, col_tot, col_rot = loadtxt(dest_dir+"/results/pigs.eng",unpack=True, usecols=[0,1,2,3])
				mean_pot   = np.sum(col_pot)/len(col_block)
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
			
				fanalyze.write(support.outputstring2(numbbeads,variable,mean_pot,mean_tot,mean_rot,error_pot,error_tot,error_rot))


				col_block, col_costheta, col_theta, col_phi = loadtxt(dest_dir+"/results/pigs.dof",unpack=True, usecols=[0,1,2,3])
				mean_costheta = np.sum(col_costheta)/len(col_block)
				mean_theta    = np.sum(col_theta)/len(col_block)
				mean_phi      = np.sum(col_phi)/len(col_block)
				x2			  = np.multiply(col_costheta, col_costheta)
				y2			  = np.multiply(col_theta, col_theta)
				z2			  = np.multiply(col_phi, col_phi)
				mean_sq_costheta = np.sum(x2)/len(col_block)
				mean_sq_theta = np.sum(y2)/len(col_block)
				mean_sq_phi   = np.sum(z2)/len(col_block)
				error_costheta   = sqrt((mean_sq_costheta - mean_costheta*mean_costheta)/len(col_block))
				error_theta   = sqrt((mean_sq_theta - mean_theta*mean_theta)/len(col_block))
				error_phi     = sqrt((mean_sq_phi - mean_phi*mean_phi)/len(col_block))
			
				fanalyze_angularDOF.write(support.outputstring2(numbbeads,variable,mean_costheta,mean_theta,mean_phi,error_costheta,error_theta,error_phi))
			except:
				print "no file ", i
				pass

if status == "submission":
	fsubmit.close()

if status == "analysis":
	fanalyze.close()
	fanalyze_angularDOF.close()
	call(["cat",file_output])
	print
	print
	call(["cat",file_output_angularDOF])
