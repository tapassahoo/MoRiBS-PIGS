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
status              = "submission"                                            
status              = "analysis"                                            

NameOfServer        = "nlogn"
NameOfPartition     = "ntapas"

#NameOfServer        = "graham"
#TypeCal             = 'PIMC'
#TypeCal             = 'PIGS'
TypeCal             = 'ENT'

#molecule            = "HFC60"                                                  
molecule            = "HF"                                                      
#molecule            = "H2"                                                    
molecule_rot        = "HF"

#print 5/(support.bconstant(molecule_rot)/0.695)
#print 7/(support.bconstant(molecule_rot)/0.695)
#exit()

numbblocks	        = 100000
numbmolecules       = 6
numbpass            = 100
beta     	        = 0.2

Rpt                 = 10.05
dipolemoment        = 1.826        #J. Chern. Phys. 73(5), 2319 (1980).

status_rhomat       = "Yes"                                                 
status_cagepot      = "No"                                                      
#RUNDIR              = "work"
RUNDIR              = "scratch"
RUNIN               = "nCPU"

loopStart           = 8
loopEnd             = 61
skip                = 2

preskip             = 1000
postskip            = 0

particleA           = int(numbmolecules/2)
ENT_TYPE = "SWAPTOUNSWAP"
#ENT_TYPE = "SWAP"
#ENT_TYPE = "BROKENPATH"
#ENT_TYPE = "REGULARPATH"
extra_file_name     = "-Passes"+str(numbpass)+"-NumTimes"
#extra_file_name     = ""

if (TypeCal == "PIGS"):
	file1_name      = "Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

if (TypeCal == "PIMC"):
	file1_name      = "PIMC-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

if (TypeCal == "ENT"):
	file1_name      = "Entanglement-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-beta"+str(beta)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads"+ENT_TYPE


file2_name          = ""                                                           #change param10
argument2           = "beads"                                                      #change param11
var                 = "tau"                                                        #change param13

src_path            = os.getcwd()
if (NameOfServer == "graham"):
    dest_path       = "/scratch/tapas/linear_rotors/"
else:
    dest_path       = "/work/tapas/linear_rotors/"
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"     

if status   == "submission":
	if RUNDIR == "scratch":
		dest_path   = "/scratch/tapas/linear_rotors/" 
			

if (NameOfServer == "graham"):
	final_path      = "/scratch/tapas/linear_rotors/"     
else:
	final_path      = "/work/tapas/linear_rotors/"             

temperature         = 1.0/beta   

#==================================== MCStep ===================================# 
if (molecule_rot == "H2"):
	#step           = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
	#step           = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
	step            = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6
	level           = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

if (molecule_rot == "HF"):
	#step          	= [2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,1.8,1.8,1.8,1.6,1.6,1.6,1.6,1.4,1.4,1.8,1.7,1.6,1.5,1.5,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4]  
	#step          	= [2.0,2.0,2.0,1.8,1.8,1.8,1.6,1.6,1.6,1.6,1.4,1.4,1.8,1.7,1.6,1.5,1.5,1.4,1.4,1.4,1.4,1.4,1.4,1.4,1.4]  # beads 41 to 61
	step           	= [1.6,1.6,1.6]  # beads 57 to 61
					# 2 HF beta 0.1 K-1 #change param6 for 10.05 Angstrom and Dipole Moment 1.86 Debye PIGS
	#step           = [0.7, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.5, 1.5, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.0, 1.0, 1.0, 1.0, 1.0]  
					# 2 HF beta 0.2 K-1 #change param6 for 10 Angstrom PIMC
	level           = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

#===============================================================================
#                                                                              |
#   compilation of linden.f to generate rotational density matrix - linden.out |
#   Yet to be generalized                                                      |
#                                                                              |
#===============================================================================
if status == "submission":
	if (NameOfServer == "graham"):
		dest_pimc = "/scratch/tapas/linear_rotors/"+file1_name+"PIMC"
	else:
		dest_pimc = "/work/tapas/linear_rotors/"+file1_name+"PIMC"
	call(["rm", "-rf",  dest_pimc])
	call(["mkdir", "-p", dest_pimc])
	call(["cp", run_file, dest_pimc])
	if status_rhomat == "Yes":
		support.compile_rotmat()
	if status_cagepot == "Yes":
		support.compile_cagepot()
		support.cagepot();
		call(["mv", "hfc60.pot", dest_pimc])

if status == "analysis":
	if (TypeCal == "PIGS"):
		file_output             = "Energy-vs-"+str(var)+"-fixed-"
		file_output            += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output            += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+".txt"

		file_output_angularDOF  = "AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+".txt"
		
		SaveEnergy              = src_path+"/ResultsOfPIGS/"+file_output
		SaveAngDOFS             = src_path+"/ResultsOfPIGS/"+file_output_angularDOF
		call(["rm", SaveEnergy, SaveAngDOFS])
	
	if (TypeCal == "PIMC"):
		file_output             = "PIMC-Energy-vs-"+str(var)+"-fixed-"
		file_output            += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output            += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+".txt"

		file_output_angularDOF  = "PIMC-AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+".txt"

		SaveEnergy              = src_path+"/ResultsOfPIMC/"+file_output
		SaveAngDOFS             = src_path+"/ResultsOfPIMC/"+file_output_angularDOF
		call(["rm", SaveEnergy, SaveAngDOFS])
	
	if (TypeCal == "ENT"):
		file_output             = "Entropy-vs-"+str(var)+"-fixed-"
		file_output            += "beta"+str(beta)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output            += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-ParticleA"+str(particleA)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+"-"+ENT_TYPE+".txt"

		SaveEntropy             = src_path+"/ResultsOfPIGSENT/"+file_output
		call(["rm", SaveEntropy])
	
	
	if (TypeCal == "ENT"):
		fanalyze             = open(SaveEntropy, "a")
		fanalyze.write(support.fmtAverageEntropy(status,var,ENT_TYPE))
	else:
		fanalyze             = open(SaveEnergy, "a")           
		fanalyze.write(support.fmtAverageEnergy(status,var))
		fanalyze_angularDOF  = open(SaveAngDOFS, "a")           
		fanalyze_angularDOF.write(support.fmtAverageOrientation(status,var))
	#support.FileOutput(status, TypeCal, var, beta, Rpt, dipolemoment, numbblocks, numbmolecules, molecule, trunc)

if (TypeCal == "ENT"):
	numbmolecules  *= 2
	loopStart       = 20

# Loop over jobs
#list_nb = [8,16,32,64,96,128]
#skip    = 2
#for i in list_nb:

iStep = 0
for i in range(loopStart, loopEnd, skip):                                                

	if (TypeCal == 'PIMC'):
		if ((i%2) == 0):
			value = i
			tau          = beta/value
			numbbeads    = value
			folder_run   = file1_name+str(numbbeads)+file2_name
			dest_dir     = dest_path + folder_run 
			print(dest_dir)

			if status   == "submission":
				support.Submission(status, RUNDIR, dest_path, folder_run, src_path, run_file, dest_dir, Rpt, numbbeads, value, step, level, temperature,numbblocks,numbpass,molecule_rot,numbmolecules,dipolemoment, status_rhomat, TypeCal, argument2, final_path, dest_pimc, RUNIN, particleA, NameOfServer, NameOfPartition, status_cagepot, iStep)

			if status == "analysis":

				variable          = tau
				try:
					if (TypeCal == "ENT"):
						fanalyze.write(support.GetAverageEntropy(numbbeads,variable,dest_dir,preskip,postskip))
					else:
						fanalyze.write(support.GetAverageEnergy(numbbeads,variable,dest_dir,preskip,postskip))
						fanalyze_angularDOF.write(GetAverageOrientation(support.numbbeads,variable,dest_dir,preskip,postskip))
				except:
					pass
	else:
		if (i % skip == 0 ):

			if ((i % 2) != 0):
				value    = i
			else:
				value    = i+1
			tau          = beta/(value-1)
			numbbeads    = value
			folder_run   = file1_name+str(numbbeads)+file2_name
			dest_dir     = dest_path + folder_run 
			print(dest_dir)

			if status   == "submission":
				support.Submission(status, RUNDIR, dest_path, folder_run, src_path, run_file, dest_dir, Rpt, numbbeads, i, step, level, temperature,numbblocks,numbpass,molecule_rot,numbmolecules,dipolemoment, status_rhomat, TypeCal, argument2, final_path, dest_pimc, RUNIN, particleA, NameOfServer, NameOfPartition, status_cagepot, iStep)

			if status == "analysis":

				variable          = tau
				try:
					if (TypeCal == "ENT"):
						fanalyze.write(support.GetAverageEntropy(numbbeads,variable,dest_dir,preskip,postskip,ENT_TYPE))
					else:
						fanalyze.write(support.GetAverageEnergy(numbbeads,variable,dest_dir,preskip,postskip))
						fanalyze_angularDOF.write(support.GetAverageOrientation(numbbeads,variable,dest_dir,preskip,postskip))
				except:
					pass
	iStep = iStep + 1

if status == "analysis":
#	support.FileClose(TypeCal)
	fanalyze.close()
	if (TypeCal == "ENT"):
		call(["cat",SaveEntropy])
	print("")
	print("")
	if (TypeCal == "PIGS" or TypeCal == "PIMC"):
		fanalyze_angularDOF.close()
		call(["cat",SaveAngDOFS])
		call(["cat",SaveEnergy])
