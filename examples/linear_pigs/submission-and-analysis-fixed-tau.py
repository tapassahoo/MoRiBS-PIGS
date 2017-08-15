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

#NameOfServer        = "graham"
NameOfServer        = "nlogn"
NameOfPartition     = "ntapas"

#TypeCal             = 'PIGS'
TypeCal             = 'ENT'

#molecule            = "HFC60"                                                 
molecule            = "HF"                                                     
#molecule            = "H2"                                                   
molecule_rot        = "HF"                                                   

numbblocks	        = 100000
numbmolecules       = 2
numbpass            = 200
tau                 = 0.005

Rpt                 = 10.05
dipolemoment        = 1.826        #J. Chern. Phys. 73(5), 2319 (1980).

status_rhomat       = "Yes"                                                      
status_cagepot      = "No"                                                      
#RUNDIR              = "work"
RUNDIR              = "scratch"
RUNIN               = "nCPU"

loopStart           = 2
loopEnd             = 51
skip                = 2

preskip             = 00000
postskip            = 0
particleA           = int(numbmolecules/2)

ENT_TYPE = "SWAPTOUNSWAP"
#ENT_TYPE = "SWAP"
#ENT_TYPE = "BROKENPATH"
#ENT_TYPE = "REGULARPATH"

extra_file_name     = "-Passes"+str(numbpass)+"-NumTimes"
#extra_file_name     = ""

if (TypeCal == "PIGS"):
	file1_name      = "Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-e0vsbeads" 

if (TypeCal == "ENT"):
	file1_name      = "Entanglement-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-tau"+str(tau)+"Kinv-Blocks"+str(numbblocks)
	file1_name     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-ParticleA"+str(particleA)+"-e0vsbeads-"+ENT_TYPE 


file2_name          = ""                                                           #change param13
argument2           = "beads"                                                      #change param14
value_min           = 1                                                            #change param15
var                 = "beta"                                                       #change param16

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

#===============================================================================
#                                                                              |
#   Analysis of output files 												   |
#                                                                              |
#===============================================================================
if status == "analysis":
	if (TypeCal == "PIGS"):
		file_output      = "Energy-vs-"+str(var)+"-fixed-"
		file_output     += "tau"+str(tau)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+".txt"
		file_output_angularDOF  = "AngularDOF-vs-"+str(var)+"-fixed-"
		file_output_angularDOF += "tau"+str(tau)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output_angularDOF += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+".txt"
		SaveEnergy              = src_path+"/ResultsOfPIGS/"+file_output
		SaveAngDOFS             = src_path+"/ResultsOfPIGS/"+file_output_angularDOF
		call(["rm", SaveEnergy, SaveAngDOFS])

	if (TypeCal == "ENT"):
		file_output      = "Entropy-vs-"+str(var)+"-fixed-"
		file_output     += "tau"+str(tau)+"Kinv-Rpt"+str(Rpt)+"Angstrom-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
		file_output     += extra_file_name+"-System"+str(numbmolecules)+str(molecule)+"-ParticleA"+str(particleA)+"-preskip"+str(preskip)+"-postskip"+str(postskip)+"-"+ENT_TYPE+".txt"
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

if (TypeCal == "ENT"):
	numbmolecules  *= 2
	loopStart       = 5

step = [1.8 for i in range(loopStart, 1000, skip)]
# Loop over jobs
for i in range(loopStart, loopEnd, skip):                                                

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
		support.Submission(status, RUNDIR, dest_path, folder_run, src_path, run_file, dest_dir, Rpt, numbbeads, i, skip, step, temperature,numbblocks,numbpass,molecule_rot,numbmolecules,dipolemoment, status_rhomat, TypeCal, argument2, final_path, dest_pimc, RUNIN, particleA, NameOfServer, NameOfPartition, status_cagepot)


	if status == "analysis":

		variable          = beta
		try:
			if (TypeCal == "ENT"):
				fanalyze.write(support.GetAverageEntropy(numbbeads,variable,dest_dir,preskip,postskip,ENT_TYPE))
			else:
				fanalyze.write(support.GetAverageEnergy(numbbeads,variable,dest_dir,preskip,postskip))
				fanalyze_angularDOF.write(GetAverageOrientation(support.numbbeads,variable,dest_dir,preskip,postskip))
		except:
			pass

if status == "analysis":
    fanalyze.close()
    if (TypeCal == "ENT"):
        call(["cat",SaveEntropy])
    print
    print
    if (TypeCal == "PIGS"):
        fanalyze_angularDOF.close()
        call(["cat",SaveAngDOFS])
        call(["cat",SaveEnergy])

