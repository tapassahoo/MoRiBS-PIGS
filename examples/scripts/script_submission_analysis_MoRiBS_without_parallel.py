#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support_without_parallel as support
import inputFile
import sys
import argparse

parser = argparse.ArgumentParser(description='It is a script file, written in Python, used to submit jobs in a queue as well as analyze output data files. Note: Module support.py consists of many functions and it is not permitted to modify without consulting the developer - Dr. Tapas Sahoo. User can easily modify module inputFile.py to generate lists of beads (see Getbeads function), step lengths for rotational and translational motions, and levels for Bisection move (see class GetStepAndLevel) as needed.')
parser.add_argument("-d", "--DipoleMoment", type=float, help="Dipole Moment of a bipolar molecule in Debye.", default = -1.0)
parser.add_argument("-g", "--gFactor", type=float, help="It defines interaction strength.", default = -1.0)
parser.add_argument("-R", "--Rpt", type=float, help="Inter molecular spacing.", default = -1.0)
parser.add_argument("variable", help="Name of a variable: either beta or tau. It must be a string. Note: for finite temperature computations only the variable tau is needed.", choices =["tau","beta"])
parser.add_argument("job", help="Type of a job: submission of new jobs or analyzing output files. It must be a string.", choices = ["submission", "analysis"])
parser.add_argument("cal", help="Type of calculation - it is a string: a) PIMC - Finite Temperature calculation by Path Integral Monte Carlo b) PIGS - Ground State Path Integral c) ENT - Entanglement by replica algorithm based on PIGS.", choices = ["PIMC", "PIGS", "ENT"])
parser.add_argument("--scal", help="subtype of calculations - must be defined as a string in case of ENT.", default = "SWAPTOUNSWAP", choices = ["SWAPTOUNSWAP", "BROKENPATH"])
parser.add_argument("--RATIO", help="subtype of calculations - must be defined as a string in case of ENT. It applies ratio trick algorithm.", action="store_true")
parser.add_argument("-N", help="Number of Molecules. It must be an integer.", type = int)
parser.add_argument("-Block", help="Number of Blocks. It must be an integer", type = int)
parser.add_argument("-Pass", help="Number of Passes. It must be an integer", type = int)
parser.add_argument("--MOVECOM", action="store_true", help="It allows translational motions of molecules or particles.")
parser.add_argument("--ROTMOVE", action="store_true", help="It allows rotational motions of molecules or particles.")
parser.add_argument("--partition", help="allows to submit jobs in a specific cpu. It is a string. User does not need it.", default = "ntapas")
parser.add_argument("Molecule", help="Name of molecular system. E.g. - H2O, FCC-H2O, H2O@c60")
parser.add_argument("--PPA", action="store_true", help="Inclussion of Pair Product Approximation. It is in the developing condition.")
parser.add_argument("-lmax", "--lmaxloop,max", help="Maximum l quantum number. It is needed for exact computations of Energy and Entropy of linear rotors in pinned to a chain.", default = 0)
parser.add_argument("-ltotalmax", "--ltotalmax", help="Maximum lmax quantum number. It is needed for exact computations of Energy and Entropy of linear rotors in pinned to a chain.", default = 0)
parser.add_argument("Rotor", help="Name of rotor. E.g. - HF, H2O. It is needed to save rotational density matrix.")
parser.add_argument("param", type=float, help="Fixed value of beta or tau.")
parser.add_argument("--preskip", type=int, help="skips # of lines from the begining of an output file. It can be needed while analysis flag is open!", default = 0)
parser.add_argument("--postskip", type=int, help="skips # of lines from the end of an output file. It can be needed while analysis flag is open!", default = 0)
parser.add_argument("-C", "--compiler", action="store_true", help="User can use it if the execution file (pimc) is already generated.")
parser.add_argument("--RESTART", action="store_true", help="It is used to restart the code.")
parser.add_argument("-NR", type = int, help="Number of blocks extended by uesr.", default = 0)
parser.add_argument("--CRYSTAL", action="store_true", help="Reads Lattice configurations and the corresponding dipolemoments.")
parser.add_argument("--Type", default = "LINEAR", help="Specify your rotor type: LINEAR or NONLINEAR.")
args = parser.parse_args()

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
variableName        = args.variable
#
TransMove           = args.MOVECOM
RotMove             = args.ROTMOVE
#
status              = args.job
#
#Request to change
#If user wish to run MoRiBs in graham.computecanada.ca, just replace "NameOfServer = "nlogn"" by "NameOfServer = "graham""
NameOfServer        = "nlogn"
#NameOfServer        = "graham"
NameOfPartition     = args.partition
#
TypeCal             = args.cal
#
molecule            = args.Molecule
molecule_rot        = args.Rotor
#
#print 5/(support.bconstant(molecule_rot)/0.695)
#print 7/(support.bconstant(molecule_rot)/0.695)
#exit()
#
numbblocks	        = args.Block
numbmolecules1      = args.N
numbpass            = args.Pass
#
if(args.Rpt):
	Rpt             	= args.Rpt
if(args.DipoleMoment):
	dipolemoment        = args.DipoleMoment
if(args.gFactor):
	gfact           = args.gFactor
#if args.Rpt:
#	support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
#exit()

#Request to change
#User should change the following 4 lines as developer already have explained in the README file. If user wish to include cage potential, "No" flag must be turned into "Yes" by the user. 

status_cagepot      = False                                                      
#RUNDIR              = "work"
RUNDIR              = "scratch"
RUNIN               = "nCPU"

preskip             = args.preskip
postskip            = args.postskip
if (args.CRYSTAL == True):
	crystal         = args.CRYSTAL
else: 
	crystal         = False
RotorType           = args.Type

ENT_TYPE 			= args.scal

extra_file_name     = extraName

src_dir             = os.getcwd()
if (variableName == "tau"):
	parameterName   = "beta"
	beta            = args.param
	parameter       = beta
	temperature     = 1.0/beta   

if (variableName == "beta"):
	parameterName   = "tau"
	tau             = args.param
	parameter       = tau

steplevel           = inputFile.GetStepAndLevel(molecule_rot,variableName)
step_COM            = steplevel.step_trans
step_rot	        = steplevel.step
level_bisection     = steplevel.level
numbblocks_Restart1 = args.NR

#Request to change
#User should change the following 5 lines as developer already have explained in the README file.
user_name           = "tapas"                   
source_dir          = "MoRiBS-PIGS/" 
out_dir             = "NameOfOutputDirectory/"
input_dir           = "/home/tapas/MoRiBS-PIGS/examples/scripts/"
final_results_path  = "/home/"+user_name+"/ResultsOf"+TypeCal+"/"
dir = os.path.dirname(final_results_path)
if not os.path.exists(dir):
	os.makedirs(dir)

source_dir_exe = "/home/"+user_name+"/"+source_dir
if status   == "submission":

	if (RUNDIR == "scratch") or (NameOfServer == "graham"):
		dir_run_job = "/scratch/"+user_name+"/"+out_dir 
		if (NameOfServer == "graham"):
			dir = os.path.dirname(dir_run_job)
			if not os.path.exists(dir):
				os.makedirs(dir)
	else:	
		dir_run_job = "/work/"+user_name+"/"+out_dir

	execution_file  = "/home/"+user_name+"/"+source_dir+"pimc"     
	if not args.compiler:
		if not args.RESTART:
			support.makeexecutionfile(src_dir,TypeCal,ENT_TYPE, source_dir_exe)

if (NameOfServer == "graham"):
	dir_output      = "/scratch/"+user_name+"/"+out_dir     
else:
	dir_output      = "/work/"+user_name+"/"+out_dir            


if (TypeCal == "ENT"):
	maxloop = int(numbmolecules1/2)
else:
	maxloop = 1

if (args.RATIO):	
	particleAList = np.arange(1, maxloop+1)
else:
	particleAList = [maxloop]

for particleA in particleAList:
	#==================================Generating files for submission================#
	file1_name = support.GetFileNameSubmission(TypeCal, molecule_rot, TransMove, RotMove, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, particleA, extra_file_name, crystal)
	#===============================================================================
	#                                                                              |
	#   compilation of linden.f to generate rotational density matrix - linden.out |
	#   Yet to be generalized                                                      |
	#                                                                              |
	#===============================================================================
	if status == "submission":
		if (NameOfServer == "graham"):
			dir_run_input_pimc = "/scratch/"+user_name+"/"+out_dir+file1_name+"PIMC"
		else:
			dir_run_input_pimc = "/work/"+user_name+"/"+out_dir+file1_name+"PIMC"

		if (os.path.isdir(dir_run_input_pimc) == False):
			call(["rm", "-rf",  dir_run_input_pimc])
			call(["mkdir", "-p", dir_run_input_pimc])

		if not args.RESTART:
			call(["cp", execution_file, dir_run_input_pimc])

		if (RotorType == "LINEAR"):
			if not args.RESTART:
				support.compile_rotmat(source_dir_exe, input_dir)

		if (status_cagepot == True):
			if not args.RESTART:
				support.compile_cagepot(source_dir_exe, input_dir)
				support.cagepot(source_dir_exe);
				call(["mv", "hfc60.pot", dir_run_input_pimc])

	if ((status == "analysis") and (TypeCal != "ENT")):
		FileAnalysis = support.GetFileNameAnalysis(TypeCal, False, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA)
		
		if (preskip >= numbblocks):
			print("")
			print("Warning!!!!!!!")
			print("============================================================================")
			print("Number of Blocks = "+str(numbblocks))
			print("Number of preskip= "+str(preskip))
			print("Error message: Number of preskip data must be less than Number of Blocks")
			print("============================================================================")
			exit()

		fanalyzeEnergy       = open(FileAnalysis.SaveEnergy, "a")           
		fanalyzeEnergy.write(support.fmtAverageEnergy(TypeCal,status,variableName))
		fanalyzeCorr         = open(FileAnalysis.SaveCorr, "a")           
		fanalyzeCorr.write(support.fmtAverageOrientation(status,variableName))
		fanalyzeTotalCorr    = open(FileAnalysis.SaveTotalCorr, "a")           
		fanalyzeXCorr        = open(FileAnalysis.SaveXCorr, "a")           
		fanalyzeYCorr        = open(FileAnalysis.SaveYCorr, "a")           
		fanalyzeZCorr        = open(FileAnalysis.SaveZCorr, "a")           
		fanalyzeXYCorr       = open(FileAnalysis.SaveXYCorr,"a")           

	if (TypeCal == "ENT"):
		numbmolecules  = 2*numbmolecules1
	else:
		numbmolecules  = numbmolecules1

	list_nb = inputFile.Getbeads(TypeCal, variableName)

	iStep = 0
	for i in list_nb:

		if (TypeCal == 'PIMC'):

			if i % 2 == 0:
				value    = i
			else:
				value    = i+1

			if (variableName == "beta"):
				beta     = tau*value
				temperature = 1.0/beta
				variable = beta
			if (variableName == "tau"):
				tau      = beta/value
				variable = tau

			numbbeads    = value
			folder_run   = file1_name+str(numbbeads)

			if status   == "submission":

				if args.PPA:
					PPA1 = True
					lmax = args.lmax
					ltotalmax = args.ltotalmax
					support.GetTwoBodyDensity(Rpt, dipolemoment, numbbeads, lmax, ltotalmax, tau, molecule_rot)
					call(["mv", "PairDensity.txt", dir_run_input_pimc])
				else: 
					PPA1 = False

				if args.RESTART:
					Restart1 = True
				else:
					Restart1 = False

				support.Submission(status,TransMove, RotMove,RUNDIR, dir_run_job, folder_run, src_dir, execution_file, Rpt, numbbeads, i, step_rot, step_COM, level_bisection, temperature, numbblocks, numbpass, molecule_rot, numbmolecules, gfact, dipolemoment, TypeCal, dir_output, dir_run_input_pimc, RUNIN, particleA, NameOfPartition, status_cagepot, iStep, PPA1, user_name, out_dir, source_dir_exe, Restart1, numbblocks_Restart1,crystal,RotorType)

			if status == "analysis":

				final_dir_in_work = dir_output+folder_run
				try:
					fanalyzeEnergy.write(support.GetAverageEnergy(TypeCal,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeCorr.write(support.GetAverageOrientation(numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeTotalCorr.write(support.GetAverageCorrelation("TotalCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeXCorr.write(support.GetAverageCorrelation("XCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeYCorr.write(support.GetAverageCorrelation("YCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeZCorr.write(support.GetAverageCorrelation("ZCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
					fanalyzeXYCorr.write(support.GetAverageCorrelation("XYCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
				except:
					pass
		else:

			if ((i % 2) != 0):
				value    = i
			else:
				value    = i+1

			if (variableName == "beta"):
				beta     = tau*(value-1)
				temperature = 1.0/beta
				variable = beta
			if (variableName == "tau"):
				tau      = beta/(value-1)
				variable = tau

			numbbeads    = value
			folder_run   = file1_name+str(numbbeads)

			if status   == "submission":

				if args.PPA:
					PPA1 = True
					lmax = args.lmax
					ltotalmax = args.ltotalmax
					support.GetTwoBodyDensity(Rpt, dipolemoment, numbbeads, lmax, ltotalmax, tau, molecule_rot)
					call(["mv", "PairDensity.txt", dir_run_input_pimc])
				else: 
					PPA1 = False

				if args.RESTART:
					Restart1 = True
				else:
					Restart1 = False

				support.Submission(status,TransMove, RotMove,RUNDIR, dir_run_job, folder_run, src_dir, execution_file, Rpt, numbbeads, i, step_rot, step_COM, level_bisection, temperature, numbblocks, numbpass, molecule_rot, numbmolecules, gfact, dipolemoment, TypeCal, dir_output, dir_run_input_pimc, RUNIN, particleA, NameOfPartition, status_cagepot, iStep, PPA1, user_name, out_dir, source_dir_exe, Restart1, numbblocks_Restart1, crystal,RotorType)

			if status == "analysis":

				final_dir_in_work = dir_output+folder_run
				try:
					if (TypeCal != "ENT"):
						fanalyzeEnergy.write(support.GetAverageEnergy(TypeCal,numbbeads,variable,final_dir_in_work,preskip,postskip))
						fanalyzeCorr.write(support.GetAverageOrientation(numbbeads,variable,final_dir_in_work,preskip,postskip))
						'''
						fanalyzeTotalCorr.write(support.GetAverageCorrelation("TotalCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
						fanalyzeXCorr.write(support.GetAverageCorrelation("XCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
						fanalyzeYCorr.write(support.GetAverageCorrelation("YCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
						fanalyzeZCorr.write(support.GetAverageCorrelation("ZCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
						fanalyzeXYCorr.write(support.GetAverageCorrelation("XYCorr", numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip))
						'''
				except:
					pass
		iStep = iStep + 1

	if ((status == "analysis") and (TypeCal != "ENT")):
		fanalyzeEnergy.close()
		fanalyzeCorr.close()
		fanalyzeTotalCorr.close()
		fanalyzeXCorr.close()
		fanalyzeYCorr.close()
		fanalyzeZCorr.close()
		fanalyzeXYCorr.close()
		call(["cat",FileAnalysis.SaveEnergy])
		print("")
		print("")
		call(["cat",FileAnalysis.SaveCorr])
		print("")
		print("")
		call(["cat",FileAnalysis.SaveTotalCorr])
#=========================File Checking===============================#
		SavedFile = FileAnalysis.SaveEnergy
		support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
		SavedFile = FileAnalysis.SaveCorr
		support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
		SavedFile = FileAnalysis.SaveTotalCorr
		support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
		SavedFile = FileAnalysis.SaveXCorr
		support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
		SavedFile = FileAnalysis.SaveYCorr
		support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
		SavedFile = FileAnalysis.SaveZCorr
		support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
		SavedFile = FileAnalysis.SaveXYCorr
		support.FileCheck(TypeCal,list_nb,variableName,SavedFile)
# END ========

if (status == "analysis") and (TypeCal == "ENT"):
	print("Final Entropy obtained by employing Ratio Trick")
	support.GetAverageEntropyRT(particleAList, TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, preskip, postskip, extra_file_name, dir_output, variable, crystal, final_results_path)
	exit()
	'''
	support.GetEntropyRT(status, maxloop, TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, preskip, postskip, extra_file_name, dir_output, variable, crystal)
	'''
