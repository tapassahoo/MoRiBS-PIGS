import argparse
import decimal
import os
import sys
import time
from os import system
from subprocess import call

import numpy as np
from numpy import *

sys.path.append("../../examples/scripts")
import inputFile
import mypkg.pkgMoribs.support_without_parallel as support

parser = argparse.ArgumentParser(prog="script_submission_analysis_MoRiBS_without_parallel.py",description="It is a script file and used to submit jobs in a queue as well as analyze output data files. Note: Module support.py consists of many functions and it is not permitted to modify without consulting the developer - Dr. Tapas Sahoo. User can easily modify module inputFile.py to generate lists of beads (see Getbeads function), step lengths for rotational and translational motions, and levels for Bisection move (see class GetStepAndLevel) as needed.",epilog="Enjoy the program! :)")
parser.add_argument("cal", help="Type of calculation - it is a string: a) PIMC - Finite Temperature calculation by Path Integral Monte Carlo b) PIGS - Ground State Path Integral c) ENT - Entanglement by replica algorithm based on PIGS.", choices=["PIMC", "PIGS", "ENT"])
parser.add_argument("job", help="Type of a job: submission of new jobs or analyzing output files. It must be a string.", choices=["submission", "analysis", "rename"])
parser.add_argument("Molecule", help="Name of molecular system. E.g. - H2O, FCC-H2O, H2O@c60")
parser.add_argument("Rotor", help="Name of rotor. E.g. - HF, H2O. It is needed to save rotational density matrix.")
parser.add_argument("param", type=float, help="Fixed value of beta or tau.")
parser.add_argument("variable", help="Name of a variable: either beta or tau. It must be a string. Note: for finite temperature computations only the variable tau is needed.", choices=["tau", "beta"])

parser.add_argument("--MOVECOM", action="store_true", help="It allows translational motions of molecules or particles.")
parser.add_argument("--ROTMOVE", action="store_true", help="It allows rotational motions of molecules or particles.")
parser.add_argument("-N", help="Number of Molecules. It must be an integer.", type=int)
parser.add_argument("-Block", help="Number of Blocks. It must be an integer", type=int)
parser.add_argument("-Pass", help="Number of Passes. It must be an integer", type=int)
parser.add_argument("--scal", help="subtype of calculations - must be defined as a string in case of ENT.", default="EXTENDED_ENSMBL", choices=["EXTENDED_ENSMBL", "BROKENPATH"])
parser.add_argument("--RATIO", help="subtype of calculations - must be defined as a string in case of ENT. It applies ratio trick algorithm.", default="WR", choices=["WR", "WOR"])
parser.add_argument("-C", "--compiler", action="store_true", help="User can use it if the execution file (pimc) is already generated.")
parser.add_argument("--preskip", type=int, help="skips # of lines from the begining of an output file. It can be needed while analysis flag is open!", default=0)
parser.add_argument("--postskip", type=int, help="skips # of lines from the end of an output file. It can be needed while analysis flag is open!", default=0)
parser.add_argument("--RESTART", action="store_true", help="It is used to restart the code.")
parser.add_argument("-NR", type=int, help="Number of blocks extended by uesr.", default=0)
parser.add_argument("--Type", default="LINEAR", help="Specify your rotor type: LINEAR or NONLINEAR.")
parser.add_argument("-spin", "--SpinIsomer", type=int, help="It hvae three values -1, 0, 1 for spinless, para and ortho isomers.", default=-1)
parser.add_argument("-IM", "--Impurity", action="store", nargs=4, help="Use it if the system comprises various types of particles, molecules. IMPURITY[0] = type of particle (ATOM, MOLECULE, PLANAR, LINEAR or NONLINEAR); IMPURITY[1] = Name of the particle; IMPURITY[2] = Number of the particle; IMPURITY[3] = Either BOSE or BOLTZMANN.")

parser.add_argument("-R", "--Rpt", type=float, help="Distance between Centre of Masses of two molecules. Unit is Angstrom.", default=-1.0)
parser.add_argument("-d", "--DipoleMoment", type=float, help="Dipole Moment of a bipolar molecule in Debye.", default=-1.0)
parser.add_argument("-g", "--gFactor", type=float, help="It defines interaction strength.", default=-1.0)


parser.add_argument("--partition", help="allows to submit jobs in a specific cpu. It is a string. User does not need it.", default="ntapas")
parser.add_argument("--PPA", action="store_true", help="Inclussion of Pair Product Approximation. It is in the developing condition.")
parser.add_argument("-lmax", "--lmaxloop,max", help="Maximum l quantum number. It is needed for exact computations of Energy and Entropy of linear rotors in pinned to a chain.", default=0)
parser.add_argument("-ltotalmax", "--ltotalmax", help="Maximum lmax quantum number. It is needed for exact computations of Energy and Entropy of linear rotors in pinned to a chain.", default=0)
parser.add_argument("--CRYSTAL", action="store_true", help="Reads Lattice configurations and the corresponding dipolemoments.")
args = parser.parse_args()

# ===============================================================================
#																			   |
#	Some parameters for submission of jobs and analysis outputs.			   |
#	Change the parameters as you requied.									   |
#																			   |
# ===============================================================================
variableName = args.variable
#
TransMove = args.MOVECOM
RotMove = args.ROTMOVE
if (args.Impurity):
	impurity=args.Impurity
else:
	impurity=[""]


status = args.job
#
myhost = os.uname()[1]
if (myhost == "gra-login1") or (myhost == "gra-login2") or (myhost == "gra-login3"):
	NameOfServer = "graham"
elif ((myhost == "cedar1.cedar.computecanada.ca") or (myhost == "cedar2.cedar.computecanada.ca") or (myhost == "cedar3.cedar.computecanada.ca") or (myhost == "cedar4.cedar.computecanada.ca") or (myhost == "cedar5.cedar.computecanada.ca")):
    NameOfServer = "graham"
else:
	NameOfServer = "nlogn"
NameOfPartition = args.partition
#
TypeCal = args.cal
#
molecule = args.Molecule
molecule_rot = args.Rotor
#
# print 5/(support.bconstant(molecule_rot)/0.695)
# print 7/(support.bconstant(molecule_rot)/0.695)
# exit()
#
numbblocks = args.Block
numbmolecules1 = args.N
numbpass = args.Pass
#
if args.Rpt:
	Rpt = args.Rpt
if args.DipoleMoment:
	dipolemoment = args.DipoleMoment
if args.gFactor:
	gfact = args.gFactor
# if args.Rpt:
#	support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
# exit()

# Request to change
# User should change the following 4 lines as developer already have explained in the README file. If user wish to include cage potential, "No" flag must be turned into "Yes" by the user.

status_cagepot = False
# RUNDIR			  = "work"
RUNDIR = "scratch"
RUNIN = "nCPU"

preskip = args.preskip
postskip = args.postskip
if args.CRYSTAL == True:
	crystal = args.CRYSTAL
else:
	crystal = False
RotorType = args.Type

ENT_TYPE = args.scal
ENT_ALGR = args.RATIO

extra_file_name = extraName

src_dir = os.getcwd()
if (variableName == "tau"):
	parameterName = "beta"
	beta = args.param
	parameter = beta
	temperature = 1.0 / beta

if (variableName == "beta"):
	parameterName = "tau"
	tau = args.param
	parameter = tau

spin_isomer = args.SpinIsomer

if (spin_isomer == -1):
	preName = ""
if (spin_isomer == 0):
	preName = "-p-"
if (spin_isomer == 1):
	preName = "-o-"
molecule = preName+molecule

#Monte Carlo step size - Block start
step_COM = [0.1*i for i in range(20)] 
level_bisection = [1 for i in range(20)] 
step_COM_impurity = [0.1*i for i in range(20)] 
level_bisection_impurity = [1 for i in range(20)] 
step_rot = [0.1*i for i in range(20)]

steplevel = inputFile.GetStepAndLevel(molecule_rot, variableName, TypeCal)
step_COM = steplevel.step_trans
level_bisection = steplevel.level
if (len(impurity)==4):
	step_COM_impurity = steplevel.step_trans1
	level_bisection_impurity = steplevel.level1
step_rot = steplevel.step
#Monte Carlo step size - Block end

numbblocks_Restart1 = args.NR

# Request to change
# User should change the following 5 lines as developer already have explained in the README file.
user_name = os.getlogin()
input_dir = os.getcwd()+"/"
home = os.path.expanduser("~")
source_dir = "MoRiBS-PIGS/"
out_dir = "NameOfOutputDirectory/"
final_results_path = home + "/ResultsOf" + TypeCal + "/"
dir = os.path.dirname(final_results_path)
if not os.path.exists(dir):
	os.makedirs(dir)

source_dir_exe = "/home/" + user_name + "/" + source_dir
if status == "submission":

	if (RUNDIR == "scratch") or (NameOfServer == "graham"):
		dir_run_job = "/scratch/" + user_name + "/" + out_dir
		if NameOfServer == "graham":
			dir = os.path.dirname(dir_run_job)
			if not os.path.exists(dir):
				os.makedirs(dir)
	else:
		dir_run_job = "/work/" + user_name + "/" + out_dir

	execution_file = "/home/" + user_name + "/" + source_dir + "pimc"
	if not args.compiler:
		if not args.RESTART:
			support.makeexecutionfile(src_dir, TypeCal, ENT_TYPE, source_dir_exe)

if NameOfServer == "graham":
	dir_output = "/scratch/" + user_name + "/" + out_dir
else:
	dir_output = "/work/" + user_name + "/" + out_dir


if (TypeCal == "ENT"):
	maxloop = int(numbmolecules1/2)
	ext_ent=args.RATIO
else:
	maxloop = 1
	ext_ent=""

if (args.RATIO == "WR"):
	particleAList = np.arange(1, maxloop+1)
else:
	particleAList = [maxloop]
# particleAList = [9]

for particleA in particleAList:
	# Generating files for submission
	file1_name = support.GetFileNameSubmission(TypeCal,molecule_rot,TransMove,RotMove,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules1,molecule,ENT_TYPE,particleA,extra_file_name,crystal,impurity,ext_ent)

	if (status == "rename"):
		numbblocks_rename = args.NR
		file2_name = support.GetFileNameSubmission(TypeCal,molecule_rot,TransMove,RotMove,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks_rename,numbpass,numbmolecules1,molecule,ENT_TYPE,particleA,extra_file_name,crystal,impurity,ext_ent)

		if NameOfServer == "graham":
			dir_run_input_pimc = "/scratch/" + user_name + "/" + out_dir + file1_name + "PIMC"
			dir_input_pimc_renamed = "/scratch/" + user_name + "/" + out_dir + file2_name + "PIMC"
		else:
			dir_run_input_pimc = "/work/" + user_name + "/" + out_dir + file1_name + "PIMC"
			dir_input_pimc_renamed = "/work/" + user_name + "/" + out_dir + file2_name + "PIMC"

	# ===============================================================================
	#																			   |
	#	compilation of linden.f to generate rotational density matrix - linden.out |
	#	Yet to be generalized													   |
	#																			   |
	# ===============================================================================
	if status == "submission":
		if NameOfServer == "graham":
			dir_run_input_pimc = "/scratch/" + user_name + "/" + out_dir + file1_name + "PIMC"
		else:
			dir_run_input_pimc = "/work/" + user_name + "/" + out_dir + file1_name + "PIMC"

		if os.path.isdir(dir_run_input_pimc) == False:
			call(["rm", "-rf", dir_run_input_pimc])
			call(["mkdir", "-p", dir_run_input_pimc])

		if not args.RESTART:
			call(["cp", execution_file, dir_run_input_pimc])

		if RotorType == "LINEAR":
			if not args.RESTART:
				support.compile_rotmat(source_dir_exe, input_dir)

		if status_cagepot == True:
			if not args.RESTART:
				support.compile_cagepot(source_dir_exe, input_dir)
				support.cagepot(source_dir_exe)
				call(["mv", "hfc60.pot", dir_run_input_pimc])

	if (status == "analysis"):
		FileAnalysis = support.GetFileNameAnalysis(TypeCal,False,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules1,molecule,ENT_TYPE,preskip,postskip,extra_file_name,final_results_path,particleA, ENT_ALGR)

		if (TypeCal != "ENT"):
			if preskip >= numbblocks:
				print("")
				print("Warning!!!!!!!")
				print("============================================================================")
				print("Number of Blocks = " + str(numbblocks))
				print("Number of preskip= " + str(preskip))
				print("Error message: Number of preskip data must be less than Number of Blocks")
				print("============================================================================")
				exit()

			fanalyzeEnergy = open(FileAnalysis.SaveEnergy, "a")
			fanalyzeEnergy.write(support.fmtAverageEnergy(TypeCal, status, variableName))
			fanalyzeCorr = open(FileAnalysis.SaveCorr, "a")
			fanalyzeCorr.write(support.fmtAverageOrderParam(status, variableName))

		if ((status == "analysis") and ((TypeCal == "ENT") and (ENT_ALGR == "WOR"))):
			fanalyzeEntropy = open(FileAnalysis.SaveEntropy, "a")
			fanalyzeEntropy.write(support.fmtAverageEntropy(status, variableName, ENT_TYPE))

	if TypeCal == "ENT":
		numbmolecules = 2 * numbmolecules1
	else:
		numbmolecules = numbmolecules1

	list_nb = inputFile.Getbeads(TypeCal, variableName)

	iStep = 0
	for i in list_nb:

		if (TypeCal == "PIMC"):

			value = i

			if (variableName == "beta"):
				beta = tau*value
				temperature=1.0/beta
				variable=beta
			if (variableName == "tau"):
				tau=beta/value
				variable=tau

			numbbeads = value
			folder_run = file1_name + str(numbbeads)

			if status == "submission":

				if args.PPA:
					PPA1 = True
					lmax = args.lmax
					ltotalmax = args.ltotalmax
					support.GetTwoBodyDensity(Rpt,dipolemoment,numbbeads,lmax,ltotalmax,tau,molecule_rot)
					call(["mv", "PairDensity.txt", dir_run_input_pimc])
				else:
					PPA1 = False

				if args.RESTART:
					Restart1 = True
				else:
					Restart1 = False

				support.Submission(NameOfServer,status,TransMove,RotMove,RUNDIR,dir_run_job,folder_run,src_dir,execution_file,Rpt,numbbeads,i,step_rot,step_COM,level_bisection,temperature,numbblocks,numbpass,molecule_rot,numbmolecules,gfact,dipolemoment,TypeCal,ENT_TYPE, ENT_ALGR, dir_output,dir_run_input_pimc,RUNIN,particleA,NameOfPartition,status_cagepot,iStep,PPA1,user_name,out_dir,source_dir_exe,Restart1,numbblocks_Restart1,crystal,RotorType,spin_isomer,impurity, step_COM_impurity, level_bisection_impurity)

			if status == "analysis":

				final_dir_in_work = dir_output + folder_run
				try:
					fanalyzeEnergy.write(support.GetAverageEnergy(TypeCal,numbbeads,variable,final_dir_in_work,preskip,postskip,numbblocks))
					fanalyzeCorr.write(support.GetAverageOrderParam(TypeCal,numbmolecules,numbbeads, variable, final_dir_in_work, preskip, postskip,numbblocks))
				except:
					pass
		else:

			if (i % 2) != 0:
				value = i
			else:
				value = i + 1

			if variableName == "beta":
				beta = tau * (value - 1)
				temperature = 1.0 / beta
				variable = beta
			if variableName == "tau":
				tau = beta / (value - 1)
				variable = tau

			numbbeads = value
			folder_run = file1_name + str(numbbeads)

			if (status=="rename"):
				folder_rename = file2_name + str(numbbeads)
				support.GetRenamingFunc(dir_run_input_pimc,dir_input_pimc_renamed,dir_output,folder_run,folder_rename,src_dir)

			if (status=="submission"):

				if args.PPA:
					PPA1 = True
					lmax = args.lmax
					ltotalmax = args.ltotalmax
					support.GetTwoBodyDensity(Rpt,dipolemoment,numbbeads,lmax,ltotalmax,tau,molecule_rot)
					call(["mv", "PairDensity.txt", dir_run_input_pimc])
				else:
					PPA1 = False

				if args.RESTART:
					Restart1 = True
				else:
					Restart1 = False

				support.Submission(NameOfServer,status,TransMove,RotMove,RUNDIR,dir_run_job,folder_run,src_dir,execution_file,Rpt,numbbeads,i,step_rot,step_COM,level_bisection,temperature,numbblocks,numbpass,molecule_rot,numbmolecules,gfact,dipolemoment,TypeCal, ENT_TYPE, ENT_ALGR, dir_output, dir_run_input_pimc, RUNIN, particleA, NameOfPartition, status_cagepot, iStep, PPA1, user_name, out_dir, source_dir_exe, Restart1, numbblocks_Restart1, crystal, RotorType, spin_isomer,impurity, step_COM_impurity, level_bisection_impurity)

			if status == "analysis":

				final_dir_in_work = dir_output + folder_run
				#support.RemoveFiles(TypeCal, numbbeads, temperature, molecule_rot, RotorType, preskip, postskip, numbblocks, final_dir_in_work)
				try:
					if ((TypeCal == "ENT") and (ENT_ALGR == "WOR")):
						fanalyzeEntropy.write(support.GetAverageEntropy(numbbeads,variable,final_dir_in_work,preskip,postskip,numbblocks,ENT_TYPE))
					if (TypeCal != "ENT"):
						fanalyzeEnergy.write(support.GetAverageEnergy(TypeCal,numbbeads,variable,final_dir_in_work,preskip,postskip,numbblocks))
						fanalyzeCorr.write(support.GetAverageOrderParam(TypeCal,numbmolecules,numbbeads, variable, final_dir_in_work, preskip, postskip, numbblocks))
				except:
					pass
		iStep = iStep + 1

	if (status == "analysis") and (TypeCal != "ENT"):
		fanalyzeEnergy.close()
		fanalyzeCorr.close()
		call(["cat", FileAnalysis.SaveEnergy])
		print("")
		print("")
		call(["cat",FileAnalysis.SaveCorr])
		print("")
		print("")
		# =========================File Checking===============================#
		SavedFile = FileAnalysis.SaveEnergy
		support.FileCheck(TypeCal, list_nb, variableName, SavedFile)
		SavedFile = FileAnalysis.SaveCorr
		support.FileCheck(TypeCal, list_nb, variableName, SavedFile)

	if ((status == "analysis") and ((TypeCal == "ENT") and (ENT_ALGR == "WOR"))):
		fanalyzeEntropy.close()
		call(["cat", FileAnalysis.SaveEntropy])
		print("")
		print("")
		SavedFile = FileAnalysis.SaveEntropy
		support.FileCheck(TypeCal, list_nb, variableName, SavedFile)

if ((status == "analysis") and ((TypeCal == "ENT") and (ENT_TYPE == "EXTENDED_ENSMBL") and (ENT_ALGR == "WR"))):
	print("Final Entropy obtained by employing Ratio Trick")
	support.GetAverageEntropyRT(particleAList,TypeCal,molecule_rot,TransMove,RotMove,variableName,Rpt,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules1,molecule,ENT_TYPE,preskip,postskip,extra_file_name,dir_output,variable,crystal,final_results_path,impurity,ext_ent)
	"""
	support.GetEntropyRT(status, maxloop, TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, preskip, postskip, extra_file_name, dir_output, variable, crystal)
	"""
