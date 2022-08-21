import argparse
import decimal
import os
import sys
import time
from os import system
from subprocess import call
import getpass

import numpy as np
from numpy import *

sys.path.append("../../examples/scripts")
import get_beads_and_mc_steps as mc
import mypkg.moribs_runner.support as support

parser = argparse.ArgumentParser(prog="submission_analysis.py",
                                 description="It is used to submit \
								 jobs in a queue and analyze output \
								 data files. Note: Module support.py \
								 consists of multiple functions \
								 and is not permitted to modify \
					             without consulting the developer - \
								 Dr Tapas Sahoo. Users should modify \
								 module /MoRiBS-PIGS/examples/scripts/get_beads_and_mc_steps.py \
								 to generate lists of beads, step lengths \
								 for rotational and translational motions, \
								 and levels for bisection move associated \
								 with each bead.",
                                 epilog="Enjoy the program! :)")
parser.add_argument("job",
					 help="Type of jobs: submission of new jobs or analyzing \
					 output files.", 
					 choices=["submission", "analysis", "rename"], 
				     type=str)
parser.add_argument("-C",
					"--compiler", 
					 action="store_true", 
					 help="Users can use it if the execution file \
					 (/MoRiBS-PIGS/pimc) is already generated.")
parser.add_argument("method", 
					 help="The methodology should be specified as \
					 one of PIMC, PIGS, and ENT. PIMC is for finite \
					 temperature calculation; PIGS is for ground state \
					 calculation; ENT is for entanglement by replica \
					 algorithm based on PIGS.", 
					 choices=["PIMC", "PIGS", "ENT"], 
					 type=str)
parser.add_argument("--scal", \
					 help="The subtype of calculations - \
					 must be defined in the case of ENT.", 
					 default="EXTENDED_ENSMBL", 
					 choices=["EXTENDED_ENSMBL", "BROKENPATH"])
parser.add_argument("--RATIO", 
				     help="The subtype of calculations - \
					 must be defined in the case of ENT. \
					 WR and WOR stand for with and without \
					 ratio trick algorithm.", 
					 default="WR", 
					 choices=["WR", "WOR"])
parser.add_argument("Molecule", 
					 help="Name of a molecular system. \
					 E.g., H2, HF, H2O, FCC-H2O, H2O@c60",
					 type=str)
parser.add_argument("Rotor", 
					 help="Name of the rotor. E.g., HF, H2O. \
					 It is needed to call the rotational density matrix file.",
					 type=str)
parser.add_argument("param", 
					type=float, 
					help="Fixed value of beta or tau.")
parser.add_argument("variable", 
					help="Name of the variable: beta (for param=tau) 
					or tau (param=beta). Note: for finite temperature 
					computations, only the variable tau is needed.", 
					choices=["tau", "beta"],
					type=str)
parser.add_argument("--MOVECOM", 
					action="store_true", 
					help="It allows translational motions of 
					molecules or particles.")
parser.add_argument("--ROTMOVE", 
					action="store_true", 
					help="It allows rotational motions of molecules.")
parser.add_argument("--Type", 
					default="LINEAR", 
					help="Users should specify the type of the rotor 
					either LINEAR or NONLINEAR.",
					choices=["LINEAR","NONLINEAR"])
parser.add_argument("-spin", 
					"--SpinIsomer", 
					type=int, help="It has three values. These are -1, 0, and 1 
					for spinless, para and ortho isomers, 
					respectively. Users should mention one of them.", 
					default=-1, 
					choices=[-1,0,1])
parser.add_argument("-N", 
					help="Number of Molecules.", 
					type=int)
parser.add_argument("-Block", 
					help="Number of Blocks required for Monte Carlo move.", 
					type=int)
parser.add_argument("-Pass", 
					help="Number of Passes required for Monte Carlo move.", 
					type=int)
parser.add_argument("-R", 
					"--Rpt", 
					type=float, 
					help="Distance between Centre of Masses of two 
					molecules. The unit is Angstrom.", 
					default=-1.0)
parser.add_argument("-d", 
					"--DipoleMoment", 
					type=float, 
					help="It defines the dipole moment of a linear 
					polar molecule in Debye.", 
					default=-1.0)
parser.add_argument("-g", 
					"--gFactor", 
					type=float, 
					help="It defines interaction strength.", 
					default=-1.0)
parser.add_argument("--RESTART", 
					action="store_true", 
					help="It is used to restart the code.")
parser.add_argument("-NR", 
					type=int, 
					help="The number of Monte Carlo blocks users 
					wish to extend further is mentioned.", 
					default=0)
parser.add_argument("--preskip", 
					type=int, help="The number of blocks in outputs 
					skipped from the beginning. It can be necessary 
					while the analysis flag is open to remove 
					pre-equilibrated data.", 
					default=0)
parser.add_argument("--postskip", 
					type=int, 
					help="The number of blocks in outputs skipped 
					from the end of the file. It can be necessary 
					to increase the error bar. Ensure the analysis 
					flag is open.", default=0)
parser.add_argument("-IM", 
					"--Impurity", 
					action="store", 
					nargs=4, 
					help="Use it if the system comprises various 
					types of particles and molecules. 
					IMPURITY[0] = type of particle (ATOM, MOLECULE, 
					PLANAR, LINEAR or NONLINEAR); 
					IMPURITY[1] = Name of the particle; 
					IMPURITY[2] = Number of the particle; 
					IMPURITY[3] = Either BOSE or BOLTZMANN.")
parser.add_argument("--partition", 
					help="It allows submitting jobs in a specific CPU. 
					The user does not need it.", 
					default="ntapas")
parser.add_argument("--CRYSTAL", 
					action="store_true", 
					help="It is required to read lattice configurations 
					and the associated dipole moments. 
					It is in developing condition.")
args = parser.parse_args()

# ===============================================================================
#																			   |
#	Some parameters for submission of jobs and analysis outputs.			   |
#	Change the parameters as you requied.									   |
#																			   |
# ===============================================================================
var_name = args.variable
#
tanslational_move = args.MOVECOM
rotational_move = args.ROTMOVE
if (args.Impurity):
	impurity=args.Impurity
else:
	impurity=[""]
print("Hi, I'm here")
exit()


status = args.job
#
myhost = os.uname()[1]
myhost = myhost[0:3]
if ((myhost == "gra") or (myhost == "ced")):
	server_name = "graham"
else:
	server_name = "nlogn"
partition_name = args.partition
#
method = args.method
#
molecule = args.Molecule
molecule_rot = args.Rotor
#
numbblocks = args.Block
numbmolecules1 = args.N
numbpass = args.Pass
#
if args.Rpt:
	rpt_val = args.Rpt
if args.DipoleMoment:
	dipolemoment = args.DipoleMoment
if args.gFactor:
	gfact = args.gFactor
# if args.Rpt:
#	support.GetrAndgFactor(molecule_rot, rpt_val, dipolemoment)
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
if (var_name == "tau"):
	parameterName = "beta"
	beta = args.param
	parameter = beta
	temperature = 1.0 / beta

if (var_name == "beta"):
	parameterName = "tau"
	tau = args.param
	parameter = tau

spin_isomer = args.SpinIsomer

if (spin_isomer == -1):
	prefix_name = ""
if (spin_isomer == 0):
	prefix_name = "-p-"
if (spin_isomer == 1):
	prefix_name = "-o-"
molecule = prefix_name+molecule

#Monte Carlo step size - Block start
step_COM = [0.1*i for i in range(20)] 
level_bisection = [1 for i in range(20)] 
step_COM_impurity = [0.1*i for i in range(20)] 
level_bisection_impurity = [1 for i in range(20)] 
step_rot = [0.1*i for i in range(20)]

steplevel = mc.GetStepAndLevel(molecule_rot, var_name, method)
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
user_name = getpass.getuser()
input_dir = os.getcwd()+"/"
home = os.path.expanduser("~")
source_dir = "MoRiBS-PIGS/"
out_dir = "NameOfOutputDirectory/"
final_results_path = home + "/ResultsOf" + method + "/"
dir = os.path.dirname(final_results_path)
if not os.path.exists(dir):
	os.makedirs(dir)

source_dir_exe = "/home/" + user_name + "/" + source_dir
if status == "submission":

	if (RUNDIR == "scratch") or (server_name == "graham"):
		dir_run_job = "/scratch/" + user_name + "/" + out_dir
		if server_name == "graham":
			dir = os.path.dirname(dir_run_job)
			if not os.path.exists(dir):
				os.makedirs(dir)
	else:
		dir_run_job = "/work/" + user_name + "/" + out_dir

	execution_file = "/home/" + user_name + "/" + source_dir + "pimc"
	if not args.compiler:
		if not args.RESTART:
			support.makeexecutionfile(src_dir, method, ENT_TYPE, source_dir_exe)

if server_name == "graham":
	dir_output = "/scratch/" + user_name + "/" + out_dir
else:
	dir_output = "/work/" + user_name + "/" + out_dir


if (method == "ENT"):
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
	file1_name = support.GetFileNameSubmission(method,molecule_rot,tanslational_move,rotational_move,rpt_val,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules1,molecule,ENT_TYPE,particleA,extra_file_name,crystal,impurity,ext_ent)

	if (status == "rename"):
		numbblocks_rename = args.NR
		file2_name = support.GetFileNameSubmission(method,molecule_rot,tanslational_move,rotational_move,rpt_val,gfact,dipolemoment,parameterName,parameter,numbblocks_rename,numbpass,numbmolecules1,molecule,ENT_TYPE,particleA,extra_file_name,crystal,impurity,ext_ent)

		if server_name == "graham":
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
		if server_name == "graham":
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
		FileAnalysis = support.GetFileNameAnalysis(method,False,molecule_rot,tanslational_move,rotational_move,var_name,rpt_val,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules1,molecule,ENT_TYPE,preskip,postskip,extra_file_name,final_results_path,particleA, ENT_ALGR)

		if (method != "ENT"):
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
			fanalyzeEnergy.write(support.fmtAverageEnergy(method, status, var_name))
			fanalyzeCorr = open(FileAnalysis.SaveCorr, "a")
			fanalyzeCorr.write(support.fmtAverageOrderParam(status, var_name))

		if ((status == "analysis") and ((method == "ENT") and (ENT_ALGR == "WOR"))):
			fanalyzeEntropy = open(FileAnalysis.SaveEntropy, "a")
			fanalyzeEntropy.write(support.fmtAverageEntropy(status, var_name, ENT_TYPE))

	if method == "ENT":
		numbmolecules = 2 * numbmolecules1
	else:
		numbmolecules = numbmolecules1

	list_nb = mc.Getbeads(method, var_name)

	iStep = 0
	for i in list_nb:

		if (method == "PIMC"):

			value = i

			if (var_name == "beta"):
				beta = tau*value
				temperature=1.0/beta
				variable=beta
			if (var_name == "tau"):
				tau=beta/value
				variable=tau

			numbbeads = value
			folder_run = file1_name + str(numbbeads)

			if status == "submission":

				if args.RESTART:
					Restart1 = True
				else:
					Restart1 = False

				support.Submission(server_name,status,tanslational_move,rotational_move,RUNDIR,dir_run_job,folder_run,src_dir,execution_file,rpt_val,numbbeads,i,step_rot,step_COM,level_bisection,temperature,numbblocks,numbpass,molecule_rot,numbmolecules,gfact,dipolemoment,method,ENT_TYPE, ENT_ALGR, dir_output,dir_run_input_pimc,RUNIN,particleA,partition_name,status_cagepot,iStep,user_name,out_dir,source_dir_exe,Restart1,numbblocks_Restart1,crystal,RotorType,spin_isomer,impurity, step_COM_impurity, level_bisection_impurity)

			if status == "analysis":

				final_dir_in_work = dir_output + folder_run
				try:
					fanalyzeEnergy.write(support.GetAverageEnergy(method,numbbeads,variable,final_dir_in_work,preskip,postskip,numbblocks))
					fanalyzeCorr.write(support.GetAverageOrderParam(method,numbmolecules,numbbeads, variable, final_dir_in_work, preskip, postskip,numbblocks))
				except:
					pass
		else:

			if (i % 2) != 0:
				value = i
			else:
				value = i + 1

			if var_name == "beta":
				beta = tau * (value - 1)
				temperature = 1.0 / beta
				variable = beta
			if var_name == "tau":
				tau = beta / (value - 1)
				variable = tau

			numbbeads = value
			folder_run = file1_name + str(numbbeads)

			if (status=="rename"):
				folder_rename = file2_name + str(numbbeads)
				support.GetRenamingFunc(dir_run_input_pimc,dir_input_pimc_renamed,dir_output,folder_run,folder_rename,src_dir)

			if (status=="submission"):

				if args.RESTART:
					Restart1 = True
				else:
					Restart1 = False

				support.Submission(server_name,status,tanslational_move,rotational_move,RUNDIR,dir_run_job,folder_run,src_dir,execution_file,rpt_val,numbbeads,i,step_rot,step_COM,level_bisection,temperature,numbblocks,numbpass,molecule_rot,numbmolecules,gfact,dipolemoment,method, ENT_TYPE, ENT_ALGR, dir_output, dir_run_input_pimc, RUNIN, particleA, partition_name, status_cagepot, iStep, user_name, out_dir, source_dir_exe, Restart1, numbblocks_Restart1, crystal, RotorType, spin_isomer,impurity, step_COM_impurity, level_bisection_impurity)

			if status == "analysis":

				final_dir_in_work = dir_output + folder_run
				#support.RemoveFiles(method, numbbeads, temperature, molecule_rot, RotorType, preskip, postskip, numbblocks, final_dir_in_work)
				try:
					if ((method == "ENT") and (ENT_ALGR == "WOR")):
						fanalyzeEntropy.write(support.GetAverageEntropy(numbbeads,variable,final_dir_in_work,preskip,postskip,numbblocks,ENT_TYPE))
					if (method != "ENT"):
						fanalyzeEnergy.write(support.GetAverageEnergy(method,numbbeads,variable,final_dir_in_work,preskip,postskip,numbblocks))
						fanalyzeCorr.write(support.GetAverageOrderParam(method,numbmolecules,numbbeads, variable, final_dir_in_work, preskip, postskip, numbblocks))
				except:
					pass
		iStep = iStep + 1

	if (status == "analysis") and (method != "ENT"):
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
		support.FileCheck(method, list_nb, var_name, SavedFile)
		SavedFile = FileAnalysis.SaveCorr
		support.FileCheck(method, list_nb, var_name, SavedFile)

	if ((status == "analysis") and ((method == "ENT") and (ENT_ALGR == "WOR"))):
		fanalyzeEntropy.close()
		call(["cat", FileAnalysis.SaveEntropy])
		print("")
		print("")
		SavedFile = FileAnalysis.SaveEntropy
		support.FileCheck(method, list_nb, var_name, SavedFile)

if ((status == "analysis") and ((method == "ENT") and (ENT_TYPE == "EXTENDED_ENSMBL") and (ENT_ALGR == "WR"))):
	print("Final Entropy obtained by employing Ratio Trick")
	support.GetAverageEntropyRT(particleAList,method,molecule_rot,tanslational_move,rotational_move,var_name,rpt_val,gfact,dipolemoment,parameterName,parameter,numbblocks,numbpass,numbmolecules1,molecule,ENT_TYPE,preskip,postskip,extra_file_name,dir_output,variable,crystal,final_results_path,impurity,ext_ent)
	"""
	support.GetEntropyRT(status, maxloop, method, molecule_rot, tanslational_move, rotational_move, var_name, rpt_val, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, preskip, postskip, extra_file_name, dir_output, variable, crystal)
	"""
