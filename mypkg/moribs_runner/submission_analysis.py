import argparse
import decimal
import os
import sys
import time
from subprocess import call
import getpass
import datetime
import numpy as np

sys.path.append("../../examples/scripts")
import get_beads_and_mc_steps as mc
import mypkg.moribs_runner.support as support

parser = argparse.ArgumentParser(
					description="It is used to submit jobs in a queue \
					and analyze output data files. \
					Note: /MoRiBS-PIGS/mypkg/moribs_runner/support.py \
					consists of multiple functions and is not permitted \
					to modify without consulting the developer - \
					Dr Tapas Sahoo. Users should modify module \
					/MoRiBS-PIGS/examples/scripts/get_beads_and_mc_steps.py \
					to generate lists of beads, step lengths for rotational \
					and translational motions, and levels for bisection move \
					associated with each bead.",
					epilog="Enjoy the program! :)")
parser.add_argument("job",
				    type=str,
					choices=["submission", "analysis", "rename"], 
					help="Type of jobs: submission of new jobs or analyzing \
					output files.")
parser.add_argument("method", 
					type=str,
					choices=["PIMC", "PIGS", "ENT"], 
					help="The methodology should be specified as \
					one of PIMC, PIGS, and ENT. PIMC is for finite \
					temperature calculation; PIGS is for ground state \
					calculation; ENT is for entanglement by replica \
					algorithm based on PIGS.")
parser.add_argument("--compiled", 
					action="store_true", 
					help="Users can use it if the execution file \
					(/MoRiBS-PIGS/pimc) is already generated.")
parser.add_argument("--ent_method", 
					type=str,
					default="EXTENDED_ENSMBL", 
					choices=["EXTENDED_ENSMBL", "BROKENPATH"],
					help="The subtype of calculations - \
					must be defined in the case of ENT.")
parser.add_argument("--algorithm", 
					type=str,
					default="WR", 
					choices=["WR", "WOR"],
				    help="The subtype of calculations - \
					must be defined in the case of ENT. \
					WR and WOR stand for with and without \
					ratio trick algorithm.")
parser.add_argument("system", 
					type=str,
					help="Name of a molecular system. \
					E.g., H2, HF, H2O, FCC-H2O, H2O@c60")
parser.add_argument("parameter_name", 
					type=str, 
					help="Name of the parameter - either beta or tau.")
parser.add_argument("parameter_value", 
					type=float,
					help="The value associated with the parameter. \
					The unit is inverse temperature.")
parser.add_argument("--rotor", 
					type=str,
					metavar="NAME",
					help="Name of the rotor. E.g., HF, H2O. \
					It is needed to call the rotational density matrix file.")
parser.add_argument("--com_move", 
					action="store_true", 
					help="It allows translational motions of \
					molecules or particles.")
parser.add_argument("--rot_move", 
					action="store_true", 
					help="It allows rotational motions of molecules.")
parser.add_argument("--rotor_type", 
					type=str,
					default="LINEAR", 
					choices=["LINEAR","NONLINEAR"],
					help="Users should specify the type of the rotor \
					either LINEAR or NONLINEAR.")
parser.add_argument("--spin_isomer", 
					type=int, 
					default=-1, 
					choices=[-1,0,1],
					help="It has three values. These are -1, 0, and 1 \
					for spinless, para and ortho isomers, respectively. \
					Users should mention one of them.")
parser.add_argument("-n", 
					"--nmolecule",
					type=int,
					metavar='NUMBER',
					help="Number of Molecules.")
parser.add_argument("--nblock", 
					type=int,
					metavar='NUMBER',
					help="Number of Blocks required for Monte Carlo move.")
parser.add_argument("--npass", 
					type=int,
					metavar='NUMBER',
					help="Number of Passes required for Monte Carlo move.")
parser.add_argument("-r", 
					"--rpt", 
					type=float, 
					metavar='VALUE',
					default=-1.0,
					help="Distance between Centre of Masses of two \
					molecules. The unit is Angstrom.")
parser.add_argument("-d", 
					"--dipole_moment", 
					type=float, 
					metavar='VALUE',
					default=-1.0,
					help="It defines the dipole moment of a linear \
					polar molecule in Debye. It is applicable only \
					for the polar linear rotors.")
parser.add_argument("-g", 
					"--gfactor", 
					type=float, 
					metavar='VALUE',
					default=-1.0,
					help="It defines interaction strength. \
					It is applicable only for the polar linear rotors.")
parser.add_argument("--restart", 
					action="store_true", 
					help="It is used to restart the code.")
parser.add_argument("--nblock_restart", 
					type=int, 
					metavar='NUMBER',
					default=0,
					help="The number of Monte Carlo blocks users \
					wish to extend further is mentioned.")
parser.add_argument("--preskip", 
					metavar='NUMBER',
					type=int, 
					default=0,
					help="The number of blocks in outputs \
					skipped from the beginning. It can be necessary \
					while the analysis flag is open to remove \
					pre-equilibrated data.")
parser.add_argument("--postskip", 
					type=int, 
					metavar='NUMBER',
					default=0,
					help="The number of blocks in outputs skipped \
					from the end of the file. It can be necessary \
					to increase the error bar. Ensure the analysis \
					flag is open.")
parser.add_argument("-im", 
					"--impurity", 
					action="store", 
					nargs=4, 
					help="Use it if the system comprises various \
					types of particles and molecules. \
					IMPURITY[0] = type of particle (ATOM, MOLECULE, \
					PLANAR, LINEAR or NONLINEAR); \
					IMPURITY[1] = Name of the particle; \
					IMPURITY[2] = Number of the particle; \
					IMPURITY[3] = Either BOSE or BOLTZMANN.")
parser.add_argument("--partition", 
					type=str,
					default="ntapas",
					metavar='CPU',
					help="It allows submitting jobs in a specific CPU. \
					The user does not need it.")
parser.add_argument("--crystal", 
					action="store_true", 
					help="It is required to read lattice configurations \
					and the associated dipole moments. \
					It is in developing condition.")
args = parser.parse_args()

status = args.job
method = args.method
parameter_name = args.parameter_name
tanslational_move = args.com_move
rotational_move = args.rot_move
if (args.impurity):
	impurity=args.impurity
else:
	impurity=[""]

print(status)
exit()
#
myhost = os.uname()[1]
myhost = myhost[0:3]
if ((myhost == "gra") or (myhost == "ced")):
	server_name = "graham"
else:
	server_name = "nlogn"
partition_name = args.partition
#
#
molecular_system = args.system
rotor = args.rotor
#
numb_blocks = args.nblock
numb_molecules1 = args.nmolecule
numb_passes = args.npass
#
rpt_val = args.rpt
dipole_moment = args.dipole_moment
gfactor = args.gfactor
#support.GetrAndgFactor(molecule_rot, rpt_val, dipolemoment)

# Request to change
# User should change the following 4 lines as developer already have explained in the README file. If user wish to include cage potential, "No" flag must be turned into "Yes" by the user.

status_cagepot = False
run_dir="scratch" # work
run_cpu="nCPU"

preskip = args.preskip
postskip = args.postskip
crystal = args.crystal
rotor_type = args.rotor_type

ent_method = args.ent_method
ent_algorithm = args.algorithm

extra_file_name = extra_name

src_dir = os.getcwd()
if (variable_name == "tau"):
	parameter_name = "beta"
	beta = args.parameter
	parameter = beta
	temperature = 1.0/beta

if (variable_name == "beta"):
	parameter_name = "tau"
	tau = args.parameter
	parameter = tau

spin_isomer = args.spin_isomer

if (spin_isomer == -1):
	prefix_name = ""
if (spin_isomer == 0):
	prefix_name = "-p-"
if (spin_isomer == 1):
	prefix_name = "-o-"
rotor_name = prefix_name+rotor

print(rotor_name)
exit()

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
out_dir = "name_of_output_directory/"
final_results_path = home + "/results-of-" + method + "/"
dir = os.path.dirname(final_results_path)
if not os.path.exists(dir):
	os.makedirs(dir)

source_dir_exe = "/home/" + user_name + "/" + source_dir
if (status == "submission"):

	if ((RUNDIR == "scratch") or (server_name == "graham")):
		dir_run_job = "/scratch/" + user_name + "/" + out_dir
		if (server_name == "graham"):
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
