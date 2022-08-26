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
import mypkg.moribs_runner.support as support
import get_beads_and_mc_steps as mc

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
parser.add_argument("--ent_algorithm",
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
					choices=["LINEAR", "NONLINEAR"],
					help="Users should specify the type of the rotor \
					either LINEAR or NONLINEAR.")
parser.add_argument("--spin_isomer",
					type=int,
					default=-1,
					choices=[-1, 0, 1],
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

# Initial setup
status = args.job
method = args.method
molecular_system = args.system
translational_move = args.com_move
rotational_move = args.rot_move
rotor = args.rotor
rotor_type = args.rotor_type
spin_isomer = args.spin_isomer
numb_molecule1 = args.nmolecule
parameter_name = args.parameter_name
parameter_value = args.parameter_value
numb_block = args.nblock
numb_pass = args.npass
#
rpt_val = args.rpt
dipole_moment = args.dipole_moment
crystal = args.crystal
ent_method = args.ent_method
ent_algorithm = args.ent_algorithm
gfactor = args.gfactor
#
if (args.impurity):
    impurity=args.impurity
else:
    impurity=[""]
#
partition_name = args.partition
preskip = args.preskip
postskip = args.postskip
#
# Request to change
# User should change the following 4 lines as developer already have
# explained in the README file. If user wish to include cage potential,
# "No" flag must be turned into "Yes" by the user.

status_cagepot = False
dir_run = "scratch"  # work
cpu_run = "nCPU"
extra_file_name = extra_name
#
myhost = os.uname()[1]
myhost = myhost[0:3]
if ((myhost == "gra") or (myhost == "ced")):
	server_name = "graham"
else:
	server_name = "nlogn"
#
submit_job_dir = os.getcwd()

if (spin_isomer == -1):
	prefix_name = ""
if (spin_isomer == 0):
	prefix_name = "-p-"
if (spin_isomer == 1):
	prefix_name = "-o-"
rotor_name = prefix_name + rotor

mc_step = mc.GetBeadStepLevel(molecular_system, parameter_name, method)
bead_list = mc_step.beads
step_com_move = mc_step.step_com
level_bisection = mc_step.level_com
step_rot_move = mc_step.step_rot
#
step_com_impurity = step_com_move
level_bisection_impurity = level_bisection

numb_block_restart = args.nblock_restart

# User should change the following 5 lines as developer already have
# explained in the README file.
user_name = getpass.getuser()
input_dir = os.getcwd() + "/"
home = os.path.expanduser("~")
source_code_dir = "MoRiBS-PIGS/"
output_file_dir = "name_of_output_directory/"
final_results_path = home + "/results-of-" + method + "/"
temp_dir = os.path.dirname(final_results_path)
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)

source_dir_exe = "/home/" + user_name + "/" + source_code_dir

if (status == "submission"):
	if ((dir_run == "scratch") or (server_name == "graham")):
		dir_run_job = "/scratch/" + user_name + "/" + output_file_dir
		if (server_name == "graham"):
			temp_dir = os.path.dirname(dir_run_job)
			if not os.path.exists(temp_dir):
				os.makedirs(temp_dir)
	else:
		dir_run_job = "/work/" + user_name + "/" + output_file_dir

	execution_file = "/home/" + user_name + "/" + source_code_dir + "pimc"
	if not args.compiled:
		if not args.restart:
			support.makeexecutionfile(
				submit_job_dir, method, ent_method, source_dir_exe)

if (server_name == "graham"):
	output_dir_path = "/scratch/" + user_name + "/" + output_file_dir
else:
	output_dir_path = "/work/" + user_name + "/" + output_file_dir

ent_algorithm = args.ent_algorithm
if (method == "ENT"):
	maxloop = int(numb_molecule1 / 2)
else:
	maxloop = 1

if (args.ent_algorithm == "WR"):
	particle_a_list = np.arange(1, maxloop + 1)
else:
	particle_a_list = [maxloop]

for particle_a in particle_a_list:
	working_dir_name = support.get_working_file(
		method,
		molecular_system,
		numb_molecule1,
		translational_move,
		rotational_move,
		parameter_name,
		parameter_value,
		rpt_val,
		dipole_moment,
		gfactor,
		numb_block,
		numb_pass,
		rotor_name,
		extra_file_name,
		crystal,
		impurity,
		ent_method,
		particle_a,
		ent_algorithm)

	'''
	if (status == "rename"):
		numbblocks_rename = args.NR
		file2_name = support.get_woring_file(
			method,
			rotor,
			translational_move,
			rotational_move,
			rpt_val,
			gfact,
			dipolemoment,
			parameterName,
			parameter,
			numbblocks_rename,
			numbpass,
			numbmolecules1,
			molecule,
			ENT_TYPE,
			particleA,
			extra_file_name,
			crystal,
			impurity,
			ext_ent)

		if (server_name == "graham"):
			dir_run_input_pimc = "/scratch/" + user_name + \
				"/" + output_file_dir + working_dir_name + "-Logs"
			dir_input_pimc_renamed = "/scratch/" + \
				user_name + "/" + output_file_dir + file2_name + "-Logs"
		else:
			dir_run_input_pimc = "/work/" + user_name + "/" + output_file_dir + working_dir_name + "-Logs"
			dir_input_pimc_renamed = "/work/" + user_name + \
				"/" + output_file_dir + file2_name + "-Logs"

	'''
	# ===============================================================================
	#																			   |
	#	compilation of linden.f to generate rotational density matrix - linden.out |
	#	Yet to be generalized													   |
	#																			   |
	# ===============================================================================
	if (status == "submission"):
		if (server_name == "graham"):
			dir_run_input_pimc = "/scratch/" + user_name + \
				"/" + output_file_dir + working_dir_name + "-Logs"
		else:
			dir_run_input_pimc = "/work/" + user_name + "/" + output_file_dir + working_dir_name + "-Logs"

		if not os.path.isdir(dir_run_input_pimc):
			call(["mkdir", "-p", dir_run_input_pimc])

		if not args.restart:
			call(["cp", execution_file, dir_run_input_pimc])

		if (rotor_type == "LINEAR"):
			if not args.restart:
				support.compile_rotmat(source_dir_exe, input_dir)

		if status_cagepot:
			if not args.restart:
				support.compile_cagepot(source_dir_exe, input_dir)
				support.cagepot(source_dir_exe)
				call(["mv", "hfc60.pot", dir_run_input_pimc])

	'''
	if (status == "analysis"):
		FileAnalysis = support.GetFileNameAnalysis(
			method,
			False,
			rotor,
			translational_move,
			rotational_move,
			var_name,
			rpt_val,
			gfact,
			dipolemoment,
			parameterName,
			parameter,
			numbblocks,
			numbpass,
			numbmolecules1,
			molecule,
			ENT_TYPE,
			preskip,
			postskip,
			extra_file_name,
			final_results_path,
			particleA,
			ENT_ALGR)

		if (method != "ENT"):
			if preskip >= numbblocks:
				print("")
				print("Warning!!!!!!!")
				print(
					"============================================================================")
				print("Number of Blocks = " + str(numbblocks))
				print("Number of preskip= " + str(preskip))
				print(
					"Error message: Number of preskip data must be less than Number of Blocks")
				print(
					"============================================================================")
				exit()

			fanalyzeEnergy = open(FileAnalysis.SaveEnergy, "a")
			fanalyzeEnergy.write(
				support.fmtAverageEnergy(
					method, status, var_name))
			fanalyzeCorr = open(FileAnalysis.SaveCorr, "a")
			fanalyzeCorr.write(support.fmtAverageOrderParam(status, var_name))

		if ((status == "analysis") and ((method == "ENT") and (ENT_ALGR == "WOR"))):
			fanalyzeEntropy = open(FileAnalysis.SaveEntropy, "a")
			fanalyzeEntropy.write(
				support.fmtAverageEntropy(
					status, var_name, ENT_TYPE))
	'''

	if (method == "ENT"):
		numb_molecule = 2 * numb_molecule1
	else:
		numb_molecule = numb_molecule1

	for index, ibead in enumerate(bead_list, start=0):

		if (method == "PIMC"):

			if (parameter_name == "tau"):
				beta = tau * value
				temperature = 1.0 / beta
				variable = beta
			if (parameter_name == "beta"):
				tau = beta / value
				variable = tau

			execution_bead_dir_name = working_dir_name + "-Trotter-Number-" + str(ibead)
			print(execution_bead_dir_name)

			'''
			if (status == "submission"):

				restart_bool = False
				if args.restart:
					restart_bool = True

				support.job_submission(
					server_name,
					status,
					translational_move,
					rotational_move,
					dir_run,
					dir_run_job,
					execution_bead_dir_name,
					input_dir,
					execution_file,
					rpt_val,
					bead_list,
					ibead,
					step_rot,
					step_COM,
					level_bisection,
					temperature,
					numb_block,
					numb_pass,
					rotor,
					numb_molecule,
					gfactor,
					dipole_moment,
					method,
					ent_method,
					ent_algorithm,
					output_dir_path,
					dir_run_input_pimc,
					RUNIN,
					particle_a,
					partition_name,
					status_cagepot,
					index,
					user_name,
					output_file_dir,
					source_dir_exe,
					restart_bool,
					numbblocks_Restart1,
					crystal,
					rotor_type,
					spin_isomer,
					impurity,
					step_COM_impurity,
					level_bisection_impurity)

			if (status == "analysis"):

				final_dir_in_work = output_dir_path + execution_bead_dir_name
				try:
					fanalyzeEnergy.write(
						support.GetAverageEnergy(
							method,
							numbbeads,
							variable,
							final_dir_in_work,
							preskip,
							postskip,
							numbblocks))
					fanalyzeCorr.write(
						support.GetAverageOrderParam(
							method,
							numbmolecules,
							numbbeads,
							variable,
							final_dir_in_work,
							preskip,
							postskip,
							numbblocks))
				except BaseException:
					pass
			'''
		else:

			if (ibead % 2) != 0:
				numb_bead = ibead
			else:
				numb_bead = ibead + 1

			if (parameter_name == "tau"):
				beta = tau * (numb_bead - 1)
				temperature = 1.0 / beta
				variable = beta
			if (parameter_name == "beta"):
				beta = parameter_value
				temperature = 1.0 / beta
				tau = beta / (numb_bead - 1)
				variable = tau

			execution_bead_dir_name = working_dir_name + "-Trotter-Number-" + str(numb_bead)

			'''
			if (status == "rename"):
				folder_rename = file2_name + str(numbbeads)
				support.GetRenamingFunc(
					dir_run_input_pimc,
					dir_input_pimc_renamed,
					output_dir_path,
					execution_bead_dir_name,
					folder_rename,
					input_dir)
			'''

			if (status == "submission"):

				if args.restart:
					restart_bool = True
				else:
					restart_bool = False

				support.job_submission(
					server_name,
					status,
					translational_move,
					rotational_move,
					dir_run,
					dir_run_job,
					execution_bead_dir_name,
					input_dir,
					execution_file,
					rpt_val,
					ibead,
					numb_bead,
					step_rot_move,
					step_com_move,
					level_bisection,
					temperature,
					numb_block,
					numb_pass,
					rotor,
					numb_molecule,
					gfactor,
					dipole_moment,
					method,
					ent_method,
					ent_algorithm,
					output_dir_path,
					dir_run_input_pimc,
					cpu_run,
					particle_a,
					partition_name,
					status_cagepot,
					user_name,
					output_file_dir,
					source_dir_exe,
					restart_bool,
					numb_block_restart,
					crystal,
					rotor_type,
					spin_isomer,
					impurity,
					step_com_impurity,
					level_bisection_impurity)

			'''
			if status == "analysis":

				final_dir_in_work = output_dir_path + execution_bead_dir_name
				#support.RemoveFiles(method, numbbeads, temperature, rotor, RotorType, preskip, postskip, numbblocks, final_dir_in_work)
				try:
					if ((method == "ENT") and (ENT_ALGR == "WOR")):
						fanalyzeEntropy.write(
							support.GetAverageEntropy(
								numbbeads,
								variable,
								final_dir_in_work,
								preskip,
								postskip,
								numbblocks,
								ENT_TYPE))
					if (method != "ENT"):
						fanalyzeEnergy.write(
							support.GetAverageEnergy(
								method,
								numbbeads,
								variable,
								final_dir_in_work,
								preskip,
								postskip,
								numbblocks))
						fanalyzeCorr.write(
							support.GetAverageOrderParam(
								method,
								numbmolecules,
								numbbeads,
								variable,
								final_dir_in_work,
								preskip,
								postskip,
								numbblocks))
				except BaseException:
					pass
			'''

	if (status == "analysis") and (method != "ENT"):
		fanalyzeEnergy.close()
		fanalyzeCorr.close()
		call(["cat", FileAnalysis.SaveEnergy])
		print("")
		print("")
		call(["cat", FileAnalysis.SaveCorr])
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

'''
if ((status == "analysis") and ((method == "ENT") and (
		ENT_TYPE == "EXTENDED_ENSMBL") and (ENT_ALGR == "WR"))):
	print("Final Entropy obtained by employing Ratio Trick")
	support.GetAverageEntropyRT(
		particleAList,
		method,
		rotor,
		translational_move,
		rotational_move,
		var_name,
		rpt_val,
		gfact,
		dipolemoment,
		parameterName,
		parameter,
		numbblocks,
		numbpass,
		numbmolecules1,
		molecule,
		ENT_TYPE,
		preskip,
		postskip,
		extra_file_name,
		output_dir_path,
		variable,
		crystal,
		final_results_path,
		impurity,
		ext_ent)
	support.GetEntropyRT(status, maxloop, method, rotor, translational_move, rotational_move, var_name, rpt_val, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, preskip, postskip, extra_file_name, output_dir_path, variable, crystal)
'''
