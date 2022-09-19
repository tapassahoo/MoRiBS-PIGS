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
					choices=["submission", "analysis"],
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
root_dir_run = "scratch"  # work
cpu_run = "nCPU"
extra_file_name = extra_name
#
myhost = os.uname()[1]
myhost = myhost[0:3]
if ((myhost == "gra") or (myhost == "ced")):
	server_name = "graham"
elif ((myhost == "feynman") or (myhost == "nlogn")):
	server_name = "nlogn"
else:
	server_name = "moribs"
#
job_submit_dir = os.getcwd()

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

source_code_dir = path_moribs_dir+"MoRiBS-PIGS/"
output_file_dir = "name_of_output_directory/"
if (method == "PIGS"):
	final_result_path = home + "/" + plot_dir_path + "final-pigs-outputs-for-plotting/"
elif (method == "PIMC"):
	final_result_path = home + "/" + plot_dir_path + "final-pimc-outputs-for-plotting/"
else:
	final_result_path = home + "/" + plot_dir_path + "final-ent-outputs-for-plotting/"
temp_dir = os.path.dirname(final_result_path)
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)

source_dir_exe = home + "/" + source_code_dir

if (status == "submission"):
	if (server_name == "graham"):
		dir_run_job = "/scratch/" + user_name + "/" + output_file_dir
	elif (server_name == "nlogn"):
		dir_run_job = "/work/" + user_name + "/" + output_file_dir
	else:
		dir_run_job = home + "/" + output_file_dir

	temp_dir = os.path.dirname(dir_run_job)
	if not os.path.exists(temp_dir):
		os.makedirs(temp_dir)

	execution_file = home + "/" + source_code_dir + "pimc"
	if not args.compiled:
		if not args.restart:
			support.get_execution_file(
				job_submit_dir, method, ent_method, source_dir_exe)

if (server_name == "graham"):
	output_dir_path = "/scratch/" + user_name + "/" + output_file_dir
elif (server_name == "nlogn"):
	output_dir_path = "/work/" + user_name + "/" + output_file_dir
else:
	output_dir_path = home + "/" + output_file_dir


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

	if (status == "submission"):
		if (server_name == "graham"):
			dir_run_input_pimc = "/scratch/" + user_name + \
				"/" + output_file_dir + working_dir_name + "-Logs"
		elif (server_name == "nlogn"):
			dir_run_input_pimc = "/work/" + user_name + "/" + output_file_dir + working_dir_name + "-Logs"
		else:
			dir_run_input_pimc = home + "/" + output_file_dir + working_dir_name + "-Logs"

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

	if (status == "analysis"):
		analysis_file_name = support.GetAnalysisFileName(
			final_result_path,
			method,
			False,
			molecular_system,
			rotor,
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
			preskip,
			postskip,
			extra_file_name,
			particle_a,
			ent_method,
			ent_algorithm)

		if (method != "ENT"):
			if preskip >= numb_block:
				print("")
				print("Warning!!!!!!!")
				print(
					"============================================================================")
				print("Number of Blocks = " + str(numb_block))
				print("Number of preskip= " + str(preskip))
				print(
					"Error message: Number of preskip data must be less than Number of Blocks")
				print(
					"============================================================================")
				exit()

			analyzed_energy_file = open(analysis_file_name.save_file_energy, "a")
			analyzed_energy_file.write(
				support.fmt_energy_data(
					method, parameter_name))
			analyzed_correlation_file = open(analysis_file_name.save_file_correlation, "a")
			analyzed_correlation_file.write(support.fmtAverageOrderParam(status, parameter_name))

		'''
		if ((status == "analysis") and ((method == "ENT") and (ent_algorithm == "WOR"))):
			analyzed_entropy_file = open(analysis_file_name.SaveEntropy, "a")
			analyzed_entropy_file.write(
				support.fmtAverageEntropy(
					status, parameter_name, ent_method))
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
					root_dir_run,
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
					analyzed_energy_file.write(
						support.get_average_energy(
							method,
							numbbeads,
							variable,
							final_dir_in_work,
							preskip,
							postskip,
							numbblocks))
					analyzed_correlation_file.write(
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
				beta = parameter_value * (numb_bead - 1)
				temperature = 1.0 / beta
				variable = beta
			if (parameter_name == "beta"):
				beta = parameter_value
				temperature = 1.0 / beta
				tau = beta / (numb_bead - 1)
				variable = tau

			execution_bead_dir_name = working_dir_name + "-Trotter-Number-" + str(numb_bead)

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
					root_dir_run,
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

			if (status == "analysis"):

				final_dir_in_work = output_dir_path + execution_bead_dir_name
				#support.RemoveFiles(method, numbbeads, temperature, rotor, RotorType, preskip, postskip, numbblocks, final_dir_in_work)
				try:
					if ((method == "ENT") and (ent_algorithm == "WOR")):
						analyzed_entropy_file.write(
							support.GetAverageEntropy(
								numb_bead,
								parameter_value,
								final_dir_in_work,
								preskip,
								postskip,
								numb_block,
								ent_method))
					if (method != "ENT"):
						analyzed_energy_file.write(
							support.get_average_energy(
								method,
								numb_bead,
								parameter_name,
								parameter_value,
								final_dir_in_work,
								preskip,
								postskip,
								numb_block))
						'''
						analyzed_correlation_file.write(
							support.GetAverageOrderParam(
								method,
								numb_molecule,
								numb_bead,
								parameter_name,
								parameter_value,
								final_dir_in_work,
								preskip,
								postskip,
								numb_block))
						'''
				except BaseException:
					pass

	if (status == "analysis") and (method != "ENT"):
		analyzed_energy_file.close()
		analyzed_correlation_file.close()
		call(["cat", analysis_file_name.save_file_energy])
		print("")
		print("")
		call(["cat", analysis_file_name.save_file_correlation])
		print("")
		print("")
		# =========================File Checking===============================#
		SavedFile = analysis_file_name.save_file_energy
		support.FileCheck(method, bead_list, parameter_name, SavedFile)
		SavedFile = analysis_file_name.save_file_correlation
		support.FileCheck(method, bead_list, parameter_name, SavedFile)

	if ((status == "analysis") and ((method == "ENT") and (ent_algorithm == "WOR"))):
		analyzed_entropy_file.close()
		call(["cat", analysis_file_name.SaveEntropy])
		print("")
		print("")
		SavedFile = analysis_file_name.SaveEntropy
		support.FileCheck(method, bead_list, parameter_name, SavedFile)

'''
if ((status == "analysis") and ((method == "ENT") and (
		ent_method == "EXTENDED_ENSMBL") and (ent_algorithm == "WR"))):
	print("Final Entropy obtained by employing Ratio Trick")
	support.GetAverageEntropyRT(
		particleAList,
		method,
		rotor,
		translational_move,
		rotational_move,
		parameter_name,
		rpt_val,
		gfactor,
		dipolemoment,
		parameterName,
		parameter,
		numbblocks,
		numbpass,
		numbmolecules1,
		molecule,
		ent_method,
		preskip,
		postskip,
		extra_file_name,
		output_dir_path,
		variable,
		crystal,
		final_result_path,
		impurity,
		ext_ent)
	support.GetEntropyRT(status, maxloop, method, rotor, translational_move, rotational_move, parameter_name, rpt_val, gfactor, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ent_method, preskip, postskip, extra_file_name, output_dir_path, variable, crystal)
'''
