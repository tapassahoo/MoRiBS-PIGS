import argparse, decimal, os, sys, subprocess, getpass
from datetime import datetime
import numpy as np
from termcolor import colored

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
parser.add_argument("--get_energy",
					action="store_true",
					help="It provides the expectation value of the total energy of the system.")
parser.add_argument("--get_op",
					action="store_true",
					help="It provides the expectation value of the order parameter of the system.")
parser.add_argument("--get_itcf",
					action="store_true",
					help="It provides the expectation value of the imaginary time correlation function of the system.")
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
get_energy = args.get_energy
get_op = args.get_op
get_itcf = args.get_itcf
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

user_name = getpass.getuser()
input_dir_path = os.getcwd()
home = os.path.expanduser("~")

print("*"*80 + "\n")
print(colored("Developer:".ljust(30),"blue") + colored("Dr. Tapas Sahoo", "yellow") + "\n")
now = datetime.now() # current date and time
date_time = now.strftime("%d/%m/%Y, %H:%M:%S")
print("date and time:".capitalize().ljust(29), date_time, "\n")
print("*"*80 + "\n")

debugging=False
if debugging:
	print(colored("File systems are given below:", "blue") + "\n")
	print("user_name: ".ljust(30) + user_name)
	print("home: ".ljust(30) + home)
	print("input_dir_path: ".ljust(30) + input_dir_path)

source_code_dir_name = os.path.join(path_moribs_dir, "MoRiBS-PIGS")
output_file_dir_name = "name_of_output_directory"
if (method == "PIGS"):
	final_result_path = os.path.join(home, plot_dir_path, "final-pigs-outputs-for-plotting")
elif (method == "PIMC"):
	final_result_path = os.path.join(home, plot_dir_path, "final-pimc-outputs-for-plotting")
else:
	final_result_path = os.path.join(home, plot_dir_path, "final-ent-outputs-for-plotting")
temp_dir = os.path.dirname(final_result_path)
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)

execution_file_path = os.path.join(home, source_code_dir_name)

if (debugging):
	print("source_code_dir_name: ".ljust(30) + source_code_dir_name)
	print("output_file_dir_name: ".ljust(30) + output_file_dir_name)
	print("final_result_path: ".ljust(30) + final_result_path)
	print("execution_file_path: ".ljust(30) + execution_file_path)

if (status == "submission"):
	if (server_name == "graham"):
		run_job_root_dir = os.path.join("/scratch", user_name, output_file_dir_name)
	elif (server_name == "nlogn"):
		run_job_root_dir = os.path.join("/work", user_name, output_file_dir_name)
	else:
		run_job_root_dir = os.path.join(home, output_file_dir_name)

	temp_dir = os.path.dirname(run_job_root_dir)
	if not os.path.exists(temp_dir):
		os.makedirs(temp_dir)

	execution_file = os.path.join(home, source_code_dir_name, "pimc")
	if (debugging):
		print("run_job_root_dir: ".ljust(30) + run_job_root_dir)
		print("source_code_dir_name: ".ljust(30) + source_code_dir_name)
		print("execution_file: ".ljust(30) + execution_file)

	if not args.compiled:
		if not args.restart:
			support.get_execution_file(method, ent_method, execution_file_path)

if (server_name == "graham"):
	output_dir_path = os.path.join("/scratch", user_name, output_file_dir_name)
elif (server_name == "nlogn"):
	output_dir_path = os.path.join("/work", user_name, output_file_dir_name)
else:
	output_dir_path = os.path.join(home, output_file_dir_name)

if (debugging):
	print("output_dir_path: ".ljust(30) + output_dir_path)

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
		slurm_script_dir = os.path.join(run_job_root_dir, working_dir_name + "-Logs")

		if not os.path.isdir(slurm_script_dir):
			subprocess.call(["mkdir", "-p", slurm_script_dir])

		if not args.restart:
			subprocess.call(["cp", execution_file, slurm_script_dir])

		if (rotor_type == "LINEAR"):
			if not args.restart:
				support.compile_rotmat(execution_file_path, input_dir_path)

		if status_cagepot:
			if not args.restart:
				support.compile_cagepot(execution_file_path, input_dir_path)
				support.cagepot(execution_file_path)
				subprocess.call(["mv", "hfc60.pot", slurm_script_dir])

	if (status == "analysis"):
		execution_for="write"
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
			ent_algorithm,
			execution_for)

		if (method != "ENT"):
			if preskip >= numb_block:
				print("*"*80)
				print("Warning!!!!!!!".center(80, " ") + "\n")
				print("Number of Blocks = " + str(numb_block))
				print("Number of preskip= " + str(preskip))
				print(
					"Error message: Number of preskip data must be less than Number of Blocks".center(80, " ") + "\n")
				print("*"*80)
				exit()

			if (get_energy):
				analyzed_energy_file = open(analysis_file_name.save_file_energy, "a")
				analyzed_energy_file.write(support.fmt_energy_data(method, parameter_name))
			if (get_op):
				analyzed_order_parameter_file = open(analysis_file_name.save_file_order_parameter, "a")
				analyzed_order_parameter_file.write(support.fmt_order_parameter(method, parameter_name))
			if (get_itcf):
				analyzed_imaginary_time_correlation_file = open(analysis_file_name.save_file_imaginary_time_correlation, "a")
				analyzed_imaginary_time_correlation_file.write(support.fmt_imaginary_time_correlation(method, parameter_name))

		"""
		if ((status == "analysis") and ((method == "ENT") and (ent_algorithm == "WOR"))):
			analyzed_entropy_file = open(analysis_file_name.SaveEntropy, "a")
			analyzed_entropy_file.write(
				support.fmtAverageEntropy(
					status, parameter_name, ent_method))
		"""

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

			dir_name_trotter_number = working_dir_name + "-Trotter-Number-" + str(ibead)
			print(dir_name_trotter_number)

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

			dir_name_trotter_number = working_dir_name + "-Trotter-Number-" + str(numb_bead)

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
					run_job_root_dir,
					dir_name_trotter_number,
					input_dir_path,
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
					slurm_script_dir,
					cpu_run,
					particle_a,
					partition_name,
					status_cagepot,
					user_name,
					output_file_dir_name,
					execution_file_path,
					restart_bool,
					numb_block_restart,
					crystal,
					rotor_type,
					spin_isomer,
					impurity,
					step_com_impurity,
					level_bisection_impurity)

			if (status == "analysis"):

				final_dir_in_work = os.path.join(output_dir_path, dir_name_trotter_number)
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
						if (get_energy):
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
						if (get_op):
							analyzed_order_parameter_file.write(
								support.get_average_order_parameter(
									debugging,
									method,
									numb_molecule,
									numb_bead,
									parameter_name,
									parameter_value,
									final_dir_in_work,
									preskip,
									postskip,
									numb_block))
						if (get_itcf):
							analyzed_imaginary_time_correlation_file.write(
								support.get_imaginary_time_correlation(
									debugging,
									final_dir_in_work,
									method,
									rotor_type,
									numb_molecule,
									numb_bead,
									parameter_name,
									parameter_value,
									numb_block,
									preskip,
									postskip))
				except BaseException:
					pass

	if (status == "analysis") and (method != "ENT"):
		if (get_energy):
			analyzed_energy_file.close()
			print("\n" + "*"*80 + "\n")
			subprocess.call(["cat", analysis_file_name.save_file_energy])
			support.is_data_exist(analysis_file_name.save_file_energy)
		if (get_op):
			analyzed_order_parameter_file.close()
			print("\n" + "*"*80 + "\n")
			subprocess.call(["cat", analysis_file_name.save_file_order_parameter])
			support.is_data_exist(analysis_file_name.save_file_order_parameter)
		if (get_itcf):
			analyzed_imaginary_time_correlation_file.close()
			print("\n" + "*"*80 + "\n")
			subprocess.call(["cat", analysis_file_name.save_file_imaginary_time_correlation])
			support.is_data_exist(analysis_file_name.save_file_imaginary_time_correlation)
		print("\n" + "*"*80 + "\n")

	if ((status == "analysis") and ((method == "ENT") and (ent_algorithm == "WOR"))):
		analyzed_entropy_file.close()
		subprocess.call(["cat", analysis_file_name.SaveEntropy])
		print("\n\n")
		support.is_data_exist(analysis_file_name.SaveEntropy)

"""
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
"""
