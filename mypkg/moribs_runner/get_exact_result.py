#!/usr/bin/python

from subprocess import call
import os
from os import system
import decimal
import numpy as np
import argparse, json
import getpass
import mypkg.moribs_runner.support as support

parser = argparse.ArgumentParser(
	description="The script is used to calculate the ground-state energy, \
				von Neumann entropy, and Rényi entropy of a bipartite many-body \
				system consisting of bipolar rotors pinned into a linear lattice \
				chain. The source code was written in Julia by Dmitri Iouchtchenko \
				from the University of Waterloo, Canada. Module support.py contains \
				many functions, so users must consult the developer before modifying \
				it. Please visit the GitHub page - https://github.com/tapassahoo/DipoleChain \
				before you carry out the computations.",
	epilog="Enjoy the program! :)")
parser.add_argument("job",
	type=str,
	choices=["submission", "analysis"],
	help="Type of jobs: submission of new jobs or analyzing \
	output files.")
parser.add_argument(
	"method",
	type=str,
	choices=["dmrg","mm"],
	help="Either dmrg or piqmc.")
parser.add_argument(
	"-d",
	"--dipole_moment",
	type=float,
	metavar='VALUE',
	default=-1.0,
	help="It defines the dipole moment of a linear \
	polar molecule in Debye. It is applicable only \
	for the polar linear rotors.")
parser.add_argument(
	"--rfactor",
	type=float,
	metavar='VALUE',
	default=-1.0,
	help="It defines interaction strength. \
	It is applicable only for the polar linear rotors.")
parser.add_argument(
	"--r_list",
	type=float,
	nargs='+',
	help="Distance between Centre of Masses of two \
	molecules. The unit is Angstrom.")
parser.add_argument(
	"-n",
	"--nmolecule",
	type=int,
	metavar='N',
	help="Number of Molecules.")
parser.add_argument(
	"--l_max",
	type=int,
	metavar='L',
	help="Local basis truncation.",
	default=0)
parser.add_argument(
	"--l_total_max",
	type=int,
	metavar='L',
	help="Many-body basis truncation (default: no truncation)",
	default=0)
parser.add_argument(
	"--rotor",
	type=str,
	metavar="NAME",
	help="Name of the rotor. E.g., HF, H2O.")
args = parser.parse_args()
#
status=args.job
method = args.method
rotor_name = args.rotor
numb_rotor = args.nmolecule
#
rpt_list = args.r_list
if (args.dipole_moment):
	dipole_moment = args.dipole_moment
if (args.rfactor):
	rfactor = args.rfactor
#
l_max = args.l_max
l_total_max = args.l_total_max
#
# File systems
myhost = os.uname()[1]
myhost = myhost[0:3]
if ((myhost == "gra") or (myhost == "ced")):
	server_name = "graham"
	root_dir_execution = "scratch"
elif ((myhost == "feynman") or (myhost == "nlogn")):
	server_name = "nlogn"
	root_dir_execution = "work"
else:
	server_name = "moribs"
	root_dir_execution = "home"
#
job_submit_dir = os.getcwd()
home = os.path.expanduser("~")
source_code_dir = home + "/" + path_dmrg_dir + "DipoleChain.jl/examples/"
final_output_dir = home + "/" + plot_dir_path + "final-dmrg-outputs-for-plotting/"
temp_dir = os.path.dirname(final_output_dir)
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)

output_file_dir = name_of_output_directory + "/"
user_name = getpass.getuser()
input_dir = os.getcwd() + "/"

if (server_name == "graham"):
	job_submission_root_dir = "/scratch/" + user_name + "/" + output_file_dir
elif (server_name == "nlogn"):
	job_submission_root_dir = "/work/" + user_name + "/" + output_file_dir
else:
	job_submission_root_dir = home + "/" + output_file_dir

temp_dir = os.path.dirname(job_submission_root_dir)
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)

if (status == "analysis"):
	final_result_dir = home + "/final-dmrg-outputs-for-plotting"
	final_output_file = final_result_dir + "/dmrg-results-of-" + str(numb_rotor) + rotor_name + "-for-energy-von-neumann-and-second-renyi-entropies-vs-intermolecular-distance-lmax" + str(l_max) + "-ltotalmax" + str(l_total_max) + ".txt"

for rpt_value in rpt_list:
	rpt_value = "{:3.2f}".format(rpt_value)
	working_file_name = support.get_dmrg_working_file(
		method,
		rotor_name,
		numb_rotor,
		float(rpt_value),
		dipole_moment,
		l_max,
		l_total_max)
	job_execution_dir = working_file_name

	if (status == "submission"):
		rfactor = support.get_gfactor_rfactor(rotor_name, float(rpt_value), dipole_moment)["R"]

		support.get_dmrg_result(
			server_name,
			root_dir_execution,
			method,
			rotor_name,
			job_submission_root_dir,
			job_execution_dir,
			input_dir,
			source_code_dir,
			numb_rotor,
			rfactor,
			l_max,
			l_total_max,
			final_output_dir)

	if (status == "analysis"):
		dmrg_output_dir = job_submission_root_dir + working_file_name
		data_file = dmrg_output_dir + "/dmrg_output_for_energy.txt"
		get_data = np.genfromtxt(data_file)
		print(get_data)
