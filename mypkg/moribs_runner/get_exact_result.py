#!/usr/bin/python

from subprocess import call
import os
from os import system
import decimal
import numpy as np
import argparse
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
	"--rpt",
	type=float,
	metavar='VALUE',
	default=-1.0,
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
method = args.method
rotor_name = args.rotor
numb_rotor = args.nmolecule
#
rpt_value = args.rpt
if (args.dipole_moment):
	dipole_moment = args.dipole_moment
if (args.rfactor):
	rfactor = args.rfactor
#
l_max = args.l_max
l_total_max = args.l_total_max
#
# File systems
root_dir_execution = "scratch"
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
home = os.path.expanduser("~")
source_code_dir = home + "/" + path_dmrg_dir + "DipoleChain.jl/examples/"
final_result_path = home + "/" + plot_dir_path + "final-dmrg-outputs-for-plotting/"
temp_dir = os.path.dirname(final_result_path)
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)

output_file_dir = name_of_output_directory + "/"
user_name = getpass.getuser()
input_dir = os.getcwd() + "/"

if (server_name == "graham"):
	dir_run_job = "/scratch/" + user_name + "/" + output_file_dir
elif (server_name == "nlogn"):
	dir_run_job = "/work/" + user_name + "/" + output_file_dir
else:
	dir_run_job = home + "/" + output_file_dir

temp_dir = os.path.dirname(dir_run_job)
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)

working_file_name = support.get_dmrg_working_file(
	method,
	rotor_name,
	numb_rotor,
	rpt_value,
	dipole_moment,
	l_max,
	l_total_max)
job_execution_dir = working_file_name
#
rfactor = support.get_gfactor_rfactor(rotor_name, rpt_value, dipole_moment)["R"]
#
support.get_dmrg_result(
	server_name,
	root_dir_execution,
	method,
	rotor_name,
	dir_run_job,
	job_execution_dir,
	input_dir,
	source_code_dir,
	numb_rotor,
	rfactor,
	l_max,
	l_total_max)

rpt_value = "{:3.2f}".format(args.rpt)
file_moved_name = final_result_path + "ground-state-energy-of-" + str(numb_molecule) + rotor_name + "-at" + str(rpt_value) + "angstrom.txt"
call(["mv", "dmrg_output_for_energy.txt", file_moved_name]) 
