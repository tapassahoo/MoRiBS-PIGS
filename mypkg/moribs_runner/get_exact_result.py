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
rotor = args.rotor
numb_molecule = args.nmolecule
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
job_submit_dir = os.getcwd()
home = os.path.expanduser("~")
source_code_dir = home + "/" + path_dmrg_dir + "DipoleChain.jl/examples/"
final_result_path = home + "/" + plot_dir_path + "final-dmrg-outputs-for-plotting/"
temp_dir = os.path.dirname(final_result_path)
if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)
#
rfactor = support.get_gfactor_rfactor(rotor, rpt_value, dipole_moment)["R"]
#
support.get_dmrg_result(
	source_code_dir,
	numb_molecule,
	rfactor,
	l_max,
	l_total_max)

rpt_value = "{:3.2f}".format(args.rpt)
file_moved_name = final_result_path + "ground-state-energy-of-" + str(numb_molecule) + rotor + "-at" + str(rpt_value) + "angstrom.txt"
call(["mv", "dmrg_output_for_energy.txt", file_moved_name]) 
