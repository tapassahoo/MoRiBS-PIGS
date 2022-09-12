#!/usr/bin/python

from subprocess import call
from os import system
import decimal
import numpy as np
import support
import argparse

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
	"--l-max",
	type=int,
	metavar='L',
	help="Local basis truncation.",
	default=0)
parser.add_argument(
	"--l-total-max",
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
print("Dr. Tapas Sahoo")
exit()
#
l_max = args.l-max
l_total_max = args.l-total-max
#
src_dir = os.getcwd()
extra_file_name = ""

particleA = int(numbmolecules / 2)

if (variableName == "tau"):
	parameterName = "beta"
	beta = args.param
	parameter = beta
	temperature = 1.0 / beta

if (variableName == "beta"):
	parameterName = "tau"
	tau = args.param
	parameter = tau

numbblocks = 100000
numbpass = 100
#
srcCodePath = "/home/tapas/DipoleChain.jl/examples/"
#srcCodePath		 = "/home/tapas/DipoleChain.jl-master/examples/"
Units = support.GetUnitConverter()
BConstant = support.GetBconst(molecule)  # in wavenumber
BConstantK = BConstant * Units.CMRECIP2KL
#
user_name = "tapas"
final_results_path = "/home/" + user_name + "/ResultsOf" + TypeCal + "/"
#
FilePlotName = support.GetFileNamePlot(
	TypeCal,
	molecule_rot,
	False,
	True,
	variableName,
	Rpt,
	gfact,
	dipolemoment,
	parameterName,
	parameter,
	numbblocks,
	numbpass,
	numbmolecules,
	molecule,
	ENT_TYPE,
	0,
	0,
	extra_file_name,
	final_results_path,
	particleA,
	10)
#
if (args.dipole_moment > 0.0):
	rfactor_list = support.get_rfactor(rotor, rpt_value, dipole_moment)
	RFactor = RFactorList[0]
if (args.gfactor > 0.0):
	RFactor = 1.0 / math.pow(gfact, 1.0 / 3.0)
exit()

support.GetEDResults(
	TypeCal,
	FilePlotName,
	srcCodePath,
	RFactor,
	numbmolecules,
	particleA,
	lmax,
	ltotalmax)
