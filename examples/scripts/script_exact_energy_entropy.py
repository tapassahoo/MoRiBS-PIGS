#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support
import inputFile
import math
import argparse

parser = argparse.ArgumentParser(description='It is a Python based script file that is used to compute energy, and von Neumann, R\'{e}nyi entropies of a bipartite many-body system consisting of bipolar rotors pinned into a linear lattice chain. The source codes are written in Julia by Dmitri Iouchtchenko. Module support.py consists of many functions and it is not permitted to modify without consulting the developer - Dr. Tapas Sahoo. User can easily modify module inputFile.py to generate lists of beads (see Getbeads function) as needed. Do not forget to load module julia in your computer or server to execute source .jl files. Before you run it, please visit https://github.com/1/DipoleChain.jl. It is only for rotational DOFs.')
parser.add_argument("-d", "--DipoleMoment", type=float, help="Dipole Moment of a bipolar molecule in Debye.", default = -1.0)
parser.add_argument("-g", "--gFactor", type=float, help="It defines interaction strength.", default = -1.0)
parser.add_argument("-R", "--Rpt", type=float, help="Inter molecular spacing.", default = -1.0)
parser.add_argument("-var", "--var", help="Name of a variable: either beta or tau. It must be a string. Note: for finite temperature computations only the variable tau is needed.", choices =["tau","beta"], default="tau")
parser.add_argument("cal", help="Type of calculation - it is a string: a) PIMC - Finite Temperature calculation by Path Integral Monte Carlo b) PIGS - Ground State Path Integral c) ENT - Entanglement by replica algorithm based on PIGS.", choices = ["PIMC", "PIGS", "ENT"])
parser.add_argument("--scal", help="subtype of calculations - must be defined as a string in case of ENT.", default = "SWAPTOUNSWAP", choices = ["SWAPTOUNSWAP", "BROKENPATH"])
parser.add_argument("-N", help="Number of Molecules. It must be an integer.", type = int)
parser.add_argument("-lmax", "--lmax", help="Maximum l quantum number. It is needed for exact computations of Energy and Entropy of linear rotors in pinned to a chain.", type = int, default = 0)
parser.add_argument("-ltotalmax", "--ltotalmax", help="Maximum lmax quantum number. It is needed for exact computations of Energy and Entropy of linear rotors in pinned to a chain.", type = int, default = 0)
parser.add_argument("--Rotor", help="Name of rotor. E.g. - HF, H2O. It is needed to save rotational density matrix.")
parser.add_argument("--param", type=float, help="Fixed value of beta or tau.", default=-1.0)
args = parser.parse_args()

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
variableName        = args.var
TypeCal             = args.cal
#
molecule_rot        = args.Rotor
molecule            = args.Rotor
numbmolecules       = args.N
#
Rpt             	= args.Rpt
if(args.DipoleMoment):
	dipolemoment	= args.DipoleMoment
if(args.gFactor):
	gfact           = args.gFactor
#
ENT_TYPE 			= args.scal
lmax 				= args.lmax
ltotalmax 			= args.ltotalmax
#
src_dir             = os.getcwd()
extra_file_name     = ""

particleA           = int(numbmolecules/2)

if (variableName == "tau"):
	parameterName   = "beta"
	beta            = args.param
	parameter       = beta
	temperature     = 1.0/beta   

if (variableName == "beta"):
	parameterName   = "tau"
	tau             = args.param
	parameter       = tau

numbblocks	        = 100000
numbpass            = 100
#
srcCodePath         = "/home/tapas/DipoleChain.jl/examples/"
#srcCodePath         = "/home/tapas/DipoleChain.jl-master/examples/"
Units               = support.GetUnitConverter()
BConstant           = support.GetBconst(molecule)  # in wavenumber
BConstantK          = BConstant*Units.CMRECIP2KL
#
user_name           = "tapas"                   
final_results_path  = "/home/"+user_name+"/ResultsOf"+TypeCal+"/"
#
FilePlotName        = support.GetFileNamePlot(TypeCal, molecule_rot, False, True, variableName, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, 0, 0, extra_file_name, final_results_path, particleA, 10)
#
if(args.DipoleMoment > 0.0):
	RFactorList         = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
	RFactor             = RFactorList[0]
if(args.gFactor > 0.0):
	RFactor             = 1.0/math.pow(gfact,1.0/3.0)

support.GetEDResults(TypeCal, FilePlotName, srcCodePath, RFactor, numbmolecules, particleA, lmax, ltotalmax)

#loop = inputFile.Getbeads(TypeCal, variableName)
#support.GetExactValues(FilePlotName, srcCodePath, RFactor, numbmolecules, loop, particleA, molecule_rot, Rpt, dipolemoment, parameter, BConstantK, variableName, TypeCal)
