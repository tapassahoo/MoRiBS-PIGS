#!/usr/bin/python
 
import time
import os
from os import system
from subprocess import call
import numpy as np
from numpy import *
import decimal
import support
import FigureGenerator

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
variableName        = "tau"
#variableName        = "beta"
#
TransMove           = False
RotMove             = True
#
#TypeCal             = 'PIMC'
#TypeCal             = 'PIGS'
TypeCal             = 'ENT'
#
#TypePlot            = "Energy"
#TypePlot            = "ChemPot"
#TypePlot            = "CorrFunc"
#TypePlot            = "S2"
TypePlot            = "GFACTOR"
#TypePlot            = "COMBINE"
#
#molecule            = "HFC60"                                                  
molecule            = "HF"                                                      
#molecule            = "H2"                                                    
molecule_rot        = "HF"
#
Rpt                 = 10.05
dipolemoment        = 1.826      #J. Chem. Phys. 73(5), 2319 (1980).
dipolemoment        = 1.0*dipolemoment

ENT_TYPE 			= "SWAPTOUNSWAP"
#ENT_TYPE 			= "SWAP"
#ENT_TYPE 			= "BROKENPATH"
#ENT_TYPE 			= "REGULARPATH"

extra_file_name     = ""

user_name           = "tapas"
final_results_path  = "/home/"+user_name+"/ResultsOf"+TypeCal+"/"
if (variableName == "tau"):
	parameterName   = "beta"
	beta            = 0.2
	parameter       = beta
	temperature     = 1.0/beta   

if (variableName == "beta"):
	parameterName   = "tau"
	tau             = 0.005
	parameter       = tau

#==================================Plotting====================================#
if (TypeCal == "ENT" and TypePlot == "GFACTOR"):
	beadsRef            = 41
	numbblocks	        = 5000
	numbpass            = 400
	preskip             = 0
	postskip            = 0

	extra_file_name     = ""
	FigureGenerator.FigureEntropyRT(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, parameterName, parameter, numbblocks, numbpass, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, TypePlot, beadsRef)

'''
if (TypeCal == "ENT" and TypePlot == "GFACTOR" or TypePlot == "S2"):
	beadsRef = 61
	numbblocks	        = 50000
	numbmolecules       = 2
	numbpass            = 200
	particleA           = int(numbmolecules/2)
	preskip             = 1000
	postskip            = 0

	FigureGenerator.FigureENT(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA,TypePlot, beadsRef)

if (TypeCal == "ENT" and TypePlot == "COMBINE"):
	beadsRef = 101
	numbblocks	        = 50000
	numbmolecules       = 2
	numbpass            = 200
	particleA           = int(numbmolecules/2)
	preskip             = 1000
	postskip            = 0
	FigureGenerator.FigureENTCOMBINE(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA,TypePlot, beadsRef)


if (TypePlot == "CorrFunc"):
	beadsRef = 101
	RefPoint = [3]
	FigureGenerator.FigureCorrelation(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA, beadsRef, RefPoint)


if (TypePlot == "Energy"):
	FigureGenerator.FigureEnergyPIGS(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA, TypePlot, beadsRef)

#End plotting ---energy

if (TypePlot == "ChemPot"):
	beadsRef = 61
	FigureGenerator.FigureChemicalPotentialPIGS(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA, beadsRef)

#End plotting ---Chemical Potential
if (TypeCal == "PIGS" and TypePlot == "GFACTOR"):
	beadsRef   = 21
	numbblocks = 10000
	numbmolecules  = 2
	numbpass   = 50
	preskip    = 5000
	postskip   = 0
	particleA           = int(numbmolecules/2)
	FigureGenerator.FigureAngleDistributionGfactor(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA,TypePlot, beadsRef)
'''

'''
if (TypeCal == "PIGS" and TypePlot == "GFACTOR"):
	beadsRef   = 41
	numbblocks = 10000
	numbpass   = 50
	FigureGenerator.FigureAngleDistribution(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA,TypePlot, beadsRef)
'''
