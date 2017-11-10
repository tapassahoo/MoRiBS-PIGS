#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support
import math

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
variableName        = "tau"
#variableName        = "beta"
#
TransMove           = "No"
RotMove             = "Yes"
#
#TypeCal             = 'PIMC'
#TypeCal             = 'PIGS'
TypeCal             = 'ENT'
#
TypePlot            = "Energy"
#TypePlot            = "ChemPot"
#TypePlot            = "CorrFunc"
#
#molecule            = "HFC60"                                                  
molecule            = "HF"                                                      
#molecule            = "H2"                                                    
molecule_rot        = "HF"
#
numbblocks	        = 100000
numbmolecules       = 6
numbpass            = 100
#
Rpt                 = 10.05
dipolemoment        = 1.826        #J. Chern. Phys. 73(5), 2319 (1980).
dipolemoment        = 1.0*dipolemoment


ENT_TYPE 			= "SWAPTOUNSWAP"
#ENT_TYPE 			= "SWAP"
#ENT_TYPE 			= "BROKENPATH"
#ENT_TYPE 			= "REGULARPATH"
particleA           = int(numbmolecules/2)

extra_file_name     = ""

src_dir             = os.getcwd()
if (variableName == "tau"):
	parameterName   = "beta"
	beta            = 0.2
	parameter       = beta
	temperature     = 1.0/beta   
	#loop = [i*10+1 for i in range(1,12)]
	loop = [41,51,61,71,81,91,101]

if (variableName == "beta"):
	parameterName   = "tau"
	tau             = 0.005
	parameter       = tau
	loop = [3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51]

preskip             = 0
postskip            = 0
#==================================Plotting====================================#
srcCodePath         = "/home/tapas/DipoleChain.jl-master/examples/"
Units               = support.GetUnitConverter()
BConstant           = support.GetBconst(molecule)  # in wavenumber
BConstantK          = BConstant*Units.CMRECIP2KL
################################################################################
beadsRef = 61 # here we don't need it but it is needed to pass it into FilePlotName to keep the same arguments
DList = [1.0+0.5*i for i in range(7)]
#DList += [4.5, 5.0, 6.0, 7.0]

for dipolemoment in DList:
	RFactor             = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
	FilePlotName        = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef)
	support.GetExactValues(FilePlotName, srcCodePath, RFactor, numbmolecules, loop, particleA, molecule_rot, Rpt, dipolemoment, parameter, BConstantK, variableName, TypeCal)
