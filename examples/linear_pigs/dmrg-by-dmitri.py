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
numbblocks	        = 50000
numbmolecules       = 2
numbpass            = 200
#
Rpt                 = 10.05
dipolemoment        = 2.0        #J. Chern. Phys. 73(5), 2319 (1980).
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

if (variableName == "beta"):
	parameterName   = "tau"
	tau             = 0.005
	parameter       = tau

preskip             = 0
postskip            = 0
#==================================Plotting====================================#
FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA)

if (TypeCal == "ENT"):
	if (variableName == "tau"):
		FileToBeSaved     = FilePlotName.SaveEntropyDMRG+".txt"

		'''
		#julia run
		#DipoleMomentList = [1.0+i*0.25 for i in range(13)]
		for dipolemoment in DipoleMomentList:
			RFactor      = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
			srcCodePath  = "/home/tapas/DipoleChain.jl-master/examples/"
			commandRun   = "julia "+srcCodePath+"diagonalization.jl -R "+str(RFactor)+" -N "+str(numbmolecules)+" --l-max 4 --A-start "+str(particleA)+" --A-size "+str(particleA)
			system(commandRun)
		'''
		#julia run
		call(["rm", "output.txt"])
		loop = [i*10+1 for i in range(1,11)]
		for numbbeads in loop:
			print(numbbeads)
			RFactor      = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
			srcCodePath  = "/home/tapas/DipoleChain.jl-master/examples/"
			Units        = support.GetUnitConverter()
			BConstant    = support.GetBconst(molecule)  # in wavenumber
			BConstantK   = BConstant*Units.CMRECIP2KL
			betaRedcd    = beta*BConstantK
			commandRun   = "julia "+srcCodePath+"path_integral.jl -R "+str(RFactor)+" -N "+str(numbmolecules)+" --l-max 8 --beta "+str(betaRedcd)+" -P "+str(numbbeads)+" --pigs --A-start "+str(particleA)+" --A-size "+str(particleA)
			system(commandRun)
		call(["mv", "output.txt", FileToBeSaved])
		exit()
