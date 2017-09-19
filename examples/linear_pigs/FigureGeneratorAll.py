#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support
import FigureGenerator

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
#variableName        = "beta"
variableName        = "tau"
#
TransMove           = "No"
RotMove             = "Yes"
#
#TypeCal             = 'PIMC'
TypeCal             = 'PIGS'
#TypeCal             = 'ENT'
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

preskip             = 10000
postskip            = 0

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

#==================================Plotting====================================#
FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA)

if (TypeCal == "ENT"):
	ExactValueFile    = "ResultsOfPIGSENT/ExactValue-Entropy-RotDOFs-Rpt10.05Angstrom-DipoleMoment1.826Debye-Dmitri.txt"
	NumberOfMolecules, ExactEntropy = np.loadtxt(ExactValueFile, usecols=(0, 1), unpack=True)
	index = numbmolecules/2 - 1
	ExactEntropyValue = ExactEntropy[index]
	FileToBePlot   	  = FilePlotName.SaveEntropy+".txt"
	FilePlot          = FilePlotName.SaveEntropy+".pdf"
	FigureGenerator.FigureENT(FileToBePlot,FilePlot,TypeCal,variableName,parameter,ExactEntropyValue,numbmolecules,molecule,Rpt,dipolemoment)

if (TypeCal == "PIGS"):
	'''
	TypeCorr          = "Total"
	if (TypeCorr == "Total"):
		FileToBePlot   	  = "ResultsOfPIGS/totalcorr.txt"
		FilePlot          = FilePlotName.SaveTotalCorr+".pdf"
		call(["rm", FilePlot])
		FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
	TypeCorr          = "XCorr"
	if (TypeCorr == "XCorr"):
		FileToBePlot   	  = "ResultsOfPIGS/xcorr.txt"
		FilePlot          = FilePlotName.SaveXCorr+".pdf"
		call(["rm", FilePlot])
		FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
	TypeCorr          = "YCorr"
	if (TypeCorr == "YCorr"):
		FileToBePlot   	  = "ResultsOfPIGS/ycorr.txt"
		FilePlot          = FilePlotName.SaveYCorr+".pdf"
		call(["rm", FilePlot])
		FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
	TypeCorr          = "ZCorr"
	if (TypeCorr == "ZCorr"):
		FileToBePlot   	  = "ResultsOfPIGS/zcorr.txt"
		FilePlot          = FilePlotName.SaveZCorr+".pdf"
		call(["rm", FilePlot])
		FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
	TypeCorr          = "XYCorr"
	if (TypeCorr == "XYCorr"):
		FileToBePlot   	  = "ResultsOfPIGS/xycorr.txt"
		FilePlot          = FilePlotName.SaveXYCorr+".pdf"
		call(["rm", FilePlot])
		FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)

	'''
#End plotting ---correlation

	FileToBePlot   	  = FilePlotName.SaveEnergy+".txt"
	FilePlot          = FilePlotName.SaveEnergy+".pdf"
	call(["rm", FilePlot])
	FigureGenerator.FigureEnergyPIGS(FileToBePlot,FilePlot,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
exit()
#==================================End Plotting====================================#
