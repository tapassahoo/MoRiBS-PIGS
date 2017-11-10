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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
#TypePlot            = "Energy"
#TypePlot            = "ChemPot"
#TypePlot            = "CorrFunc"
#TypePlot            = "S2"
TypePlot            = "GFACTOR"
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
dipolemoment        = 1.826      #J. Chern. Phys. 73(5), 2319 (1980).
dipolemoment        = 1.0*dipolemoment

preskip             = 1000
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
if (TypeCal == "ENT"):
	beadsRef = 61
	FigureGenerator.FigureENT(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA,TypePlot, beadsRef)

if (TypeCal == "PIGS"):
	beadsRef = 61
	FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef)
	if (TypePlot == "CorrFunc"):
		RefPoint          = 0
		TypeCorr          = "Total"
		if (TypeCorr == "Total"):
			FileToBePlot   	  = FilePlotName.SaveTotalCorr+".txt"
			FilePlot          = FilePlotName.SaveTotalCorr+".pdf"
			call(["rm", FilePlot])
			FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment,RefPoint)
		TypeCorr          = "XCorr"
		if (TypeCorr == "XCorr"):
			FileToBePlot   	  = FilePlotName.SaveXCorr+".txt"
			FilePlot          = FilePlotName.SaveXCorr+".pdf"
			call(["rm", FilePlot])
			FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment,RefPoint)
		TypeCorr          = "YCorr"
		if (TypeCorr == "YCorr"):
			FileToBePlot   	  = FilePlotName.SaveYCorr+".txt"
			FilePlot          = FilePlotName.SaveYCorr+".pdf"
			call(["rm", FilePlot])
			FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment,RefPoint)
		TypeCorr          = "ZCorr"
		if (TypeCorr == "ZCorr"):
			FileToBePlot   	  = FilePlotName.SaveZCorr+".txt"
			FilePlot          = FilePlotName.SaveZCorr+".pdf"
			call(["rm", FilePlot])
			FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment,RefPoint)
		TypeCorr          = "XYCorr"
		if (TypeCorr == "XYCorr"):
			FileToBePlot   	  = FilePlotName.SaveXYCorr+".txt"
			FilePlot          = FilePlotName.SaveXYCorr+".pdf"
			call(["rm", FilePlot])
			FigureGenerator.FigureCorrelationPIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment,RefPoint)

#End plotting ---correlation

	if (TypePlot == "Energy"):
		FigureGenerator.FigureEnergyPIGS(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, TypePlot, beadsRef)

#End plotting ---energy

	if (TypePlot == "ChemPot"):
		beadsRef = 61
		FigureGenerator.FigureChemicalPotentialPIGS(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef)

#End plotting ---Chemical Potential
