#!/usr/bin/python

import os
from os import system
import support_without_parallel as support
from subprocess import call

#stringName1  = "test-no-ratio-trick"
stringName1  = "test-ratio-trick"
fileName1    = "script_submission_analysis_MoRiBS_without_parallel.py"
fileName2    = "script_submission_analysis_MoRiBS1.py"
support.replace("NameOfOutputDirectory", stringName1, fileName1, fileName2)

#stringName2   = '"test-"'
stringName2   = '""'
#stringName2   = '"NoRatio-Trick-"'
#stringName2   = '"Ratio-Trick-"'
fileName3    = "script_submission_analysis_MoRiBS-"+stringName1+".py"
support.replace("extraName", stringName2, fileName2, fileName3)
call(["rm", fileName2])

# Informations about the system
molecule     = "HF"
nMolecule    = 8
RCOM         = 10.05
if (nMolecule == 2):
	#gFactorList  = [0.5+0.1*i for i in range(76)]
	gFactorList  = [1.0+1.0*i for i in range(20)]
	#gFactorList  = [9.0+1.0*i for i in range(12)]
	#gFactorList  = [21.0+1.0*i for i in range(10)]
if (nMolecule == 4):
	#gFactorList  = [0.5+0.1*i for i in range(31)]
	gFactorList = [3.5+0.25*i for i in range(11)]
if (nMolecule == 8):
	gFactorList  = [0.5+0.1*i for i in range(16)]
	#gFactorList += [2.1+0.1*i for i in range(10)]
if (nMolecule == 16):
	gFactorList  = [0.5+0.1*i for i in range(11)]
if (nMolecule == 32):
	gFactorList  = [0.75+0.05*i for i in range(10)]
	#gFactorList  = [0.5+0.05*i for i in range(15)]
	#gFactorList  = [0.7]
#print(gFactorList)
#exit()
#beta         = 0.32
beta         = 0.2
Simulation   = "ENT"

for gFactor in gFactorList:
	#DipoleMoment = support.GetDipoleMomentFromGFactor(molecule, RCOM, gFactor)
	#support.GetrAndgFactor(molecule, RCOM, DipoleMoment)
	#print("")
	#printMessage = "Dipole Moment of "+molecule+" is "+str(DipoleMoment)+" Debye"
	#print(printMessage)
	#print("")
	#output  = '{:1.4f}'.format(DipoleMoment)
	gFactor = '{:03.2f}'.format(gFactor)

######For SWAPUNSWAP with Ratio Trick######
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 21000 -Pass 100 --ROTMOVE tau submission -C "+Simulation+" HF HF "+str(beta)+" --RATIO"
	command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 20000 -Pass 100 --ROTMOVE tau submission -C "+Simulation+" HF HF "+str(beta)+" --RESTART -NR 20000 --RATIO"
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 21000 -Pass 100 --ROTMOVE tau rename "+Simulation+" HF HF "+str(beta)+" -NR 21000 --RATIO"
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 20000 -Pass 100 --ROTMOVE tau analysis   --preskip 0 "+Simulation+" HF HF "+str(beta)+" --RATIO"
######For BROKENPATH######
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 50000 -Pass 200 --ROTMOVE tau submission -C "+Simulation+" HF HF "+str(beta)+" --scal BROKENPATH"
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 50000 -Pass 200 --ROTMOVE tau submission -C "+Simulation+" HF HF "+str(beta)+" --RESTART -NR 50000 --scal BROKENPATH"
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 50000 -Pass 100 --ROTMOVE tau rename         "+Simulation+" HF HF "+str(beta)+" -NR 100000 --scal BROKENPATH"
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 50000 -Pass 200 --ROTMOVE tau analysis   --preskip 40000 "+Simulation+" HF HF "+str(beta)+" --scal BROKENPATH"
######For SWAPUNSWAP without Ratio Trick######
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 50000 -Pass 200 --ROTMOVE tau submission -C "+Simulation+" HF HF "+str(beta)
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 50000 -Pass 200 --ROTMOVE tau submission -C "+Simulation+" HF HF "+str(beta)+" --RESTART -NR 50000"
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 50000 -Pass 200 --ROTMOVE tau analysis --preskip 0 "+Simulation+" HF HF "+str(beta)
	system(command_line)
call(["rm", fileName3])
