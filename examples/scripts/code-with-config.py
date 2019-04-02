#!/usr/bin/python

import os
from os import system
import support_with_config as support
from subprocess import call

stringName1  = "test-config"
fileName1    = "script_submission_analysis_MoRiBS_with_config.py"
fileName2    = "script_submission_analysis_MoRiBS1.py"
support.replace("NameOfOutputDirectory", stringName1, fileName1, fileName2)

#stringName2   = '"test-"'
#stringName2   = '""'
stringName2   = '"NoRatio-Trick-"'
fileName3    = "script_submission_analysis_MoRiBS-"+stringName1+".py"
support.replace("extraName", stringName2, fileName2, fileName3)
call(["rm", fileName2])

# Informations about the system
molecule     = "HF"
nMolecule    = 8
RCOM         = 10.05
if (nMolecule == 2):
	gFactorList  = [0.5+0.1*i for i in range(76)]
if (nMolecule == 4):
	gFactorList  = [0.5+0.1*i for i in range(31)]
if (nMolecule == 8):
	gFactorList  = [0.5+0.1*i for i in range(16)]
if (nMolecule == 16):
	gFactorList  = [0.5+0.1*i for i in range(11)]
if (nMolecule == 32):
	gFactorList  = [0.5+0.1*i for i in range(7)]
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
	gFactor = '{:03.1f}'.format(gFactor)

	command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 20000 -Pass 100 --ROTMOVE tau submission -C "+Simulation+" HF HF "+str(beta)+" --RIC"
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 20000 -Pass 100 --ROTMOVE tau submission -C "+Simulation+" HF HF "+str(beta)+" --RESTART -NR 20000"
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 20000 -Pass 100 --ROTMOVE tau analysis   --preskip 15000 "+Simulation+" HF HF "+str(beta)+" --RATIO"
	system(command_line)
call(["rm", fileName3])
