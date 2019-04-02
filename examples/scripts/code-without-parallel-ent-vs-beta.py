#!/usr/bin/python

import os
from os import system
import support_without_parallel as support
from subprocess import call

stringName1      = "test-ratio-trick"
fileName1        = "script_submission_analysis_MoRiBS_without_parallel.py"
fileName2        = "script_submission_analysis_MoRiBS1.py"
support.replace("NameOfOutputDirectory", stringName1, fileName1, fileName2)

#stringName2     = '"test-"'
stringName2      = '"Ratio-Trick-"'
#stringName2      = '""'
fileName3        = "script_submission_analysis_MoRiBS-"+stringName1+".py"
support.replace("extraName", stringName2, fileName2, fileName3)
call(["rm", fileName2])

# Informations about the system
molecule         = "HF"
nMolecule        = 32
RCOM             = 10.05
if (nMolecule == 2):
	gFactorList  = [1.0, 2.0, 4.0, 6.0, 8.0]
if (nMolecule == 8):
	gFactorList  = [1.0, 1.5, 2.0]
if (nMolecule == 16):
	gFactorList  = [1.0, 1.3]
if (nMolecule == 32):
	gFactorList  = [1.0, 1.1]
#tau              = 0.005
#tau              = 0.02
tau              = 0.04
Simulation       = "ENT"

for gFactor in gFactorList:
	gFactor = '{:03.1f}'.format(gFactor)

	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 20000 -Pass 100 --ROTMOVE beta submission -C "+Simulation+" HF HF "+str(tau)+" --RATIO"
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 20000 -Pass 100 --ROTMOVE beta submission -C "+Simulation+" HF HF "+str(tau)+" --RATIO --RESTART -NR 20000"
	command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 20000 -Pass 100 --ROTMOVE beta analysis   --preskip 15000 "+Simulation+" HF HF "+str(tau)+" --RATIO"
	system(command_line)
#-------------------------------#
call(["rm", fileName3])
