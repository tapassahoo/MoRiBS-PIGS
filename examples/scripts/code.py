#!/usr/bin/python

import os
from os import system
import support

molecule     = "HF"
RCOM         = 10.05
gFactorList  = [0.5, 0.75, 1.0]
for gFactor in gFactorList:
	DipoleMoment = support.GetDipoleMomentFromGFactor(molecule, RCOM, gFactor)
	#support.GetrAndgFactor(molecule, RCOM, DipoleMoment)
	print("")
	printMessage = "Dipole Moment of "+molecule+" is "+str(DipoleMoment)+" Debye"
	print(printMessage)
	print("")
	output  = '{:1.4f}'.format(DipoleMoment)

	command_line = "python script_submission_analysis_MoRiBS.py -d "+str(output)+" -R 10.05 --scal SWAPTOUNSWAP -N 8 -Block 10000 -Pass 200 --ROTMOVE tau submission ENT HF HF 0.2"
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 --scal SWAPTOUNSWAP -N 2 -Block 200000 -Pass 200 --ROTMOVE --RESTART -NR 400000 tau submission ENT HF HF 0.2"
	command_line = "python script_submission_analysis_MoRiBS.py -d "+str(output)+" -R 10.05 --scal SWAPTOUNSWAP -N 8 -Block 10000 -Pass 200 --ROTMOVE --preskip 0 tau analysis ENT HF HF 0.2"
	system(command_line)
#-------------------------------#
	#command_line = "python script_exact_energy_entropy.py -d "+str(DipoleMoment)+" -R 10.05 -N 4 --ROTMOVE tau PIGS HF HF 0.2"
