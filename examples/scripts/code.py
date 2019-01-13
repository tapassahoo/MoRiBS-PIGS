#!/usr/bin/python

import os
from os import system
import support
from subprocess import call

stringName1   = "linear-rotors"
fileName1    = "script_submission_analysis_MoRiBS.py"
fileName2    = "script_submission_analysis_MoRiBS1.py"
support.replace("NameOfOutputDirectory", stringName1, fileName1, fileName2)

stringName2   = '""'
fileName3    = "script_submission_analysis_MoRiBS-"+stringName1+".py"
support.replace("extraName", stringName2, fileName2, fileName3)
call(["rm", fileName2])

#
molecule     = "HF"
nMolecule    = 16
RCOM         = 10.05
gFactorList  = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]
beta         = 0.2
Erot         = support.GetAvgRotEnergy(molecule,beta)
Simulation   = "PIGS"
print("")
print("Avg. Rotational Energy of "+str(nMolecule)+" free "+molecule+" is "+str(nMolecule*Erot)+" Kelvin")
print("")
for gFactor in gFactorList:
	DipoleMoment = support.GetDipoleMomentFromGFactor(molecule, RCOM, gFactor)
	#support.GetrAndgFactor(molecule, RCOM, DipoleMoment)
	print("")
	printMessage = "Dipole Moment of "+molecule+" is "+str(DipoleMoment)+" Debye"
	print(printMessage)
	print("")
	output  = '{:1.4f}'.format(DipoleMoment)

	#command_line = "python "+fileName3+" -d "+str(output)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 10000 -Pass 50 --ROTMOVE tau submission "+Simulation+" HF HF "+str(beta)
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 --scal SWAPTOUNSWAP -N 2 -Block 200000 -Pass 200 --ROTMOVE --RESTART -NR 400000 tau submission ENT HF HF 0.2"
	command_line = "python "+fileName3+" -d "+str(output)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 10000 -Pass 50 --ROTMOVE tau analysis   "+Simulation+" HF HF "+str(beta)
	system(command_line)
#-------------------------------#
	#command_line = "python script_exact_energy_entropy.py -d "+str(DipoleMoment)+" -R 10.05 -N 4 --ROTMOVE tau PIGS HF HF 0.2"
call(["rm", fileName3])
