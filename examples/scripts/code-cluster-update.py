#!/usr/bin/python

import os
from os import system
import support
from subprocess import call

stringName1   = "test-cluster-update"
fileName1    = "script_submission_analysis_MoRiBS.py"
fileName2    = "script_submission_analysis_MoRiBS1.py"
support.replace("NameOfOutputDirectory", stringName1, fileName1, fileName2)

stringName2   = '""'
#stringName2   = '"cluster-update-"'
fileName3    = "script_submission_analysis_MoRiBS-"+stringName1+".py"
support.replace("extraName", stringName2, fileName2, fileName3)
call(["rm", fileName2])

# Informations about the system
molecule     = "HF"
nMolecule    = 2
RCOM         = 10.05
'''
if (nMolecule == 4):
	gFactorList  = [0.5+0.25*i for i in range(15)]   # 4 HF 
if (nMolecule == 2):
	gFactorList  = [0.5+0.25*i for i in range(21)]  # 2 HF
	gFactorList += [5.75+0.25*i for i in range(14)] # 2 HF
if (nMolecule >= 16):
	gFactorList  = [0.4+0.1*i for i in range(7)]  # 2 HF
gFactorList  = [5.0]
beta         = 0.05
Simulation   = "PIMC"
for gFactor in gFactorList:
	#DipoleMoment = support.GetDipoleMomentFromGFactor(molecule, RCOM, gFactor)
	#support.GetrAndgFactor(molecule, RCOM, DipoleMoment)
	#print("")
	#printMessage = "Dipole Moment of "+molecule+" is "+str(DipoleMoment)+" Debye"
	#print(printMessage)
	#print("")
	#output  = '{:1.4f}'.format(DipoleMoment)

	command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 5 -Pass 100 --ROTMOVE tau submission -C "+Simulation+" HF HF "+str(beta)
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 --scal SWAPTOUNSWAP -N 2 -Block 200000 -Pass 200 --ROTMOVE --RESTART -NR 400000 tau submission ENT HF HF 0.2"
	#command_line = "python "+fileName3+" -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 5000 -Pass 100 --ROTMOVE tau analysis   --preskip 0 "+Simulation+" HF HF "+str(beta)
	system(command_line)
	#command_line = "python script_exact_energy_entropy.py -d "+str(DipoleMoment)+" -R 10.05 -N 4 --ROTMOVE tau PIGS HF HF 0.2"
call(["rm", fileName3])
#-------------------------------#
'''
beta         = 0.05
Simulation   = "PIMC"
for i in range(1):
	dm = 2.0+i*2.0
	dm = 1.86
	dm = '{:1.2f}'.format(dm)

	#command_line = "python "+fileName3+" -d "+str(dm)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 50000 -Pass 100 --ROTMOVE tau submission -C "+Simulation+" HF HF "+str(beta)
	command_line = "python "+fileName3+" -d "+str(dm)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 50000 -Pass 100 --ROTMOVE tau analysis "+Simulation+" HF HF "+str(beta)
	system(command_line)
call(["rm", fileName3])
