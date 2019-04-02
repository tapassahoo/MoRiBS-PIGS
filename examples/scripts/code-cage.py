#!/usr/bin/python

import os
from os import system
import support
from subprocess import call

stringName1   = "cage-pot"
fileName1    = "script_submission_analysis_MoRiBS.py"
fileName2    = "script_submission_analysis_MoRiBS1.py"
support.replace("NameOfOutputDirectory", stringName1, fileName1, fileName2)

stringName2   = '""'
fileName3    = "script_submission_analysis_MoRiBS-"+stringName1+".py"
support.replace("extraName", stringName2, fileName2, fileName3)
call(["rm", fileName2])

#
molecule     = "HF"
nMolecule    = 1
beta         = 0.2
Simulation   = "PIGS"

command_line = "python "+fileName3+" -N "+str(nMolecule)+" -Block 1000 -Pass 50 --ROTMOVE tau submission -C "+Simulation+" HF@C60 HF "+str(beta)
#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 --scal SWAPTOUNSWAP -N 2 -Block 200000 -Pass 200 --ROTMOVE --RESTART -NR 400000 tau submission ENT HF HF 0.2"
#command_line = "python "+fileName3+" -N "+str(nMolecule)+" -Block 1000 -Pass 50 --MOVECOM --ROTMOVE tau analysis   --preskip 0 "+Simulation+" HF@C60 HF "+str(beta)
system(command_line)
#-------------------------------#
	#command_line = "python script_exact_energy_entropy.py -d "+str(DipoleMoment)+" -R 10.05 -N 4 --ROTMOVE tau PIGS HF HF 0.2"
call(["rm", fileName3])
