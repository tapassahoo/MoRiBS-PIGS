#!/usr/bin/python

import os
from os import system
import support
from subprocess import call

stringName1   = "harmonic-potential"
fileName1    = "script_submission_analysis_MoRiBS.py"
fileName2    = "script_submission_analysis_MoRiBS1.py"
support.replace("NameOfOutputDirectory", stringName1, fileName1, fileName2)

stringName2   = "harmonic-potential"
fileName3    = "script_submission_analysis_MoRiBS-"+stringName1+".py"
support.replace("extraName", stringName2, fileName2, fileName3)
call(["rm", fileName2])

molecule     = "HF"
nMolecule    = 1
beta         = 0.2
Simulation   = "PIGS"

#-------------------------------#
command_line = "python "+fileName3+" -N "+str(nMolecule)+" -Block 10000 -Pass 50 --MOVECOM -C tau submission "+Simulation+" HF HF "+str(beta) 
#command_line = "python "+fileName3+" -N "+str(nMolecule)+" -Block 10000 -Pass 50 --MOVECOM tau analysis   "+Simulation+" HF HF "+str(beta)
system(command_line)
#-------------------------------#
call(["rm", fileName3])
