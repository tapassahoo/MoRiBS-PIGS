#!/usr/bin/python

import os
from os import system

#DList = [1.0+0.5*i for i in range(2)]
DList = [1.8548]
print(DList)
for DipoleMoment in DList:
##For PIMC:
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 6.0 -N 2 -Block 100000 -Pass 100 --ROTMOVE tau submission PIMC H2O H2O 0.0333333333"
	command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 6.0 -N 2 -Block 100000 -Pass 100 --ROTMOVE --preskip 10000 tau analysis PIMC H2O H2O 0.0333333333"
#-------------------------------#
##For PIGS:
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 -N 2 -Block 2000 -Pass 20 --ROTMOVE tau submission PIGS HF HF 0.2"
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 -N 2 -Block 50000 -Pass 400 --ROTMOVE --preskip 10000 tau analysis PIGS HF HF 0.2"
#-------------------------------#
	system(command_line)
