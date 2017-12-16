#!/usr/bin/python

import os
from os import system

DList = [1.0+0.5*i for i in range(11)]
print(DList)
code_name = "submission-and-analysis-fixed-beta.py"
for DipoleMoment in DList:
	command_line = "python "+code_name+" -d "+str(DipoleMoment)+" -R 10.05 -N 4 -Block 50000 -Pass 50 --ROTMOVE --preskip 1000 tau analysis PIGS HF HF 0.2"
	#command_line = "python "+code_name+" -d "+str(DipoleMoment)+" -R 10.05 -N 2 -Block 50000 -Pass 50 --ROTMOVE tau submission PIGS HF HF 0.2"
	system(command_line)
