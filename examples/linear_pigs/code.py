#!/usr/bin/python

import os
from os import system

DList = [1.0+0.5*i for i in range(9)]
DList += [5.5, 6.0, 6.5, 7.0, 7.5, 8.0]
print(DList)
for DipoleMoment in DList:
	command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 -N 2 -Block 50000 -Pass 100 --ROTMOVE --preskip 10000 tau analysis PIGS HF HF 0.2"
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 -N 2 -Block 50000 -Pass 100 --ROTMOVE tau submission PIGS HF HF 0.2"
	system(command_line)
	#command_line = "python script_exact_energy_entropy.py -d "+str(DipoleMoment)+" -R 10.05 -N 4 --ROTMOVE tau PIGS HF HF 0.2"
'''
#command_line = "python script_submission_analysis_MoRiBS.py -d 1.0 -N 1 -Block 100 -Pass 50 --MOVECOM -C tau submission PIGS HF HF 1.0"
#system(command_line)

command_line = "python script_PairDensityGenerator.py -d 1.0 -R 10.05 -P 20 -l-max 2 tau HF 0.2"
system(command_line)
'''
