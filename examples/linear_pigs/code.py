#!/usr/bin/python

import os
from os import system

#DList = [1.0+0.5*i for i in range(2)]
DList = [6.0]
print(DList)
for DipoleMoment in DList:
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 -N 2 -Block 500 -Pass 200 --ROTMOVE tau analysis PIMC HF HF 0.2"
#-------------------------------#
##For PIGS:
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 -N 2 -Block 2000 -Pass 20 --ROTMOVE tau submission PIGS HF HF 0.2"
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 -N 2 -Block 2000 -Pass 20 --ROTMOVE --RESTART -NR 4000 tau submission PIGS HF HF 0.2"
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 -N 2 -Block 50000 -Pass 400 --ROTMOVE --preskip 10000 tau analysis PIGS HF HF 0.2"
#-------------------------------#
##For ENT
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 --scal SWAPTOUNSWAP -N 2 -Block 200000 -Pass 200 --ROTMOVE tau submission ENT HF HF 0.2"
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 --scal SWAPTOUNSWAP -N 2 -Block 200000 -Pass 200 --ROTMOVE --RESTART -NR 400000 tau submission ENT HF HF 0.2"
	command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 --scal SWAPTOUNSWAP -N 2 -Block 20000 -Pass 200 --ROTMOVE --preskip 0 tau analysis ENT HF HF 0.2"
#-------------------------------#
	system(command_line)
	#command_line = "python script_exact_energy_entropy.py -d "+str(DipoleMoment)+" -R 10.05 -N 4 --ROTMOVE tau PIGS HF HF 0.2"
'''
command_line = "python script_submission_analysis_MoRiBS.py -d 1.0 -R 10.05 -N 2 -Block 100 -Pass 50 --ROTMOVE -lmax 2 -C tau submission PIGS HF PPA HF 0.2"
system(command_line)
'''

'''
command_line = "python script_PairDensityGenerator.py -d 1.0 -R 10.05 -P 20 -l-max 2 tau HF 0.2"
system(command_line)
'''
