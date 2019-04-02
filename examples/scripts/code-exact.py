#!/usr/bin/python

import os
from os import system
import support
from subprocess import call

# Informations about the system
nMolecule    = 4
RCOM         = 10.05
lmax         = 10
ltotalmax    = 10
Rotor        = " HF "
Simulation   = " ENT "
#
if (nMolecule == 2):
	gFactorList  = [0.5+0.1*i for i in range(76)]
if (nMolecule == 4):
	gFactorList  = [0.5+0.1*i for i in range(31)]
#
for gFactor in gFactorList:
	gFactor = '{:03.1f}'.format(gFactor)
	command_line = "python script_exact_energy_entropy.py -g "+str(gFactor)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -lmax "+str(lmax)+" -ltotalmax "+str(ltotalmax)+" --Rotor"+Rotor+Simulation
	system(command_line)
