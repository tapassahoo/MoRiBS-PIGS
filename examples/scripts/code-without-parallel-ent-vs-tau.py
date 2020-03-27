import os
from subprocess import call
import numpy as np

import mypkg.pkgMoribs.support_without_parallel as support
import mypkg.pkgMoribs

module_path = mypkg.pkgMoribs.__file__
module_path=module_path.replace('__init__.py', '')

# Informations about the system
simType = "ENT"
simType1="submission "
#simType1 = "analysis"

molecule = "HF"
rotor = "HF"
SpinIsomer = -1

#var = "beta" # for fixed tau
#param = 0.002 # for fixed tau

var = "tau"  # for fixed beta
param = 0.2 # for fixed beta

nMolecule = 2
nblocks = 80
npass = 100

stringName2 = '""'
#stringName2 = '"TIP4P-2005-"'

if simType1 == "analysis":
	cmd1 = "--preskip 0"
else:
	cmd1 = ""

stringName1 = "linear-rotors"
rcom = 10.05
rcom="{:3.2f}".format(rcom)

if (nMolecule == 2):
	#gFactorList  = [0.5+0.5*i for i in range(16)]
	gFactorList  = [9.0+1.0*i for i in range(12)]
if (nMolecule == 4):
	#gFactorList  = [0.5+0.1*i for i in range(31)]
	gFactorList = [3.5+0.25*i for i in range(11)]
if (nMolecule == 8):
	gFactorList  = [0.5+0.1*i for i in range(16)]
	#gFactorList += [2.1+0.1*i for i in range(10)]
if (nMolecule == 16):
	gFactorList  = [0.5+0.1*i for i in range(11)]
if (nMolecule == 32):
	gFactorList  = [0.75+0.05*i for i in range(10)]
	#gFactorList  = [0.5+0.05*i for i in range(15)]
	#gFactorList  = [0.7]
print(gFactorList)

for gFactor in gFactorList:
	#DipoleMoment = support.GetDipoleMomentFromGFactor(molecule, RCOM, gFactor)
	#support.GetrAndgFactor(molecule, RCOM, DipoleMoment)
	#print("")
	#printMessage = "Dipole Moment of "+molecule+" is "+str(DipoleMoment)+" Debye"
	#print(printMessage)
	#print("")
	#output  = '{:1.4f}'.format(DipoleMoment)
	gFactor = '{:03.2f}'.format(gFactor)

	space=" "

	fileName1 = module_path+"script_submission_analysis_MoRiBS_without_parallel.py"
	fileName2 = "script_submission_analysis_MoRiBS1.py"
	support.replace("NameOfOutputDirectory", stringName1, fileName1, fileName2)

	fileName3 = "script_submission_analysis_MoRiBS-" + stringName1 + ".py"
	support.replace("extraName", stringName2, fileName2, fileName3)
	call(["rm", fileName2])

	cmd_run = (
		"python"+space
		+ fileName3+space
		+ "-R"+space
		+ str(rcom)+space
		+ "-g"+space
		+ str(gFactor)+space
		+ "-N"+space
		+ str(nMolecule)+space
		+ "-Block"+space
		+ str(nblocks)+space
		+ "-Pass"+space
		+ str(npass)+space
		+ "--ROTMOVE"+space
		#+ "--MOVECOM"+space
		+ cmd1+space
		+ "--Type LINEAR"+space
		+ "-spin"+space
		+ str(SpinIsomer)+space
		#+ "-IM"+space+"ATOM"+space+"H2"+space+"1"+space
		+ simType+space
		+ simType1+space
		+ molecule+space
		+ rotor+space
		+ str(param)+space
		+ var+space
		+ " --RATIO WOR"+space
		#+ " --scal BROKENPATH"+space
	) 

	print(cmd_run)
	os.system(cmd_run)
	call(["rm", fileName3])
