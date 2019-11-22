import os
from subprocess import call
import numpy as np

import mypkg.pkgMoribs.support_without_parallel as support
import mypkg.pkgMoribs

module_path = mypkg.pkgMoribs.__file__
module_path=module_path.replace('__init__.py', '')

for i in range(16):

	space=" "

	stringName1 = "nonlinear-rotors"
	fileName1 = module_path+"script_submission_analysis_MoRiBS_without_parallel.py"
	fileName2 = "script_submission_analysis_MoRiBS1.py"
	support.replace("NameOfOutputDirectory", stringName1, fileName1, fileName2)

	#stringName2 = '""'
	stringName2 = '"TIP4P-2005-"'
	fileName3 = "script_submission_analysis_MoRiBS-" + stringName1 + ".py"
	support.replace("extraName", stringName2, fileName2, fileName3)
	call(["rm", fileName2])

	# Informations about the system
	simType = "PIGS"

	#simType1="submission "
	simType1 = "analysis"

	molecule = "H2O"
	rotor = "H2O"
	SpinIsomer = 0

	var = "beta" # for fixed tau
	param = 0.001 # for fixed tau

	#var = "tau"  # for fixed beta
	#param = 0.2 # for fixed beta

	rcom = 2.5+0.5*i
	rcom="{:3.1f}".format(rcom)
	#field_strength = 20.0 # Unit inverse of Kelvin
	nMolecule = 11
	nblocks = 10000
	npass = 100

	if simType1 == "analysis":
		cmd1 = "--preskip 0"
	else:
		cmd1 = ""

	cmd_run = (
		"python"+space
		+ fileName3+space
		+ "-R"+space
		+ str(rcom)+space
		#+ "-d"+space
		#+ str(field_strength)+space
		+ "-N"+space
		+ str(nMolecule)+space
		+ "-Block"+space
		+ str(nblocks)+space
		+ "-Pass"+space
		+ str(npass)+space
		+ "--ROTMOVE"+space
		+ cmd1+space
		+ "--Type NONLINEAR"+space
		+ var+space
		+ simType1+space
		+ simType+space
		+ molecule+space
		+ rotor+space
		+ str(param)+space
		+ "-spin"+space
		+ str(SpinIsomer)+space
	) 
	print(cmd_run)
	os.system(cmd_run)
	call(["rm", fileName3])
