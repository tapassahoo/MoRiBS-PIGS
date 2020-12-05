import os
from subprocess import call
import numpy as np

import mypkg.pkgMoribs.support_without_parallel as support
import mypkg.pkgMoribs

module_path = mypkg.pkgMoribs.__file__
module_path=module_path.replace('__init__.py', '')

# Informations about the system
simType = "PIMC"
#simType1="submission"
simType1 = "analysis"
#simType1="rename "

molecule = "H2O"
rotor = "H2O"
SpinIsomer = int(-1)

#field_strength = 20.0 # Unit inverse of Kelvin
nMolecule = 1
nblocks = 50000
npass = 200

#var = "beta" # for fixed tau
#param = 0.001 # for fixed tau

var = "tau"  # for fixed beta
param = 1.0/3.0 # for fixed beta

#nbeads = 20
#param = 0.005*nbeads


if (simType1 == "submission"):
	rmin = 4.0
	rmax = 4.0
	dr = 0.1
	nr = int(((rmax-rmin)+dr*0.5)/dr)
	nr = nr+1
	print(nr)
	RList = [rmin+dr*i for i in range(nr)]
	print(RList)

if (simType1 == "analysis"):
	rmin = 4.0
	rmax = 4.0
	dr = 0.1
	nr = int(((rmax-rmin)+dr*0.5)/dr)
	nr = nr+1
	print(nr)
	RList = [rmin+dr*i for i in range(nr)]
	print(RList)


'''
if (simType1 == "analysis"):
	if ((param == 0.1) and (var == "tau")):
		rmin = 2.5
		rmax = 2.7
		dr = 0.02
		nr = int(((rmax-rmin)+dr*0.5)/dr)
		nr += 1
		RList = [rmin+dr*i for i in range(nr)]
		RList += [2.75]
		rmin = 2.8
		rmax = 5.0
		dr = 0.1
		nr = int(((rmax-rmin)+dr*0.5)/dr)
		nr += 1
		RList += [rmin+dr*i for i in range(nr)]
		rmin = 5.2
		rmax = 7.0
		dr = 0.2
		nr = int(((rmax-rmin)+dr*0.5)/dr)
		nr += 1
		print(nr)
		RList += [rmin+dr*i for i in range(nr)]
		print(RList)

	if ((param == 0.2) and (var == "tau")):
		rmin = 10.0
		rmax = 10.0
		dr = 0.2
		nr = int(((rmax-rmin)+dr*0.5)/dr)
		nr += 1
		RList = [rmin+dr*i for i in range(nr)]
'''

#stringName2 = '""'
#stringName2 = '"TIP4P-2005-"'
#stringName2 = '"qTIP4P-"'
stringName2 = '"qTIP4P-one-rotor-fixed-cost1-rzero-10e-1-"'
#stringName2 = '"qTIP4P-one-rotor-fixed-cost1-"'
#stringName2 = '"qSPCFw-"'
#stringName2 = '"qTIP4P-thread4-"'
#stringName2 = '"qTIP4P-test-"'

if simType1 == "analysis":
	cmd1 = "--preskip 0"
else:
	cmd1 = ""


for rcom in RList:

	space=" "

	stringName1 = "nonlinear-rotors"
	fileName1 = module_path+"script_submission_analysis_MoRiBS_without_parallel.py"
	fileName2 = "script_submission_analysis_MoRiBS1.py"
	support.replace("NameOfOutputDirectory", stringName1, fileName1, fileName2)

	fileName3 = "script_submission_analysis_MoRiBS-" + stringName1 + ".py"
	support.replace("extraName", stringName2, fileName2, fileName3)
	call(["rm", fileName2])

	rcom="{:3.2f}".format(rcom)
	#field_strength = 20.0 # Unit inverse of Kelvin

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
		#+ "--MOVECOM"+space
		+ cmd1+space
		+ "--Type NONLINEAR"+space
		+ "-spin"+space
		+ str(SpinIsomer)+space
		#+ "-IM"+space+"ATOM"+space+"H2"+space+"1"+space
		+ simType+space
		+ simType1+space
		+ molecule+space
		+ rotor+space
		+ str(param)+space
		+ var+space
		#+ " -C --RESTART"+space
		#+ " -NR 20000"+space
	) 

	print(cmd_run)
	os.system(cmd_run)
	call(["rm", fileName3])
