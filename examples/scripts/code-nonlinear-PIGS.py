import os
from subprocess import call

import support_without_parallel as support

space=" "
stringName1 = "nonlinear-rotors"
fileName1 = "script_submission_analysis_MoRiBS_without_parallel.py"
fileName2 = "script_submission_analysis_MoRiBS1.py"
support.replace("NameOfOutputDirectory", stringName1, fileName1, fileName2)

stringName2 = '""'
fileName3 = "script_submission_analysis_MoRiBS-" + stringName1 + ".py"
support.replace("extraName", stringName2, fileName2, fileName3)
call(["rm", fileName2])

# Informations about the system
simType = "PIGS"
simType1="submission -C"
#simType1 = "analysis"
molecule = "H2O"
rotor = "H2O"
var = "beta"

rcom = 10.05
nMolecule = 2
param = 0.01
nblocks = 20000
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
) 
print(cmd_run)
os.system(cmd_run)
call(["rm", fileName3])
