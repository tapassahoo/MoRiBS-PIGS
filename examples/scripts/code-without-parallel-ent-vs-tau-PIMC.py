import os
from os import system
import support_without_parallel as support
from subprocess import call

stringName1   = "linear-rotors-PIMC"
fileName1    = "script_submission_analysis_MoRiBS_without_parallel.py"
fileName2    = "script_submission_analysis_MoRiBS1.py"
support.replace("NameOfOutputDirectory", stringName1, fileName1, fileName2)

stringName2   = '""'
fileName3    = "script_submission_analysis_MoRiBS-"+stringName1+".py"
support.replace("extraName", stringName2, fileName2, fileName3)
call(["rm", fileName2])

#
molecule     = "HF"
nMolecule    = 3
RCOM         = 10.05
gFactorList  = [1.0+1.0*i for i in range(4)]
print(gFactorList)
beta         = 0.05
Erot         = support.GetAvgRotEnergy(molecule,beta)
Simulation   = "PIMC"
print("")
print("Avg. Rotational Energy of "+str(nMolecule)+" free "+molecule+" is "+str(nMolecule*Erot)+" Kelvin")
print("")
for gFactor in gFactorList:
	DipoleMoment = gFactor #support.GetDipoleMomentFromGFactor(molecule, RCOM, gFactor)
	#support.GetrAndgFactor(molecule, RCOM, DipoleMoment)
	print("")
	printMessage = "Dipole Moment of "+molecule+" is "+str(DipoleMoment)+" Debye"
	print(printMessage)
	print("")
	output  = float(support.dropzeros(DipoleMoment))
	print(output)

	#command_line = "python "+fileName3+" -d "+str(output)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 50000 -Pass 200 --ROTMOVE tau submission "+Simulation+" HF HF "+str(beta)
	#command_line = "python script_submission_analysis_MoRiBS.py -d "+str(DipoleMoment)+" -R 10.05 --scal SWAPTOUNSWAP -N 2 -Block 200000 -Pass 200 --ROTMOVE --RESTART -NR 400000 tau submission ENT HF HF 0.2"
	command_line = "python "+fileName3+" -d "+str(output)+" -R "+str(RCOM)+" -N "+str(nMolecule)+" -Block 50000 -Pass 200 --ROTMOVE tau analysis   --preskip 0 "+Simulation+" HF HF "+str(beta)
	system(command_line)
#-------------------------------#
	#command_line = "python script_exact_energy_entropy.py -d "+str(DipoleMoment)+" -R 10.05 -N 4 --ROTMOVE tau PIGS HF HF 0.2"
call(["rm", fileName3])
