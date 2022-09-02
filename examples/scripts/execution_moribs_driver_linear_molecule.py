import os
from subprocess import call
import numpy as np

import mypkg
import mypkg.moribs_runner as runner
import mypkg.moribs_runner.support as support

module_path = mypkg.moribs_runner.__file__
module_path = module_path.replace('__init__.py', '')

# job_type is two types - "submission" and "analysis"
job_type    = "submission"
method      = "PIGS"

system      = "HF"
rotor       = "HF"
spin_isomer = int(-1)

parameter_name  = "beta"
parameter_value = 0.2

nmolecule   = 1
nblock      = 20000
npass       = 200

if (job_type == "submission"):
	rlist = np.arange(3.0, 10.01, 0.2, dtype=float)

if (job_type == "analysis"):
	if (parameter_name == "beta"):
		rlist = np.arange(3.0, 10.1, 1.0, dtype=float)

	if (parameter_name == "tau"):
		rlist = np.arange(3.0, 10.01, 1.0, dtype=float)

extra_name = '""'

blank_space = " "
dir_name = "linear-rotors"
original_file_name = module_path + "submission_analysis.py"
temp_file_name = "script_submission_analysis.py"
temp_file_name1 = "script_submission_analysis-" + dir_name + ".py"
for rcom in rlist:
	support.replace(
		"name_of_output_directory",
		dir_name,
		original_file_name,
		temp_file_name)

	support.replace(
		"extra_name", 
		extra_name, 
		temp_file_name, 
		temp_file_name1)

	call(["rm", temp_file_name])

	rcom="{:3.2f}".format(rcom)
	field_strength = 20.0 # Unit inverse of Kelvin

	cmd_run = (
		"python" + blank_space + temp_file_name1 + blank_space
		+ job_type + blank_space
		+ method + blank_space
		+ system + blank_space
		+ parameter_name + blank_space
		+ str(parameter_value) + blank_space
		+ "--rotor" + blank_space + rotor +blank_space
		+ "--rot_move" + blank_space
		+ "--nmolecule" + blank_space + str(nmolecule) + blank_space
		+ "--spin_isomer" + blank_space + str(spin_isomer) + blank_space
		+ "--rpt" + blank_space + str(rcom) + blank_space
		+ "--dipole_moment" + blank_space + str(field_strength) + blank_space
		+ "--nblock" + blank_space + str(nblock) + blank_space
		+ "--npass" + blank_space + str(npass) + blank_space
	)

	print(cmd_run)
	os.system(cmd_run)
	call(["rm", temp_file_name1])
