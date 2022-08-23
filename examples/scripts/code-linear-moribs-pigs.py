import os
from subprocess import call
import numpy as np

import mypkg
import mypkg.moribs_runner as runner
import mypkg.moribs_runner.support as support

module_path = mypkg.moribs_runner.__file__
module_path = module_path.replace('__init__.py', '')

job_type    = "submission"
method      = "PIGS"

system      = "HF"
rotor       = "HF"
spin_isomer = 0

parameter_name  = "beta"
parameter_value = 0.2

nmolecule   = 1
nblock      = 20000
npass       = 200
# field_strength = 20.0 # Unit inverse of Kelvin


if (job_type == "submission"):
	rlist = np.arange(2.5, 10.1, 10.1, dtype=float)

if (job_type == "analysis"):
	if (variable_name == "beta"):
		rlist = np.arange(3.0, 10.1, 11.0, dtype=float)

	if (variable_name == "tau"):
		if (parameter_value == 0.1):
			rlist = np.arange(2.5, 10.1, 0.1, dtype=float)
		if (parameter_value == 0.2):
			rlist = np.arange(2.5, 10.1, 0.1, dtype=float)

extra_name = '""'

space = " "
for rcom in rlist:
	dir_name = "linear-rotors"
	original_file_name = module_path + "submission_analysis.py"
	temp_file_name = "script_submission_analysis.py"
	support.replace(
		"name_of_output_directory",
		dir_name,
		original_file_name,
		temp_file_name)

	temp_file_name1 = "script_submission_analysis-" + dir_name + ".py"
	support.replace("extra_name", extra_name, temp_file_name, temp_file_name1)
	call(["rm", temp_file_name])

	rcom="{:3.2f}".format(rcom)
	field_strength = 20.0 # Unit inverse of Kelvin

	cmd_run = (
		"python" + space + temp_file_name1 + space
		+ job_type + space
		+ method + space
		+ system + space
		+ parameter_name + space
		+ str(parameter_value) + space
		+ "--rotor" + space + rotor +space
		+ "--rot_move" + space
		+ "--nmolecule" + space + str(nmolecule) + space
		+ "--spin_isomer" + space + str(spin_isomer) + space
		+ "--rpt" + space + str(rcom) + space
		+ "--dipole_moment" + space + str(field_strength) + space
		+ "--nblock" + space + str(nblock) + space
		+ "--npass" + space + str(npass) + space
	)

	print(cmd_run)
	os.system(cmd_run)
	call(["rm", temp_file_name1])
