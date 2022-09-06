import os
from subprocess import call
import numpy as np

import mypkg
import mypkg.moribs_runner as runner
import mypkg.moribs_runner.support as support

module_path = mypkg.moribs_runner.__file__
module_path = module_path.replace('__init__.py', '')

# User needs to modify once
dir_name = "linear-rotors"
# /home/tapas/MoRiBS-PIGS/examples/scripts
dir_moribs = '""'
# /home/tapas/academic-project/MoRiBS-PIGS/examples/scripts
dir_moribs = '"academic-project/"'
extra_name = '""'
plot_dir_path = '"academic-project/outputs/"'
blank_space = " "

# job_type is two types - "submission" and "analysis"
#job_type = "submission"
#job_type = "analysis"
job_type = "plot"
method = "PIGS"

system = "HF"
rotor = "HF"
spin_isomer = int(-1)

parameter_name = "tau"
parameter_value = 0.001

numb_molecule = 2
numb_block = 20000
numb_pass = 200
numb_preskip = 10000

if (numb_molecule > 1):
	dipole_moment = 1.827
elif (numb_molecule == 1):
	dipole_moment = 20.0  # It refer to the field strength and the unit inverse of Kelvin

if (job_type == "plot"):
	rlist = np.arange(3.0, 10.41, 1.0, dtype=float)

if (job_type == "submission"):
	rlist = np.arange(3.0, 3.41, 0.2, dtype=float)

if (job_type == "analysis"):
	if (parameter_name == "beta"):
		rlist = np.arange(3.0, 10.1, 0.2, dtype=float)

	if (parameter_name == "tau"):
		rlist = np.arange(3.0, 10.01, 1.0, dtype=float)

# No need to change the below three lines
original_file_name = module_path + "submission_analysis.py"
temp_file_name = "script_submission_analysis-" + dir_name + "-temp.py"
temp_file_name1 = "script_submission_analysis-" + dir_name + "-temp1.py"

for rcom in rlist:
	support.replace(
		"name_of_output_directory",
		dir_name,
		original_file_name,
		temp_file_name)

	support.replace(
		"path_moribs_dir",
		dir_moribs,
		temp_file_name,
		temp_file_name1)
	call(["mv", temp_file_name1, temp_file_name])

	support.replace(
		"extra_name",
		extra_name,
		temp_file_name,
		temp_file_name1)
	call(["mv", temp_file_name1, temp_file_name])

	support.replace(
		"plot_dir_path",
		plot_dir_path,
		temp_file_name,
		temp_file_name1)
	call(["mv", temp_file_name1, temp_file_name])

	rcom = "{:3.2f}".format(rcom)

	cmd_run = (
		"python" + blank_space + temp_file_name + blank_space
		+ job_type + blank_space
		+ method + blank_space
		+ system + blank_space
		+ parameter_name + blank_space
		+ str(parameter_value) + blank_space
		+ "--rotor" + blank_space + rotor + blank_space
		+ "--rot_move" + blank_space
		+ "--nmolecule" + blank_space + str(numb_molecule) + blank_space
		+ "--spin_isomer" + blank_space + str(spin_isomer) + blank_space
		+ "--rpt" + blank_space + str(rcom) + blank_space
		+ "--dipole_moment" + blank_space + str(dipole_moment) + blank_space
		+ "--nblock" + blank_space + str(numb_block) + blank_space
		+ "--npass" + blank_space + str(numb_pass) + blank_space
		+ "--preskip" + blank_space + str(numb_preskip) + blank_space
	)

	print(cmd_run)
	os.system(cmd_run)
	call(["rm", temp_file_name])
