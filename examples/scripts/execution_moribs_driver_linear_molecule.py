import os, subprocess
import numpy as np

import mypkg
import mypkg.moribs_runner as runner
import mypkg.moribs_runner.support as support

module_path = mypkg.moribs_runner.__file__
module_path = module_path.replace('__init__.py', '')

# User needs to modify once
submission_root_dir_name = "linear-rotors"
# /home/tapas/MoRiBS-PIGS/examples/scripts
dir_moribs = '""'
plot_dir_path = '""'

# /home/tapas/academic-project/MoRiBS-PIGS/examples/scripts
"""
dir_moribs = '"academic-project/"'
plot_dir_path = '"academic-project/outputs/"'
"""

extra_name = '""'
blank_space = " "

# job_type is two types - "submission" and "analysis"
#job_type = "submission"
job_type = "analysis"
method = "PIGS"

system = "HF"
rotor = "HF"
spin_isomer = int(-1)

parameter_name = "beta"
parameter_value = 0.2

numb_molecule=2
numb_block=20000
numb_pass=200
numb_preskip=0

if (numb_molecule > 1):
	dipole_moment = 1.827
elif (numb_molecule == 1):
	dipole_moment = 20.0  # It refer to the field strength and the unit inverse of Kelvin

if (job_type == "submission"):
	rlist = np.arange(3.0, 3.41, 0.2, dtype=float)

if (job_type == "analysis"):
	if (parameter_name == "beta"):
		if (parameter_value == 0.2):
			rlist = np.arange(5.0, 5.01, 0.2, dtype=float)
		if (parameter_value == 0.1):
			rlist = np.arange(10.0, 10.01, 0.2, dtype=float)

	if (parameter_name == "tau"):
		rlist = np.arange(3.0, 10.01, 1.0, dtype=float)

# No need to change the below three lines
original_file_name = module_path + "submission_analysis.py"
if (job_type=="analysis"):
	temp_file_name = "temp_0_" + submission_root_dir_name + "_" + job_type + "_n" + str(numb_molecule) + "_" + parameter_name + str(parameter_value) + "kelvin_inv_preskip" + str(numb_preskip) +".py"
	temp_file_name1 = "temp_1_" + submission_root_dir_name + "_" + job_type + "_n" + str(numb_molecule) + "_" + parameter_name + str(parameter_value) + "kelvin_inv_preskip" + str(numb_preskip) +".py"
if (job_type=="submission"):
	temp_file_name = "temp_0_" + submission_root_dir_name + "_" + job_type + "_n" + str(numb_molecule) + "_" + parameter_name + str(parameter_value) + "kelvin_inv.py"
	temp_file_name1 = "temp_1_" + submission_root_dir_name + "_" + job_type + "_n" + str(numb_molecule) + "_" + parameter_name + str(parameter_value) + "kelvin_inv.py"

for rcom in rlist:
	support.replace(
		"name_of_output_directory",
		submission_root_dir_name,
		original_file_name,
		temp_file_name)

	support.replace(
		"path_moribs_dir",
		dir_moribs,
		temp_file_name,
		temp_file_name1)
	subprocess.call(["mv", temp_file_name1, temp_file_name])

	support.replace(
		"extra_name",
		extra_name,
		temp_file_name,
		temp_file_name1)
	subprocess.call(["mv", temp_file_name1, temp_file_name])

	support.replace(
		"plot_dir_path",
		plot_dir_path,
		temp_file_name,
		temp_file_name1)
	subprocess.call(["mv", temp_file_name1, temp_file_name])

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
		+ "--restart" + blank_space
		+ "--nblock_restart" + blank_space + str(numb_block) + blank_space
		#+ "--get_energy" + blank_space
		#+ "--get_op" + blank_space
		+ "--get_itcf" + blank_space
	)

	print(cmd_run)
	os.system(cmd_run)
	subprocess.call(["rm", temp_file_name])
