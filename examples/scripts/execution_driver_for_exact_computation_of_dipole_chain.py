import os
from subprocess import call
import numpy as np

import mypkg
import mypkg.moribs_runner as runner
import mypkg.moribs_runner.support as support

module_path = mypkg.moribs_runner.__file__
module_path = module_path.replace('__init__.py', '')

# User needs to modify once
dir_name = '"linear-rotors"'
 
# /home/tapas/MoRiBS-PIGS/examples/scripts
#project_dir = '""'
#plot_dir_path = '""'

# /home/tapas/academic-project/MoRiBS-PIGS/examples/scripts
project_dir = '"academic-project/"'
plot_dir_path = '"academic-project/outputs/"'

extra_file_name = '""'
blank_space = " "
#
# Informations about the system
method="dmrg"
rotor_name="HF"
numb_molecule=2
rlist = np.arange(3.0, 10.01, 0.2, dtype=float)
dipole_moment = 1.827
l_max=10
l_total_max=l_max
#
# No need to change the below three lines
original_file_name = module_path + "get_exact_result.py"
temp_file_name = "get-dmrg-temp.py"
temp_file_name1 = "get-dmrg-temp1.py"

for rpt in rlist:
	rpt_value = "{:3.2f}".format(rpt)

	support.replace(
		"name_of_output_directory",
		dir_name,
		original_file_name,
		temp_file_name)

	support.replace(
		"path_dmrg_dir",
		project_dir,
		temp_file_name,
		temp_file_name1)
	call(["mv", temp_file_name1, temp_file_name])

	support.replace(
		"extra_name",
		extra_file_name,
		temp_file_name,
		temp_file_name1)
	call(["mv", temp_file_name1, temp_file_name])

	support.replace(
		"plot_dir_path",
		plot_dir_path,
		temp_file_name,
		temp_file_name1)
	call(["mv", temp_file_name1, temp_file_name])

	cmd_run = (
		"python" + blank_space + temp_file_name + blank_space
		+ method + blank_space
		+ "--dipole_moment" + blank_space + str(dipole_moment) + blank_space
		+ "--rpt" + blank_space + rpt_value + blank_space
		+ " --rotor" + blank_space + rotor_name + blank_space
		+ " --nmolecule" + blank_space + str(numb_molecule) + blank_space
		+ "--l_max" + blank_space + str(l_max) + blank_space
		+ "--l_total_max" + blank_space + str(l_total_max) + blank_space
	)
	print(cmd_run)
	os.system(cmd_run)
	call(["rm", temp_file_name])
