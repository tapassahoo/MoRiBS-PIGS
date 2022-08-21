import os
from subprocess import call
import numpy as np

import mypkg
import mypkg.moribs_runner as runner
import mypkg.moribs_runner.support as support

module_path = mypkg.__file__
module_path = module_path.replace('__init__.py', '')

job_type = "submission"
method = "PIGS"

system = "HF"
rotor = "HF"
spin_isomer = int(-1)

# field_strength = 20.0 # Unit inverse of Kelvin
nmolecule = 1
nblock = 20000
npass = 200

parameter_value = 0.2
variable_name = "tau"

if (job_type == "submission"):
    rlist = np.arange(2.5, 10.1, 0.1, dtype=float)

if (job_type == "analysis"):
    if (variable_name == "beta"):
        rlist = np.arange(3.0, 10.1, 1.0, dtype=float)
    if (variable_name == "tau"):
        if (parameter_value == 0.1):
            rlist = np.arange(2.5, 10.1, 0.1, dtype=float)
        if (parameter_value == 0.2):
            rlist = np.arange(2.5, 10.1, 0.1, dtype=float)

string_layer2 = '""'

space = " "
for rcom in rlist:
    dir_layer1 = "linear-rotors"
    path_file_name = module_path + "submission_analysis.py"
    temp_file_name = "script_submission_analysis.py"
    print(path_file_name)
    print(temp_file_name)
    support.replace(
        "name_of_output_directory",
        dir_layer1,
        path_file_name,
        temp_file_name)

    #file_name3 = "script_submission_analysis-" + string_layer1 + ".py"
    #support.replace("extraName", stringName2, fileName2, fileName3)
    #call(["rm", fileName2])
#
#	rcom="{:3.2f}".format(rcom)
#	field_strength = 20.0 # Unit inverse of Kelvin
#	exit()
#
#	cmd_run = (
#		"python" + space
#		+ fileName3 + space
#		+ "-R" + space
#		+ str(rcom) + space
#		+ "-d" + space
#		+ str(field_strength) + space
#		+ "-N" + space
#		+ str(nMolecule) + space
#		+ "-Block" + space
#		+ str(nblocks) + space
#		+ "-Pass" + space
#		+ str(npass) + space
#		+ "--ROTMOVE" + space
#		#+ "--MOVECOM"+space
#		+ nskip + space
#		+ "--Type LINEAR" + space
#		+ "-spin" + space
#		+ str(SpinIsomer) + space
#		#+ "-IM"+space+"ATOM"+space+"H2"+space+"1"+space
#		+ simType + space
#		+ simType1 + space
#		+ molecule + space
#		+ rotor + space
#		+ str(param) + space
#		+ variable_name + space
#		#+ " -C --RESTART"+space
#		#+ " -NR 20000"+space
#	)

    # print(cmd_run)
    # os.system(cmd_run)
    #call(["rm", file_name3])
