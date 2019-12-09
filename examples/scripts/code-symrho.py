import os
from subprocess import call

import mypkg.pkgSymrho

module_path = mypkg.pkgSymrho.__file__
module_path = module_path.replace("__init__.py", "")
space = " "
file_sim = module_path + "script_submission_symrho.py"

# Informations about the system
simType = "A"
rotor = "CH3F"
param = "tau"
value = 0.001
#param = "beta"
#value = 0.025
nbeads = [4, 6, 8, 10, 12, 14, 16, 18]#, 20]#, 24, 30, 34, 40, 44, 50, 60, 70, 80, 90, 100]
#nbeads = [20]
Jmax = 70
spin = 3

# Run the "cmd_run" in python interpreter
for bead in nbeads:
    cmd_run = (
        "python"
        + space
        + file_sim
        + space
        + rotor
        + space
        + param
        + space
        + str(value)
        + space
        + "-P"
        + space
        + str(bead)
        + space
        + "-iodevn"
        + space
        + str(spin)
        + space
        + "-J"
        + space
        + str(Jmax)
        + space
        + "-job"
        + space
        + simType
        + space
    )

    os.system(cmd_run)
