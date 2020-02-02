import os
from subprocess import call

import mypkg.pkgAsymrho

module_path = mypkg.pkgAsymrho.__file__
module_path = module_path.replace("__init__.py", "")
space = " "
file_sim = module_path + "script_submission_asymrho.py"

# Informations about the system
simType = "S"
rotor = "H2O"
param = "tau"
value = 0.001
#param = "beta"
#value = 0.2
#nbeads = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
nbeads =  [160]
#nbeads = [120, 140, 160, 180, 200]
#nbeads += [200, 250, 300, 350]
#nbeads += [400]
Jmax = 70
spin = 0#int(-1)

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
