import os
from subprocess import call

import mypkg.pkgAsymrho

module_path = mypkg.pkgAsymrho.__file__
module_path = module_path.replace("__init__.py", "")
space = " "
file_sim = module_path + "script_submission_asymrho.py"

# Informations about the system
simType = "A" # "S" - submits 181 jobs in a queue for a particular temperature (or beta) for all theta values from 0-180 degrees. "A"- executes compile.x to collect all the values of all thetas." 
rotor = "H2O"
#param = "tau" 
#value = 0.0005
param = "beta" # beta=1/T K^-1
value = 0.2
nbeads = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200]
Jmax = 66
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
