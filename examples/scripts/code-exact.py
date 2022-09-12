import os
import numpy as np
import support

# Informations about the system
rotor_name="HF"
numb_molecule=2
rlist = np.arange(3.0, 10.01, 0.2, dtype=float)
dipole_moment = 1.827
l_max=10
l_total_max=l_max
#
blank_space = " "

for rpt in rlist:
	rpt_value = "{:3.2f}".format(rpt)
	command_line = (
		+ str(rfactor)
		+ "-rpt"
		+ rpt_value
		+ " -N "
		+ str(numb_molecule)
		+ "--l-max "
		+ str(lmax)
		+ "--l-total-max"
		+ str(l_total_max)
		+ " --rotor"
		+ rotor_name
	)
	os.system(command_line)
