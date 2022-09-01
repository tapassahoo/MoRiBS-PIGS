import numpy as np
"""
Users can modify this code to run the executable, namely "pimc" 
generated in the MoRiBS-PIGS directory. 

This code provides Monte Carlo steps and levels.

Note that always even beads are used here for both the finite 
temperature and ground state simulations.
"""

class GetBeadStepLevel:
	"""
	inputs are molecule_rot1, variable1, simulation1
	outputs are self.beads, self.level_trans, self.step_trans, self.step_rot
	"""
	def __init__(self, molecule1, parameter_name1, simulation1):
		self.molecule = molecule1
		self.parameter_name = parameter_name1
		self.simulation = simulation1

		if (self.simulation == "PIGS"):
			if (self.parameter_name == "tau"):
				self.beads=np.arange(10,110,10,dtype=int)
				if (self.molecule == "HF"):
					self.step_com  = {10:1.5, 20:1.5, 30:1.5, 40:1.5, 50:1.5, 60:1.5,
									   70:1.5, 80:1.5, 90:1.5, 100:1.5, 110:1.5, 120:1.5}
					self.step_rot  = {10:1.5, 20:1.5, 30:1.5, 40:1.5, 50:1.5, 60:1.5,
									   70:1.5, 80:1.5, 90:1.5, 100:1.5, 110:1.5, 120:1.5}
					self.level_com = {10:1, 20:1, 30:1, 40:1, 50:1, 60:1,
									   70:1, 80:1, 90:1, 100:1, 110:1, 120:1}
			if (self.parameter_name == "beta"):
				self.beads=np.arange(10,110,10,dtype=int)
				if (self.molecule == "HF"):
					self.step_com  = {10:1.5, 20:1.5, 30:1.5, 40:1.5, 50:1.5, 60:1.5,
									   70:1.5, 80:1.5, 90:1.5, 100:1.5, 110:1.5, 120:1.5}
					self.step_rot  = {10:2.0, 20:2.0, 30:1.8, 40:1.4, 50:1.2, 60:1.2,
									   70:1.1, 80:1.0, 90:1.0, 100:0.8}
					self.level_com = {10:1, 20:1, 30:1, 40:1, 50:1, 60:1,
									   70:1, 80:1, 90:1, 100:1, 110:1, 120:1}
