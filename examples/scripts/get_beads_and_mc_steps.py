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
	def __init__(self, molecule1, parameter1, simulation1):
		self.molecule = molecule1
		self.parameter = parameter1
		self.simulation = simulation1

		if (self.simulation == "PIGS"):
			if (self.parameter == "tau"):
				self.beads=np.arange(10,130,10,dtype=int)
				if (self.molecule == "HF"):
					self.level_com = np.array([1 for i in range(100)],dtype=int)
					self.step_com = np.array([0.3 for i in range(100)],dtype=float)
					self.step_rot   = np.array([2.0 for i in range(100)],dtype=float)
			if (self.parameter == "beta"):
				self.beads=np.arange(10,110,10,dtype=int)
				if (self.molecule == "HF"):
					self.step_com = {10:2.0, 20:2.0, 30:2.0, 40:2.0, 50:1.5, 60:1.5,
									   70:1.5, 80:1.5, 90:1.0, 100:1.0, 110:1.0, 120:1.0}
					self.step_rot   = {10:2.0, 20:2.0, 30:2.0, 40:2.0, 50:1.5, 60:1.5,
									   70:1.5, 80:1.5, 90:1.0, 100:1.0, 110:1.0, 120:1.0}
					self.level_com = {10:1, 20:1, 30:1, 40:1, 50:1, 60:1,
									   70:1, 80:1, 90:1, 100:1, 110:1.0, 120:1.0}
