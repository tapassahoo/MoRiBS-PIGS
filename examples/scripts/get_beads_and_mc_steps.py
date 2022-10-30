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
				beads1=np.array([4, 6, 8, 14, 24])
				beads2=np.arange(10,110,10,dtype=int)
				self.beads=np.sort(np.append(beads1, beads2))
				if (self.molecule == "HF"):
					self.step_com  = {4:1.5, 6:1.5, 8:1.5, 10:1.5, 14:1.5, 20:1.5, 24:1.5, 30:1.5, 40:1.5, 50:1.5, 60:1.5, 70:1.5, 80:1.5, 90:1.5, 100:1.5}
					self.step_rot  = {4:0.6, 6:0.6, 8:0.6, 10:0.6, 14:0.6, 20:0.6, 24:0.6, 30:0.6, 40:0.6, 50:0.6, 60:0.6, 70:0.6, 80:0.6, 90:0.6, 100:0.6}
					self.level_com = {4:1, 6:1, 8:1, 10:1, 14:1, 20:1, 24:1, 30:1, 40:1, 50:1, 60:1, 70:1, 80:1, 90:1, 100:1}
			if (self.parameter_name == "beta"):
				self.beads=np.arange(10,110,10,dtype=int)
				if (self.molecule == "HF"):
					self.step_com  = {10:1.5, 20:1.5, 30:1.5, 40:1.5, 50:1.5, 60:1.5, 70:1.5, 80:1.5, 90:1.5, 100:1.5}
					self.step_rot  = {10:0.7, 20:0.8, 30:0.9, 40:0.9, 50:0.9, 60:0.9, 70:0.9, 80:0.9, 90:0.9, 100:0.9}
					self.level_com = {10:1, 20:1, 30:1, 40:1, 50:1, 60:1, 70:1, 80:1, 90:1, 100:1}
