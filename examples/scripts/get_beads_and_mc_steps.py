import numpy as np

class GetBeadStepLevel:
	def __init__(self, molecule_rot1, variable1, simulation1):
		self.molecule_rot = molecule_rot1
		self.variable = variable1
		self.simulation = simulation1

		if (self.simulation == "PIGS"):
			if (self.variable == "beta"):
				self.beads=np.arange(10,130,10,dtype=int)
				if (self.molecule_rot == "HF"):
					self.level_trans= np.array([1 for i in range(100)],dtype=int)
					self.step_trans = np.array([0.3 for i in range(100)],dtype=float)
					self.step_rot   = np.array([2.0 for i in range(100)],dtype=float)
			if (self.variable == "tau"):
				self.beads=np.arange(10,110,10,dtype=int)
				if (self.molecule_rot == "HF"):
					self.step_trans = {10:2.0, 20:2.0, 30:2.0, 40:2.0, 50:1.5, 60:1.5,
									   70:1.5, 80:1.5, 90:1.0, 100:1.0, 110:1.0, 120:1.0}
					self.step_rot   = {10:2.0, 20:2.0, 30:2.0, 40:2.0, 50:1.5, 60:1.5,
									   70:1.5, 80:1.5, 90:1.0, 100:1.0, 110:1.0, 120:1.0}
					self.level_trans= np.array([1 for i in range(100)],dtype=int)

if __name__ == '__main__':
	molecule_rot1="HF"
	variable1="tau"
	simulation1="PIGS"
	print(GetBeadStepLevel(molecule_rot1,variable1,simulation1).beads[0])
	print(GetBeadStepLevel(molecule_rot1,variable1,simulation1).step_rot[10])
