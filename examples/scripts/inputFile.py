import time
from subprocess import call
import os
import decimal
import numpy as np
from numpy import *
import math

def Getbeads(TypeCal, variableName):
	if (TypeCal == "ENT"):
		if (variableName == "tau"):
			list_nb = [10, 20, 30, 40, 50, 60, 80, 100, 120]#, 140]#, 160, 180, 200]
		if (variableName == "beta"):
			list_nb = [4, 6, 8, 10, 14, 20, 24, 30, 40, 50, 60, 70, 80, 90, 100]

	if (TypeCal == "PIMC"):
		if (variableName == "tau"):
			list_nb  = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]#, 120, 140, 160, 180, 200]

	if (TypeCal == "PIGS"):
		if (variableName == "tau"):
			#list_nb  = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]#, 120, 140, 160, 180, 200]
			list_nb  = [50,100]
		if (variableName == "beta"):
			list_nb  = [4, 6, 8, 10, 14, 20, 30, 40, 50, 60, 70, 80, 90, 100]

	return list_nb

class GetStepAndLevel:
	def __init__(self, molecule_rot1, variableName1, TypeCal1):
		self.molecule_rot = molecule_rot1
		self.variableName = variableName1
		self.TypeCal = TypeCal1
		if ((self.variableName == "tau") and (self.TypeCal == "ENT")):
			if (self.molecule_rot == "H2"):
				self.step_trans = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				self.step       = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6
				self.level      = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

			if (self.molecule_rot == "HF"):
				self.step_trans = [1.0,1.0,0.1,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				self.step       = [2.0, 2.0, 2.0, 2.0, 1.5, 1.5, 1.5, 1.2, 1.1, 1.1] 
				self.level      = [1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

			if (self.molecule_rot == "H2O"):
				self.step_trans = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
				self.step       = [0.5,0.25,0.2,0.16,0.14,0.13,0.12,0.11,0.1,0.08,0.08,0.08,0.08]
				self.level      = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

		if ((self.variableName == "tau") and (self.TypeCal == "PIGS")):
			if (self.molecule_rot == "H2"):
				self.step_trans = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				self.step       = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6
				self.level      = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

			if (self.molecule_rot == "HF"):
				self.step_trans = [1.0,1.0,0.1,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				self.step       = [2.0, 2.0, 2.0, 2.0, 1.5] #list_nb = [4, 8, 16, 32, 64] beta = 0.05; 
				self.level      = [1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

			if (self.molecule_rot == "H2O"):
				self.step_trans = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]
				self.level      = [  2,  2,  3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  5,  5,  5, 5, 5,5,5,5,5,5,5,5,5,5,5,5,5,5]
				#self.step       = [0.13,0.13,0.13,0.13,0.13,0.13,0.12,0.12,0.12,0.11,0.1,0.08,0.08,0.08,0.07,0.06,0.06,0.06,0.05] 
				self.step       = [0.25,0.25] 

			if (self.molecule_rot == "CH3F"):
				self.step_trans1= [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.8,0.8,0.8,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
				self.level1     = [  1,  1,  2,  2,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  3,  3,  3,  3,  3,  3,  3,  3]
				self.step_trans = [0.8,0.6,0.5,0.5,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3]
				self.level      = [  2,  2,  3,  3,  4,  4,  4,  4,  4,  4,  5,  5,  5,  5,  5,  4,  4,  4,  4,  4,  4,  4,  4]
				self.step       = [0.2,0.15,0.12,0.10,0.10,0.10,0.08,0.07,0.06,0.06,0.06,0.05,0.04,0.04,0.04,0.04]


		if ((self.variableName == "tau") and (self.TypeCal == "PIMC")):
			if (self.molecule_rot == "H2"):
				self.step_trans = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				self.step       = [1.5,2.0,2.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6
				self.level      = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

			if (self.molecule_rot == "HF"):
				self.step_trans = [1.0,1.0,0.1,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				self.step       = [2.0, 1.8, 1.4, 1.0, 0.6] #list_nb = [4, 8, 16, 32, 64] beta = 0.05 Kelvin^-1; 
				self.level      = [1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

		if (self.molecule_rot == "H2O"):
				self.step_trans = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]
				self.level      = [  2,  2,  3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  5,  5,  5, 5, 5,5,5,5,5,5,5,5,5,5,5,5,5,5]
				self.step       = [0.1,0.08,0.08,0.06,0.06,0.06,0.04,0.04,0.04,0.04,0.1,0.08,0.08,0.08,0.07,0.06,0.06,0.06,0.05] 

		if ((self.variableName == "beta") and (self.TypeCal == "PIGS")):
			if (self.molecule_rot == "H2"):
				self.step_trans  = [0.3 for i in range(100)]
				self.step        = [1.6 for i in range(100)]  
				self.level       = [1   for i in range(100)]

			if (self.molecule_rot == "HF"):
				self.step_trans  = [0.3 for i in range(100)]
				self.step        = [2.0 for i in range(100)]  
				self.level       = [1   for i in range(100)]

			if (self.molecule_rot == "H2O"):
				self.step_trans  = [2.0, 0.6, 0.5, 0.5, 0.4, 0.4, 0.3, 0.2, 0.2, 0.2, 0.2, 0.15, 0.12, 0.12]
				self.level       = [1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
				self.step        = [0.15 for i in range(20)]

			if (self.molecule_rot == "CH3F"):
				self.step_trans1 = [0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]
				self.level1      = [  1,   2,   2,   2,   2,   2,   2,   2,   2,   2, 2, 2]
				self.step_trans  = [0.4, 0.6, 0.5, 0.4, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3, 0.3]
				self.level       = [  1,   3,   3,   3,   3,   3,   3,   3,   3,   3, 2, 2]
				self.step        = [0.15 for i in range(20)]
				#self.step_trans1 = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05, 0.2, 0.2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
				#self.level1      = [   3,    4,    4,    5,    5,    5,    6,   6,   6,   6,   6,   6,   6,   6,   6]
				#self.step_trans  = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
				#self.level       = [   3,    4,    4,    5,    5,    5,    6,   6,   6,   6,   6,   6,   6,   6,   6]
				#self.step        = [0.035 for i in range(40)]

		if ((self.variableName == "beta") and (self.TypeCal == "PIMC")):
			if (self.molecule_rot == "H2O"):
				self.step_trans  = [0.3 for i in range(100)]
				self.step        = [0.1 for i in range(100)]
				self.level       = [1 for i in range(100)]

		if ((self.variableName == "beta") and (self.TypeCal == "ENT")):
			if (self.molecule_rot == "H2O"):
				self.step_trans  = [0.3 for i in range(100)]
				self.step        = [0.2 for i in range(100)]
				self.level       = [1 for i in range(100)]
