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
			list_nb = [4,8,12,16,20,24,28,32,36,40]
		if (variableName == "beta"):
			list_nb = [4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40]

	if (TypeCal == "PIMC"):
		if (variableName == "tau"):
			list_nb  = [4, 8, 16, 32, 64]
		if (variableName == "beta"):
			list_nb  = [4, 8, 16, 32, 64, 128]

	if (TypeCal == "PIGS"):
		if (variableName == "tau"):
			list_nb  = [4, 8, 16, 32, 64, 128, 256]
		if (variableName == "beta"):
			list_nb  = [4, 8, 16, 32, 64, 128, 256]

	return list_nb

class GetStepAndLevel:
	def __init__(self, molecule_rot1, variableName1, TypeCal1):
		self.molecule_rot = molecule_rot1
		self.variableName = variableName1
		self.TypeCal = TypeCal1
		if ((self.variableName == "tau") and (self.TypeCal == "ENT")):
			if (self.molecule_rot == "H2"):
				self.step_trans = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				#step           = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
				#step           = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
				self.step       = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6
				self.level      = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

			if (self.molecule_rot == "HF"):
				self.step_trans = [1.0,1.0,0.1,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				#self.step       = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.9, 1.7, 1.7, 1.7, 1.6, 1.5] #list_nb = [4, 10, 14, 20, 24, 30, 34, 40, 44, 50, 54, 60] beta = 0.2; g=[0 - 8.0]
				self.step       = [0.5, 0.8, 1.2, 1.2, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1] #list_nb = [4, 10, 14, 20, 24, 30, 34, 40, 44, 50, 54, 60] beta = 0.2; g = [9.0 - 20.0]
				#self.step       = [0.3, 0.5, 0.5, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8] #list_nb = [4, 10, 14, 20, 24, 30, 34, 40, 44, 50, 54, 60] beta = 0.2; g = [21.0 - 30.0]
				#self.step       = [2.0, 2.0, 1.9, 1.7, 1.6] #list_nb = [4, 14, 24, 34, 44, 54] beta = 0.2
				#self.step       = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0] #list_nb = [4, 10, 14, 20, 24, 30] beta = 0.2
				#self.step       = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0] #list_nb = [5, 10, 20, 30, 40, 50, 60] beta = 0.32
				self.level      = [1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

		if ((self.variableName == "tau") and (self.TypeCal == "PIGS")):
			if (self.molecule_rot == "H2"):
				self.step_trans = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				#step           = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
				#step           = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
				self.step       = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6
				self.level      = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

			if (self.molecule_rot == "HF"):
				self.step_trans = [1.0,1.0,0.1,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				#self.step       = [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.9, 1.7, 1.7, 1.7, 1.6, 1.5] #list_nb = [4, 10, 14, 20, 24, 30, 34, 40, 44, 50, 54, 60] beta = 0.2; g=[0 - 8.0]
				self.step       = [2.0, 2.0, 2.0, 2.0, 1.5] #list_nb = [4, 8, 16, 32, 64] beta = 0.05; 
				self.level      = [1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

			if (self.molecule_rot == "H2O"):
				self.step_trans = [1.0,1.0,0.1,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				#self.step       = [1.8, 1.8, 0.7, 0.4, 0.25, 0.15] #list_nb = [4, 8, 16, 32, 64, 128] beta = 0.32; 
				#self.step       = [2.0, 2.0, 2.0, 0.4, 0.25, 0.15, 0.1, 0.08] #list_nb = [4, 8, 16, 32, 64, 128, 256, 512] beta = 0.256 for 1 p-H2O V=-20*cos(theta) 
				#self.step       = [2.0, 2.0, 0.3, 0.2, 0.15, 0.12, 0.08, 0.05] #list_nb = [4, 8, 16, 32, 64, 128, 256, 512] beta = 0.128 for 1 p-H2O V=-20*cos(theta)
				self.step       = [2.0, 2.0, 2.0, 0.4, 0.25, 0.15, 0.1, 0.08] #list_nb = [4, 8, 16, 32, 64, 128, 256, 512] beta = 0.256 for 2 p-H2O 
				self.level      = [1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]


		if ((self.variableName == "tau") and (self.TypeCal == "PIMC")):
			if (self.molecule_rot == "H2"):
				self.step_trans = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				#step           = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
				#step           = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
				self.step       = [1.5,2.0,2.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6
				self.level      = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

			if (self.molecule_rot == "HF"):
				self.step_trans = [1.0,1.0,0.1,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				self.step       = [2.0, 1.8, 1.4, 1.0, 0.6] #list_nb = [4, 8, 16, 32, 64] beta = 0.05 Kelvin^-1; 
				self.level      = [1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

		if self.variableName == "beta":
			if (self.molecule_rot == "H2"):
				self.step_trans  = [0.3 for i in range(100)]
				self.step        = [1.6 for i in range(100)]  
				self.level       = [1   for i in range(100)]

			if (self.molecule_rot == "HF"):
				self.step_trans  = [0.3 for i in range(100)]
				self.step        = [2.0 for i in range(100)]  
				self.level       = [1   for i in range(100)]

			if (self.molecule_rot == "H2O"):
				self.step_trans  = [0.3 for i in range(100)]
				#self.step        = [0.4 for i in range(100)]  # for tau = 0.01
				#self.step        = [0.15 for i in range(100)]  # for tau = 0.002 for 1 p-H2O
				self.step        = [0.2 for i in range(100)]  # for tau = 0.004 for 2 p-H2O
				self.level       = [1   for i in range(100)]
