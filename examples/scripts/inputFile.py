#!/usr/bin/python

import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import math

def Getbeads(TypeCal, variableName):
	if (variableName == "beta"):
		if (TypeCal == "ENT"):
			list_nb = [4,8,12,16,20,24,28,32,36,40]
		else:
			list_nb = [4,8,12,16,20,24,28,32,36,40]

	if (variableName == "tau"):
		list_nb  = [4, 10, 20, 40, 60]
		#list_nb  = [24, 30]

	if (TypeCal != "ENT"):
		print(" ")
		print("====================================== ")
		print(" List of quantum beads ")
		print(list_nb)
		print(" ")
		print("====================================== ")
		print(" ")
	return list_nb

class GetStepAndLevel:
	def __init__(self, molecule_rot1, variableName1):
		self.molecule_rot = molecule_rot1
		self.variableName = variableName1
		if self.variableName == "tau":
			if (self.molecule_rot == "H2"):
				self.step_trans = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				#step           = [1.5,3.0,3.0,3.0,3.0,2.6,2.3,2.5,2.02] #temp 10K             #change param6
				#step           = [1.5,3.0,3.0,2.5,1.5,1.0,0.7,2.5,2.02] #temp 50K             #change param6
				self.step       = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02] #temp 100K            #change param6
				self.level      = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

			if (self.molecule_rot == "HF"):
				self.step_trans = [1.0,1.0,0.1,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
				self.step       = [2.0, 2.0, 2.0, 1.7, 1.5] #list_nb = [4, 10, 20, 40, 60] beta = 0.2
				#self.step       = [2.0, 2.0, 2.0]            #list_nb = [14, 24, 30] beta = 0.2
				self.level      = [1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

		if self.variableName == "beta":
			if (self.molecule_rot == "H2"):
				self.step_trans  = [0.3 for i in range(1000)]
				self.step        = [1.6 for i in range(1000)]  
				self.level       = [1   for i in range(1000)]

			if (self.molecule_rot == "HF"):
				self.step_trans  = [0.3 for i in range(1000)]
				self.step        = [2.0 for i in range(1000)]  
				self.level       = [1   for i in range(1000)]
