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
	list_nb = [4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,44]
	return list_nb

class GetStepAndLevel:
	def __init__(self, molecule_rot1, variableName1):
		self.step_trans  = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.10,1.20,1.30,1.40,1.50]
		self.step        = [1.5,3.0,3.0,2.0,1.0,0.7,0.5,2.5,2.02]
		self.level       = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
