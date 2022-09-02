import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
import math
from pathlib import Path
import get_beads_and_mc_steps as mc
import mypkg.asymrho_runner.support as asym
import mypkg.symrho_runner.support as sym

class get_analysis_file_name:
	def __init__(
			self,
			method1,
			method2,
			molecule_rot1,
			translational_move1,
			rotational_move1,
			rpt_val1,
			gfacor1,
			dipole_moment1,
			parameter_name1,
			parameter_value1,
			numb_block1,
			numb_pass1,
			numb_molecule1,
			molecular_system1,
			ent_method1,
			preskip1,
			postskip1,
			extra_file_name1,
			input_dir1,
			particle_a1,
			ent_algorithm):
		self.method = method1
		self.rotor = molecule_rot1
		self.translational_move = translational_move1
		self.rotational_move = rotational_move1
		self.rpt_val = rpt_val1
		self.dipole_moment = dipole_moment1
		self.gfactor = gfacor1
		self.parameter_value = parameter_value1
		self.parameter_name = parameter_name1
		self.numb_block = numb_block1
		self.numb_pass = numb_pass1
		self.numb_molecule = numb_molecule1
		self.molecular_system = molecular_system1
		self.ent_method = ent_method1
		self.preskip = preskip1
		self.postskip = postskip1
		self.extra_file_name = extra_file_name1
		self.input_dir = input_dir1
		self.particle_a = particle_a1

		if (self.method == "ENT"):
			front_layer = "ENT-" + self.extra_file_name
			add1 = "-ParticleA" + str(self.particle_a)
			add2 = "-" + self.ent_method + "-" + ent_algorithm
		else:
			front_layer = self.method + "-" + self.extra_file_name
			add1 = ""
			add2 = ""

		if ((self.translational_move) and (self.rotational_move)):
			front_layer += "TransAndRotDOFs-"
		if ((self.translational_move == False) and (self.rotational_move)):
			front_layer += "RotDOFs-"
		if ((self.translational_move) and (self.rotational_move == False)):
			front_layer += "TransDOFs-"

		if (self.rpt_val >= 0.0):
			name_rpt = "Rpt" + str(self.rpt_val) + "Angstrom-"
		else:
			name_rpt = ""

		if (self.dipole_moment >= 0.0):
			if (self.numb_molecule > 1):
				name_dipole_moment = "Dipole-Moment" + \
					str(self.dipole_moment) + "Debye-"
			else:
				name_dipole_moment = "Field" + \
					str(self.dipole_moment) + "Kinv-"
		else:
			name_dipole_moment = ""

		if (self.gfactor >= 0.0):
			name_gfactor = "gFactor" + str(self.gfactor) + "-"
		else:
			name_gfactor = ""

		if (self.parameter_name == "beta"):
			variable_name = "tau"
		if (self.parameter_name == "tau"):
			variable_name = "beta"

		name_layer1 = "vs-" + str(variable_name) + "-fixed-" + self.parameter_name + \
			str(self.parameter_value) + "Kinv-Blocks" + str(self.numb_block)
		name_layer1 += "-Passes" + str(self.numb_pass) + "-System" + str(self.numb_molecule) + str(
			self.molecular_system) + add1 + add2 + "-preskip" + str(self.preskip) + "-postskip" + str(self.postskip)

		file_output1 = front_layer + name_rpt + \
			name_dipole_moment + name_gfactor + "Energy-"
		file_output2 = front_layer + name_rpt + \
			name_dipole_moment + name_gfactor + "correlation-"
		file_output3 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "total-correlation-function-"
		file_output4 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "X-component-correlation-function-"
		file_output5 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "Y-component-correlation-function-"
		file_output6 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "Z-component-correlation-function-"
		file_output7 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "XandY-component-correlation-function-"
		file_output8 = front_layer + name_rpt + \
			name_dipole_moment + name_gfactor + "Entropy-"

		self.save_file_energy = self.input_dir + file_output1 + name_layer1 + ".txt"
		self.save_file_correlation = self.input_dir + file_output2 + name_layer1 + ".txt"
		self.SaveEntropy = self.input_dir + file_output8 + name_layer1 + ".txt"

		if (method2 == False):
			if os.path.exists(self.SaveEntropy):
				os.remove(self.SaveEntropy)
			if os.path.exists(self.save_file_energy):
				os.remove(self.save_file_energy)
			if os.path.exists(self.save_file_correlation):
				os.remove(self.save_file_correlation)

		if (self.method != "ENT"):
			print("****************** Important message *******************")
			print(f'Name of the file where the potential and the total energies computed for several {variable_name} values are stored are given below:')
			print(file_output1 + name_layer1 + ".txt")
			print("")
			print("Final analyzed results are stored in - ")
			print(self.input_dir)
			print("")
			print("********************************************************")

		if (self.method == "ENT"):
			name_layer1RT = "vs-" + str(self.variable_name) + "-fixed-" + self.parameter_name + str(
				self.parameter) + "Kinv-Blocks" + str(self.numb_block)
			name_layer1RT += "-Passes" + str(self.numb_pass) + "-System" + str(self.numb_molecule) + str(
				self.molecule) + add1 + add2 + "-preskip" + str(self.preskip) + "-postskip" + str(self.postskip)

			self.SaveEntropyRT = self.input_dir + file_output8 + name_layer1RT + ".txt"
			if (ent_algorithm != "WR"):
				# print(self.SaveEntropyRT)
				print(self.input_dir)
				print("")
				print("Final results - Entropy vs " + str(self.variable_name))
				print(self.SaveEntropy)

				print(
					"#------------------------------------------------------------------------#")


class GetFileNamePlot:
	def __init__(
			self,
			method1,
			molecule_rot1,
			translational_move1,
			rotational_move1,
			variable_name1,
			rpt_val1,
			gfacor1,
			dipolemoment1,
			parameter_name1,
			parameter1,
			numbblocks1,
			numbpass1,
			numbmolecules1,
			molecule1,
			ent_method1,
			preskip1,
			postskip1,
			extra1,
			input_dir1,
			particleA1,
			var1):
		self.method = method1
		self.rotor = molecule_rot1
		self.translational_move = translational_move1
		self.rotational_move = rotational_move1
		self.variable_name = variable_name1
		self.var = var1
		self.rpt_val = rpt_val1
		self.dipole_moment = dipolemoment1
		self.gfactor = gfacor1
		self.parameter = parameter1
		self.parameter_name = parameter_name1
		self.numb_block = numbblocks1
		self.numb_pass = numbpass1
		self.numb_molecule = numbmolecules1
		self.molecule = molecule1
		self.ent_method = ent_method1
		self.preskip = preskip1
		self.postskip = postskip1
		self.extra = extra1
		self.input_dir = input_dir1
		self.particle_a = particleA1

		if (self.method == "ENT"):
			front_layer = "ENT-" + self.extra
			add1 = "-ParticleA" + str(self.particle_a)
			if (self.ent_method):
				add2 = "-" + self.ent_method
			else:
				add2 = ""
		else:
			front_layer = self.method + "-" + self.extra
			add1 = ""
			add2 = ""

		if ((self.translational_move) and (self.rotational_move)):
			front_layer += "TransAndRotDOFs-"
		if ((self.translational_move == False) and (self.rotational_move)):
			front_layer += "RotDOFs-"
		if ((self.translational_move) and (self.rotational_move == False)):
			front_layer += "TransDOFs-"

		if (self.rpt_val >= 0.0):
			name_rpt = "rpt_val" + str(self.rpt_val) + "Angstrom-"
		else:
			name_rpt = ""

		if (self.dipole_moment >= 0.0):
			name_dipole_moment = "dipole_moment" + \
				str(self.dipole_moment) + "Debye-"
		else:
			name_dipole_moment = ""

		if (self.gfactor >= 0.0):
			name_gfactor = "gFactor" + str(self.gfactor) + "-"
		else:
			name_gfactor = ""

		name_layer1 = "vs-" + str(self.variable_name) + "-fixed-" + self.parameter_name + \
			str(self.parameter) + "Kinv-Blocks" + str(self.numb_block)
		name_layer1 += "-Passes" + str(self.numb_pass) + "-System" + str(self.numb_molecule) + str(
			self.molecule) + add1 + "-preskip" + str(self.preskip) + "-postskip" + str(self.postskip) + add2

		file_output1 = front_layer + name_rpt + \
			name_dipole_moment + name_gfactor + "Energy-"
		file_output2 = front_layer + name_rpt + \
			name_dipole_moment + name_gfactor + "correlation-"
		file_output3 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "total-correlation-function-"
		file_output4 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "X-component-correlation-function-"
		file_output5 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "Y-component-correlation-function-"
		file_output6 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "Z-component-correlation-function-"
		file_output7 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "XandY-component-correlation-function-"
		file_output8 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "Chemical-Potential-"
		file_output9 = front_layer + name_rpt + \
			name_dipole_moment + name_gfactor + "Entropy-"
		file_output10 = front_layer + name_dipole_moment + name_gfactor + "Energy-vs-R-"
		file_output11 = front_layer + name_dipole_moment + \
			name_gfactor + "Order-parameters-vs-R-"
		name_layer1FitvsR = "fixed-" + self.parameter_name + \
			str(self.parameter) + "Kinv-Blocks" + str(self.numb_block)
		name_layer1FitvsR += "-Passes" + str(self.numb_pass) + "-System" + str(self.numb_molecule) + str(
			self.molecule) + add1 + "-preskip" + str(self.preskip) + "-postskip" + str(self.postskip) + add2

		self.save_file_energy = self.input_dir + file_output1 + name_layer1
		self.save_file_correlation = self.input_dir + file_output2 + name_layer1
		self.SaveChemPot = self.input_dir + file_output8 + name_layer1
		self.SaveEntropy = self.input_dir + file_output9 + name_layer1
		self.save_file_energyFitvsR = self.input_dir + file_output10 + name_layer1FitvsR
		self.save_file_correlationFitvsR = self.input_dir + file_output11 + name_layer1FitvsR

#---------------------------------------------------------------------------#
#	   special cases														   #
#---------------------------------------------------------------------------#
		'''
				name_layer1CP		= "vs-number-of-"+str(self.molecule)+"-fixed-"+self.parameter_name+str(self.parameter)+"Kinv-Blocks"+str(self.numb_block)
				name_layer1CP	   += "-Passes"+str(self.numb_pass)+add1+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)
				name_layer1CONV	  = "vs-beta-and-tau-Blocks"+str(self.numb_block)
				name_layer1CONV	 += "-Passes"+str(self.numb_pass)+"-System"+str(self.numb_molecule)+str(self.molecule)+add1+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)
				self.SaveChemPot	  = self.input_dir+file_output8+name_layer1CP
				self.SaveEntropyCONV  = self.input_dir+"/ResultsOfPIGSENT/"+file_output1+name_layer1CONV+"-"+self.ent_method

				'''
#
		name_layer1GFAC = "vs-gFactor-of-" + str(self.molecule) + "-fixed-" + self.parameter_name + str(
			self.parameter) + "Kinv-numb_bead" + str(self.var) + "-Blocks" + str(self.numb_block)
		name_layer1GFAC += "-Passes" + \
			str(self.numb_pass) + "-preskip" + str(self.preskip) + \
			"-postskip" + str(self.postskip) + add2
		name_layer1GFACFit = "vs-gFactor-of-" + \
			str(self.molecule) + "-fixed-" + self.parameter_name + \
			str(self.parameter) + "Kinv-Blocks" + str(self.numb_block)
		name_layer1GFACFit += "-Passes" + \
			str(self.numb_pass) + "-preskip" + str(self.preskip) + \
			"-postskip" + str(self.postskip) + add2 + "-Fit"
#
		name_layer1RFAC = "vs-RFactor-of-" + str(self.molecule) + "-fixed-" + self.parameter_name + str(
			self.parameter) + "Kinv-numb_bead" + str(self.var) + "-Blocks" + str(self.numb_block)
		name_layer1RFAC += "-Passes" + \
			str(self.numb_pass) + "-preskip" + str(self.preskip) + \
			"-postskip" + str(self.postskip) + add2
		name_layer1MM = "vs-" + \
			str(self.variable_name) + "-fixed-" + \
			self.parameter_name + str(self.parameter) + "Kinv"
		name_layer1MM += "-System" + \
			str(self.numb_molecule) + str(self.molecule) + add1
		name_layer1COMBO = "vs-gFactor-fixed-" + self.parameter_name + \
			str(self.parameter) + "Kinv-Blocks" + str(self.numb_block)
		name_layer1COMBO += "-Passes" + str(self.numb_pass) + "-System" + str(self.numb_molecule) + str(
			self.molecule) + add1 + "-preskip" + str(self.preskip) + "-postskip" + str(self.postskip)
		name_layer1ED = "of-System" + \
			str(self.numb_molecule) + str(self.molecule) + add1
#
		self.SaveEntropyGFAC = self.input_dir + front_layer + \
			name_rpt + "Entropy-" + name_layer1GFAC
		self.SaveEntropyRFAC = self.input_dir + front_layer + \
			name_rpt + "Entropy-" + name_layer1RFAC
		self.save_file_energyGFAC = self.input_dir + front_layer + \
			name_rpt + "Energy-" + name_layer1GFAC
		self.save_file_energyRFAC = self.input_dir + front_layer + \
			name_rpt + "Energy-" + name_layer1RFAC
		self.save_file_energyED = self.input_dir + file_output1 + name_layer1ED + add2 + "-ED"
		self.SaveEntropyED = self.input_dir + file_output9 + name_layer1ED + add2 + "-ED"
		self.save_file_correlationGFAC = self.input_dir + front_layer + \
			name_rpt + "correlation-" + name_layer1GFAC
		self.save_file_correlationRFAC = self.input_dir + front_layer + \
			name_rpt + "correlation-" + name_layer1RFAC
		self.save_file_energyMM = self.input_dir + file_output1 + name_layer1MM + add2 + "-MM"
		self.SaveEntropyMM = self.input_dir + file_output9 + name_layer1MM + add2 + "-MM"
		self.SaveEntropyCOMBO = self.input_dir + front_layer + \
			name_rpt + "Entropy-" + name_layer1GFACFit + "-COMBINE"
		self.SaveEntropyGFACFit = self.input_dir + front_layer + \
			name_rpt + "Entropy-" + name_layer1GFACFit
		self.save_file_energyGFACFit = self.input_dir + front_layer + \
			name_rpt + "Energy-" + name_layer1GFACFit

		if (self.method == "ENT"):
			name_layer1RT = "vs-" + str(self.variable_name) + "-fixed-" + self.parameter_name + str(
				self.parameter) + "Kinv-Blocks" + str(self.numb_block)
			name_layer1RT += "-Passes" + str(self.numb_pass) + "-System" + str(self.numb_molecule) + str(
				self.molecule) + "-preskip" + str(self.preskip) + "-postskip" + str(self.postskip) + add2

			self.SaveEntropyRT = self.input_dir + file_output9 + name_layer1RT
