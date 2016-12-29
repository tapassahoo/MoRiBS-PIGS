#!/usr/bin/python
 
####
#
# Set bconst
#
####
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *

def dropzeros(number):
    mynum          = decimal.Decimal(number).normalize()
    # e.g 22000 --> Decimal('2.2E+4')
    return mynum.__trunc__() if not mynum % 1 else float(mynum)

def bconstant():
	'''
	This function calculates rotational Bconstant for linear rotor
	'''
	autocminverse  = 2.1947463137e+5
	energyj0       = -36117.5942855
	energyj1       = -35999.1009407
	bconst         = 0.5*(energyj1-energyj0)     # in cm^-1
	#bconst		   = 6.0;
	return bconst

def replace(string_old, string_new, file1, file2):
	'''
	This function replaces old string by new string
	'''
	f1             = open(file1, 'r')
	f2             = open(file2, 'w')
	for line in f1:
		f2.write(line.replace(string_old, string_new))
	f1.close()
	f2.close()

def beads(tau,beta):
	'''
	This function determins number of beads
	'''
	numbbeads1     =beta/tau+1
	numbbeads2     = int(round(numbbeads1,0))
	if (numbbeads2 % 2 == 0):
		numbbeads2 = numbbeads2 + 1
	return numbbeads2

def jobstring(file_name,value,thread):
	'''
	This function creats jobstring for #PBS script
	'''
	job_name       = "job_"+str(file_name)+str(value)
	walltime       = "200:00:00"
	processors     = "nodes=1:ppn="+str(thread)
	command_pimc_run = "./pimc"
	omp_thread     = str(thread)

	job_string     = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
####PBS -q medium
#PBS -l %s
#PBS -o %s.out
#PBS -e %s.err
export OMP_NUM_THREADS=%s
cd $PBS_O_WORKDIR
%s""" % (job_name, walltime, processors, job_name, job_name, omp_thread, command_pimc_run)
	print job_string
	return job_string

def outputstring1(numbbeads,tau,temperature):
	'''
	This function gives us the exact values of the agruments
	'''
	argu1          = "%7d"   % numbbeads
	argu2          = "%20.15f" % tau
	argu3          = "%7.5f" % temperature
	output ="numbbeads = "+argu1+", tau = "+argu2+", temperature = "+argu3+"\n"
	return output

def outputstring2(numbbeads,tau,mean_pot,mean_tot,mean_rot,error_pot,error_tot,error_rot):
	'''
	This function gives us the output 
	'''
	argu1          = "%5d"   % numbbeads
	argu2          = "%11.6f" % tau
	argu3          = "%10.5f" % mean_pot
	argu4          = "%10.5f" % mean_tot
	argu5          = "%10.5f" % mean_rot
	argu6          = "%10.5f" % error_pot
	argu7          = "%10.5f" % error_tot
	argu8          = "%10.5f" % error_rot
	output  = " "+argu1+" "+argu2+"   "+argu3+"     "+argu4+"     "+argu5
	output += "        "+argu6+"          "+argu7+"          "+argu8+"\n"
	return output


def formatting(status):
	'''
	This function gives us the output 
	'''
	if status == "analysis":
		output     = "#  Beads    Tau     Avg. Potential   Avg. Total   Avg. rotational  Error of Potential     Error of Total    Error of Rotational \n"
		output    += "#          (1/K)      Energy (K)     Energy (K)      Energy (K)        Energy (K)           Energy (K)          Energy (K) \n"
		output    += "#==============================================================================================================================\n"
		return output

def formatting1(status):
	'''
	This function gives us the output 
	'''
	if status == "analysis":
		output     = "#  Beads    Tau       Avg. CosTheta    Avg. Theta     Avg. Phi   Error of CosTheta    Error of Theta     Error of Phi  \n"
		output    += "#          (1/K)         (Radian)       (Radian)      (Radian)       (Radian)            (Radian)           (Radian)   \n"
	output    += "#==========================================================================================================================\n"
	return output


def modify_input(temperature,numbbeads,numbblocks,numbmolecules,distance):
	'''
	This function modifies parameters in qmc_run.input
	'''
	replace("temperature_input", str(temperature), "qmc_run.input", "qmc2.input")
	replace("numbbeads_input", str(numbbeads), "qmc2.input", "qmc3.input")
	replace("numbblocks_input", str(numbblocks), "qmc3.input", "qmc4.input")
	replace("numbmolecules_input", str(numbmolecules), "qmc4.input", "qmc5.input")
	replace("distance_input", str(distance), "qmc5.input", "qmc.input")
	call(["rm", "qmc_run.input", "qmc2.input", "qmc3.input", "qmc4.input", "qmc5.input"])


def rotmat(molecule,temperature,numbbeads):
	'''
	This function generates rotational density matrix - linden.dat
	'''
	#temperature1    = dropzeros(temperature)
	temperature1    = "%5.3f" % temperature
	numbbeads1		= numbbeads - 1
	command_linden_run = "../../../linear_prop/linden.x "+str(temperature)+" "+str(numbbeads1)+" "+str(bconstant())+" 1500 -1"
	print command_linden_run
	system(command_linden_run)
	file_rotdens    = molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	call(["mv", "linden.out", file_rotdens])
