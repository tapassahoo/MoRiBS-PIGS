#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *

def compile_rotmat():
	path_enter_linden = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/linear_prop/"
	os.chdir(path_enter_linden)
	call(["make", "clean"])
	call(["make"])
	path_exit_linden  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"
	os.chdir(path_exit_linden)

def jackknife(mean,data):
	ai            = [((np.sum(data) - data[j])/(len(data) - 1.0)) for j in range(len(data))]
	deviation     = ai - mean
	devsq         = np.multiply(deviation,deviation)
	error         = sqrt(np.sum(devsq)*(len(data)-1.0)/len(data))
	return error

def levels(number):
	for j in range (10):
		jj=pow(2,j)

		if jj <= (number-1):
			level = j
		else:
			break
	return level

def dropzeros(number):
    mynum          = decimal.Decimal(number).normalize()
    # e.g 22000 --> Decimal('2.2E+4')
    return mynum.__trunc__() if not mynum % 1 else float(mynum)

def bconstant(molecule_rot):
	'''
	This function calculates rotational Bconstant for linear rotor
	'''
	autocminverse  = 2.1947463137e+5
	energyj0       = -36117.5942855
	energyj1       = -35999.1009407
	bconst         = 0.5*(energyj1-energyj0)     # in cm^-1
	if (molecule_rot == "HF"):
		bconst	   = 20.9561                     # in cm^-1  and it is  taken from http://webbook.nist.gov/cgi/inchi?ID=C7664393&Mask=1000#Diatomic
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
##PBS -q medium
#PBS -l %s
#PBS -o %s.out
#PBS -e %s.err
export OMP_NUM_THREADS=%s
cd $PBS_O_WORKDIR
%s""" % (job_name, walltime, processors, job_name, job_name, omp_thread, command_pimc_run)
	print job_string
	return job_string

def inputstr(numbbeads,tau,temperature):
	'''
	This function gives us the exact values of the agruments
	'''
	argu1          = "%7d"   % numbbeads
	argu2          = "%20.15f" % tau
	argu3          = "%7.5f" % temperature
	output ="numbbeads = "+argu1+", tau = "+argu2+", temperature = "+argu3+"\n"
	return output

def outputstr_energy(numbbeads,tau,dest_dir,trunc):
	'''
	This function gives us the output 
	'''
	col_block, col_pot, col_tot, col_rot = genfromtxt(dest_dir+"/results/pigs.eng",unpack=True, usecols=[0,1,2,3], skip_header=5000, max_rows=trunc)
	print len(col_tot)
	
	mean_pot      = np.mean(col_pot)
	mean_tot      = np.mean(col_tot)
	mean_rot      = np.mean(col_rot)

	error_pot     = jackknife(mean_pot,col_pot)
	error_tot     = jackknife(mean_tot,col_tot)
	error_rot     = jackknife(mean_rot,col_rot)
	#print i, len(col_block)

	output  = '{:10d}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}'.format(numbbeads, tau, mean_pot, mean_tot, mean_rot, error_pot, error_tot, error_rot)
	output  += "\n"
	return output

def outputstr_angle(numbbeads,tau,dest_dir,trunc):
	'''
	This function gives us the output 
	'''
	col_block, col_costheta, col_theta, col_costheta1, col_theta1 = genfromtxt(dest_dir+"/results/pigs.dof",unpack=True, usecols=[0,1,2,3,4], skip_header=5000,max_rows=trunc)

	mean_costheta  = np.mean(col_costheta)
	mean_theta     = np.mean(col_theta)
	mean_costheta1 = np.mean(col_costheta1)
	mean_theta1    = np.mean(col_theta1)

	error_costheta = jackknife(mean_costheta,col_costheta)
	error_theta    = jackknife(mean_theta,col_theta)
	error_costheta1= jackknife(mean_costheta1,col_costheta1)
	error_theta1   = jackknife(mean_theta1,col_theta1)

	output  = '{:10d}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}'.format(numbbeads, tau, mean_costheta, mean_theta, error_costheta, error_theta, mean_costheta1, mean_theta1, error_costheta1, error_theta1)
	output  += "\n"
	return output

def outputstr_angle1(numbbeads,tau,dest_dir,trunc):
	'''
	This function gives us the output 
	'''
	col_block, col_costheta, col_theta, col_costheta1, col_theta1 = genfromtxt(dest_dir+"/results/pigs.dof",unpack=True, usecols=[0,5,6,7,8], skip_header=5000, max_rows=trunc)

	mean_costheta  = np.mean(col_costheta)
	mean_theta     = np.mean(col_theta)
	mean_costheta1 = np.mean(col_costheta1)
	mean_theta1    = np.mean(col_theta1)

	error_costheta = jackknife(mean_costheta,col_costheta)
	error_theta    = jackknife(mean_theta,col_theta)
	error_costheta1= jackknife(mean_costheta1,col_costheta1)
	error_theta1   = jackknife(mean_theta1,col_theta1)

	output  = '{:10d}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}'.format(numbbeads, tau, mean_costheta, mean_theta, error_costheta, error_theta, mean_costheta1, mean_theta1, error_costheta1, error_theta1)
	output  += "\n"
	return output

def fmt_energy(status,variable):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('Beads', variable, 'Avg. Potential', 'Avg. Total', 'Avg. rotational', 'Error of Potential', 'Error of Total', 'Error of Rotational')
		output    +="\n"
		output    +="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('', (str(unit)), 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)', 'Energy (K)')
		output    +="\n"
		output    += '{:=<155}'.format('#')
		output    +="\n"
		return output

def fmt_angle(status,variable):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('Beads', variable, 'Avg. CosTheta', 'Avg. Theta', 'Error of CosTheta', 'Error of Theta', 'Avg. CosTheta1', 'Avg. Theta1', 'Error of CosTheta1', 'Error of Theta1')
		output    +="\n"
		output    +="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('', (str(unit)), '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)')
		output    +="\n"
		output    += '{:=<195}'.format('#')
		output    +="\n"

	return output

def fmt_angle1(status,variable):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('Beads', variable, 'Avg. CosTheta', 'Avg. Theta', 'Error of CosTheta', 'Error of Theta', 'Avg. CosTheta1', 'Avg. Theta1', 'Error of CosTheta1', 'Error of Theta1')
		output    +="\n"
		output    +="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('', (str(unit)), '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)', '(Radian)')
		output    +="\n"
		output    += '{:=<195}'.format('#')
		output    +="\n"

	return output



def modify_input(temperature,numbbeads,numbblocks,molecule_rot,numbmolecules,distance,level,step,dipolemoment):
	'''
	This function modifies parameters in qmc_run.input
	'''
	replace("temperature_input", str(temperature), "qmc_run.input", "qmc2.input")
	replace("numbbeads_input", str(numbbeads), "qmc2.input", "qmc3.input")
	replace("numbblocks_input", str(numbblocks), "qmc3.input", "qmc4.input")
	replace("numbmolecules_input", str(numbmolecules), "qmc4.input", "qmc5.input")
	replace("distance_input", str(distance), "qmc5.input", "qmc6.input")
	replace("molecule_input", str(molecule_rot), "qmc6.input", "qmc7.input")
	replace("level_input", str(level), "qmc7.input", "qmc8.input")
	replace("dstep_input", str(step), "qmc8.input", "qmc9.input")
	replace("dipolemoment_input", str(dipolemoment), "qmc9.input", "qmc.input")
	call(["rm", "qmc2.input", "qmc3.input", "qmc4.input", "qmc5.input", "qmc6.input", "qmc7.input", "qmc8.input", "qmc9.input"])


def rotmat(molecule,temperature,numbbeads):
	'''
	This function generates rotational density matrix - linden.dat
	'''
	#temperature1    = dropzeros(temperature)
	temperature1    = "%5.3f" % temperature
	numbbeads1		= numbbeads - 1
	command_linden_run = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/linear_prop/linden.x "+str(temperature)+" "+str(numbbeads1)+" "+str(bconstant(molecule))+" 15000 -1"
	print command_linden_run
	system(command_linden_run)
	file_rotdens    = molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	call(["mv", "linden.out", file_rotdens])

def rotmat1(molecule,temperature,numbbeads):
	'''
	This function generates rotational density matrix - linden.dat
	'''
	#temperature1    = dropzeros(temperature)
	temperature1    = "%5.3f" % temperature
	numbbeads1		= numbbeads - 1
	command_linden_run = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/linear_prop/linden.x "+str(temperature)+" "+str(numbbeads1)+" "+str(bconstant(molecule))+" 3000 -1"
	print command_linden_run
	system(command_linden_run)

def jobstring_scratch(file_name, value, thread, run_dir, molecule, temperature, numbbeads, final_dir):
	'''
	This function creats jobstring for #PBS script
	'''
	job_name       = "job_"+str(file_name)+str(value)
	walltime       = "200:00:00"
	processors     = "nodes=1:ppn="+str(thread)
	command_pimc_run = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"
	omp_thread     = str(thread)
	output_dir     = run_dir+"/results"
	temperature1    = "%5.3f" % temperature
	file_rotdens    = molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
	logpath        = final_dir+"/"

	input_file    = "qmc_"+file_name+str(value)+".input"
	print input_file

	job_string     = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
##PBS -q medium
#PBS -l %s
#PBS -o %s%s.out
#PBS -e %s%s.err
export OMP_NUM_THREADS=%s
rm -rf %s
mkdir -p %s
cd $PBS_O_WORKDIR
cp IhRCOMC60.xyz %s
mv %s %s
mv %s %s 
cd %s
cp %s qmc.input
cp /home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc %s
./pimc
mv %s /work/tapas/linear_rotors
""" % (job_name, walltime, processors, logpath, job_name, logpath, job_name, omp_thread, run_dir, output_dir, run_dir, input_file, run_dir, file_rotdens, run_dir, run_dir, input_file, run_dir, run_dir)
	print job_string
	return job_string

def outputstr_entropy(numbbeads,tau,dest_dir,trunc):
	'''
	This function gives us the output 
	'''
	col_block, col_NM, col_DM, col_EN = genfromtxt(dest_dir+"/results/pigs.rden",unpack=True, usecols=[0,1,2,3], max_rows=trunc)
	print len(col_NM)
	
	mean_NM      = np.mean(col_NM)
	mean_DM      = np.mean(col_DM)
	mean_EN      = np.mean(col_EN)

	error_NM     = jackknife(mean_NM,col_NM)
	error_DM     = jackknife(mean_DM,col_DM)
	error_EN     = jackknife(mean_EN,col_EN)
	#print i, len(col_block)

	output  = '{:10d}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}{:20.5f}'.format(numbbeads, tau, mean_NM, mean_DM, mean_EN, error_NM, error_DM, error_EN)
	output  += "\n"
	return output

def fmt_entropy(status,variable):
	'''
	This function gives us the output 
	'''
	if variable == "Rpt":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output     ="#"
		output    += '{:^15}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}{:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Avg. Entropy', 'Error of Nm', 'Error of Dm', 'Error of Entropy')
		output    +="\n"
		output    += '{:=<155}'.format('#')
		output    +="\n"
		return output

