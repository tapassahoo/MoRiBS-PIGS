#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from os import system
from sys import argv

def plothistogram(srcfile, numbmolecules, molecule, dipolemoment, beta, tau,Rpt):
	col_costheta1, col_theta1 = loadtxt(srcfile, unpack=True, usecols=[5,6])
	col_costheta2, col_theta2 = loadtxt(srcfile, unpack=True, usecols=[7,8])

	hist1, bins1 = np.histogram(col_costheta1,nbins,density=True)
	hist2, bins2 = np.histogram(col_costheta2,nbins,density=True)
	print len(hist1), len(bins1)
	print np.sum(hist1)
	print np.sum(hist1*np.diff(bins1))

	plt.hist( col_costheta1, bins='auto', normed = 1)  # plt.hist passes it's arguments to np.histogram
	plt.hist( col_costheta2, bins='auto', normed = 1)  # plt.hist passes it's arguments to np.histogram
	plt.grid(True)


	plt.suptitle('Parameters: System '+str(numbmolecules)+" "+str(molecule)+", "+r'$\mu$ = '+str(dipolemoment)+' Debye, '+r'$\beta$ = '+str(beta)+' '+r'$K^{-1}$, '+r'$\tau$ = '+str(tau)+' '+r'$K^{-1}$'+' R = '+str(Rpt)+r'$\AA$', fontsize =10)
	plt.xlabel(r'$\cos(\theta)$', fontsize = 20)
	plt.ylabel('Density', fontsize = 20)
	plt.xlim(-1.01,1.01)
	plt.grid(True)
	plt.show()

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
#molecule            = "HF-C60"                                                    #change param1
molecule            = "HF"                                                         #change param1
#molecule            = "H2"                                                        #change param1
molecule_rot        = "HF"                                                         #change param2
numbblocks	        = 20000                                                        #change param3
numbmolecules       = 2                                                            #change param4
numbbeads           = 129                                                          #change param5
beta     	        = 0.128                                                        #change param6
dipolemoment        = 0.45                                                         #change param7
dRpt                = 0.5                                                          #change param7

nrange              = 20		  						                           #change param11

temperature         = 1.0/beta   
tau                 = beta/(numbbeads-1)

file1_name          = "beta"+str(beta)+"Kinv-tau"+str(tau)+"Kinv-DipoleMoment"+str(dipolemoment)+"Debye-Blocks"+str(numbblocks)
file1_name         += "-System"+str(numbmolecules)+str(molecule)+"-Beads"+str(numbbeads)+"-e0vsRpt" 

file2_name          = "Angstrom"                                                   #change param10
argument2           = "Rpt"                                                        #change param11
value_min           = 0.5                                                            #change param12
var                 = "Rpt"                                                        #change param13

src_path            = os.getcwd()
dest_path           = "/work/tapas/linear_rotors/"                                 #change param13
run_file            = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/pimc"                   #change param14

trunc               = 20000
nbins = 100          #change

# Loop over jobs
for i in range(nrange):                                                  #change param13
 
	value        = i*dRpt + value_min

	Rpt          = '{:2.1f}'.format(value)

	folder_run   = file1_name+str(Rpt)+file2_name
	dest_dir     = dest_path + folder_run + "/results"
	src_file     = dest_dir+"/pigs.dof"

	plothistogram(src_file, numbmolecules, molecule, dipolemoment, beta, tau, Rpt)
