#!/usr/bin/python
# Example PBS cluster job submission in Python
 
import time
from subprocess import call
from os import system
import os
import decimal

def dropzeros(number):
	mynum          = decimal.Decimal(number).normalize()
	# e.g 22000 --> Decimal('2.2E+4')
	return mynum.__trunc__() if not mynum % 1 else float(mynum)

 
def replace(string_old, string_new, file1, file2):
	f1 = open(file1, 'r')
	f2 = open(file2, 'w')
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

def exact_value(tau,numbbeads):
	'''
	This function determins exact value of beta
	'''
	beta_exact     = tau*(numbbeads - 1)
	return beta_exact

#initial parameters for qmc.input
numbblocks       = 10000
temperature	     = 10.0
beta             = 1.0/temperature               

numbbeads        = 50
dnumbbeads       = 3

src_path  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"
nrange           = numbbeads                                                                 #change
displacement     = dnumbbeads                                                                #change
file1_name       = "temp"+str(temperature)+"blocks"+str(numbblocks)+"e0vsbeads"              #change
file2_name       = ""                                                                        #change

fw = open("Energy-vs-tau-2-molecules-beta"+str(beta)+"-blocks"+str(numbblocks)+".txt", "a")  #change
value_min        = 0.0                                                                       #change

# Loop over your jobs
for i in range(1, nrange+1): 
 
	value = value_min + i*displacement
	numbbeads    = dropzeros(value)

	beta         = 1.0/temperature
	tau          = beta/(value-1)

	fldr         = file1_name+str(numbbeads)+file2_name
	folder_run = fldr
	call(["cp", "analyze.f", folder_run+"/results"])
	path_enter = src_path+folder_run+"/results"
	os.chdir(path_enter)

	file_pigs = "pigs"
	num_lines = sum(1 for line in open('pigs.eng'))
	print num_lines
	replace("nn_input", str(numbblocks), "analyze.f", "analyze1.f")
	replace("analyze_input", file_pigs+".eng", "analyze1.f", "analyze2.f")
	replace("param2", str(tau), "analyze2.f", "analyze3.f")
	call(["mv", "analyze3.f", "analyze.f"])
	call(["rm", "analyze1.f", "analyze2.f"])

	command_analyze_compile = "gfortran analyze.f -o analyze.x"
	command_analyze_run = "./analyze.x"
	system(command_analyze_compile)
	system(command_analyze_run)

	fr = open("average_pigs.eng", "r+")
	str1 = fr.read();

	os.chdir(src_path)

####################################################
#	beta         = value #1.0/temperature                                                  #change
#	temperature  = 1.0/beta                                                                #change

#	numbbeads    = beads(tau,beta)                                                         #change
#	beta_exact   = exact_value(tau,numbbeads)
#	temperature_exact  = 1.0/beta_exact
####################################################

	argu1        = "%7d" % numbbeads
	argu2        = "%7.5f" % temperature
	argu3        = "%7.5f" % beta
	fw.write(argu1+"   "+"   "+argu2+"   "+argu3+"  "+str1)
	fr.close()

fw.close()
