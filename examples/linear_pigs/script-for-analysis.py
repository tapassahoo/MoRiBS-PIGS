#!/usr/bin/python
# Example PBS cluster job submission in Python
 
import time
from subprocess import call
from os import system
import os
 
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
numbblocks       = 8000
temperature	     = 10.0
tau 		     = 0.001

ntau    	     = 10
dtau    	     = 0.001

rmin             = 3.0
rmax             = 10.0
nr               = 70
dr               = (rmax-rmin)/nr

dbeta            = 0.004
nbeta            = 50

src_path  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"
nrange           = nbeta                                                                    #change
displacement     = dbeta                                                                    #change
file1_name       = "tau"+str(tau)+"blocks"+str(numbblocks)+"e0vsbeta"                       #change
file2_name       = "K-1"                                                                    #change
fw = open("Energy-vs-beta-2-molecules-tau"+str(tau)+"blocks"+str(numbblocks)+".txt", "a")   #change
value_min        = 0.0                                                                      #change

# Loop over your jobs
for i in range(1, nrange+1): 
 
	value = value_min + i*displacement                                
	r_dis        = "%5.3f" % value 

	fldr         = file1_name+r_dis+file2_name
	folder_run = fldr
	call(["cp", "analyze.f", folder_run+"/results"])
	path_enter = src_path+folder_run+"/results"
	os.chdir(path_enter)

	file_pigs = "pigs"
	num_lines = sum(1 for line in open('pigs.eng'))
	print num_lines
	replace("nn_input", str(numbblocks), "analyze.f", "analyze1.f")
	replace("analyze_input", file_pigs+".eng", "analyze1.f", "analyze2.f")
	replace("param2", str(r_dis), "analyze2.f", "analyze3.f")
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
	beta         = value #1.0/temperature                                                  #change
	temperature  = 1.0/beta                                                                #change

	numbbeads    = beads(tau,beta)                                                         #change
	beta_exact   = exact_value(tau,numbbeads)
	temperature_exact  = 1.0/beta_exact
####################################################

	argu1        = "%7d" % numbbeads
	argu2        = "%7.5f" % temperature_exact
	argu3        = "%7.5f" % beta
	fw.write(argu1+"   "+"   "+argu2+"   "+argu3+"  "+str1)
	fr.close()

fw.close()
