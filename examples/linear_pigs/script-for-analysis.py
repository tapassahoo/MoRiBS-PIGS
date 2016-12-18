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

#initial parameters for qmc.input
ntau    	     = 20
dtau    	     = 0.01

rmin             = 0.0
rmax             = 5.0
nr               = 20
dr               = (rmax-rmin)/nr

dbeta            = 0.02
nbeta            = 20

src_path  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"
nrange           = ntau                           #change
displacement     = dtau                           #change
file1_name       = "e0vstau"                      #change
file2_name       = "K-1"                          #change
fw = open("Energy-vs-tau.txt", "a")               #change

# Loop over your jobs
for i in range(1, nrange+1): 
 
	value = i*displacement          
	r_dis        = "%5.3f" % value 

	fldr         = file1_name+r_dis+file2_name
	folder_run = fldr
	call(["cp", "analyze.f", folder_run+"/results"])
	path_enter = src_path+folder_run+"/results"
	os.chdir(path_enter)

	file_pigs = "pigs"
	num_lines = sum(1 for line in open('pigs.eng'))
	print num_lines
	replace("nn_input", str(num_lines), "analyze.f", "analyze1.f")
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

	fw.write(str1)
	fr.close()

fw.close()
