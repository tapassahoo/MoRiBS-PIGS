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
molecule = "H2"
numbbeads = 513
numbblocks = 4000
ntemp   = 10
tempmin = 0.0
tempmax = 2.5
dtemp   = (tempmax - tempmin)/ntemp

fldr = "run"+str(numbbeads)+"beads"+str(numbblocks)+"blocks"
src_path  = "/home/tapas/Moribs-pigs/MoRiBS-PIMC/examples/linear_pigs/"

# Loop over your jobs
for i in range(1, ntemp+1):
 
	temp = i*dtemp	

	folder_run = fldr+"%3.2fK" % temp
	call(["cp", "analyze.f", folder_run+"/results"])
	path_enter = src_path+folder_run+"/results"
	os.chdir(path_enter)

	file_pigs = "pigs%3.2fK" % temp
	replace("analyze_input", str(file_pigs)+".eng", "analyze.f", "analyze1.f")
	call(["mv", "analyze1.f", "analyze.f"])

	command_analyze_compile = "gfortran analyze.f -o analyze.x"
	command_analyze_run = "./analyze.x"
	system(command_analyze_compile)
	#system(command_analyze_run)

	os.chdir(src_path)

ntype_blocks = 5
for i in range(1, ntemp+1):
	temp = i*dtemp	
	file_pigs = "pigs%3.2fK" % temp
	for j in range(1, ntype_blocks+1):
		numbblocks = j*1000	
		fldr = "run"+str(numbbeads)+"beads"+str(numbblocks)+"blocks"

		folder_run = fldr+"%3.2fK" % temp
		fileread = src_path+folder_run+"/results/"+str(file_pigs)

		f1="blockavg%3.2fK"%temp
		filewrite = src_path+str(f1)

		print i, j, f1
		#fr = open(fileread, "r+")
		#fw = open("foo.txt", "a")
		#str = fr.read();
		#print "Read String is : ", str
		#fw.write( "%s") % str;
		#fw.close()
		#fo.close()

