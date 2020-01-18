import decimal
import glob
import math
import os
import time
from subprocess import call

import numpy as np


def dropzeros(number):
	mynum = decimal.Decimal(number).normalize()
	# e.g 22000 --> Decimal('2.2E+4')
	return mynum.__trunc__() if not mynum % 1 else float(mynum)


def GetDirNameSubmission(rotor, temperature, numbbeads, iodevn, jmax):
	if iodevn == 0:
		spinstate = "-para-"
	elif iodevn == 1:
		spinstate = "-ortho-"
	elif iodevn == -1:
		spinstate = "-spin-less-"

	file_name = (
		"rot-dens-matrix-of"
		+ spinstate
		+ rotor
		+ "-Temp"
		+ str(temperature)
		+ "Kelvin-Beads"
		+ str(numbbeads)
	)

	return file_name


def MakeExecutable(src_dir, source_dir_exe):
	"""It copies Makefile-copy as Makefile and run make clean && make to generate executable file asymrho.x.
	
	Arguments:
		src_dir {str} -- The path of this script.
		source_dir_exe {str} -- The path of the asymrho.f code.
	"""
	print("")
	print("#-----------------------Compilation started--------------------#")
	print("")
	os.chdir(source_dir_exe)
	call(["make", "clean"])
	call(["cp", "Makefile-Copy", "Makefile"])
	call(["make"])
	print("")
	print("#-----------------------Compilation ended here--------------------#")
	print("")
	os.chdir(src_dir)


def GetABCconsts(rotor):
	if rotor == "H2O":
		aconst = 27.877
		bconst = 14.512
		cconst = 9.285
	abc = {"aconst": aconst, "bconst": bconst, "cconst": cconst}
	return abc


def Submission(dir_job, script_dir, execution_file, numbbeads, temperature, molecule, dir_input, dir_output, NameOfPartition, user_name, out_dir, dir_name, iodevn, jmax, NameOfServer):
	for theta in range(181):
		folder_run = "theta" + str(theta)
		folder_run_path = dir_job + "/" + folder_run
		fname = dir_input + "/job-for-theta" + str(theta)
		logfile = "t" + str(theta) + "P" + str(numbbeads) + "T" + str(temperature)
		print(theta)

		fwrite = open(fname, "w")
		fwrite.write( jobstring(logfile,folder_run_path,dir_input,dir_output,temperature,numbbeads,iodevn,theta,molecule,jmax,folder_run, NameOfServer))

		fwrite.close()

		os.chdir(dir_input)
		if NameOfPartition == "tapas":
			call(["sbatch", "-p", "tapas", fname])
		else:
			call(["sbatch", fname])

		os.chdir(script_dir)


def jobstring(logfile,folder_run_path,dir_input,dir_output,temperature,numbbeads,iodevn,theta,molecule,jmax,folder_run,NameOfServer):
	"""
	This function creats jobstring for #SBATCH script
	"""
	logpath = dir_input + "/" + logfile
	exe_file = dir_input + "/asymrho.x"

	aconst = GetABCconsts(molecule)["aconst"]
	bconst = GetABCconsts(molecule)["bconst"]
	cconst = GetABCconsts(molecule)["cconst"]

	command_execution = ("./asymrho.x  " + str(temperature) + " " + str(numbbeads) + " " + str(iodevn) + " " + str(theta) + " " + str(theta) + " " + str(aconst) + " " + str(bconst) + " " + str(cconst) + " " + str(jmax))

	mvfiles = "mv " + dir_output + "/" + folder_run + "/* " + dir_input
	rmfolder = "rm -rf " + dir_output + "/" + folder_run
	if NameOfServer == "graham":
		CommandForMove = " "
		account = "#SBATCH --account=rrg-pnroy"
	else:
		account = ""

	job_string = """#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --output=%s.out
#SBATCH --time=1-03:00
%s
#SBATCH --mem-per-cpu=512mb
#SBATCH --cpus-per-task=1
export OMP_NUM_THREADS=1
rm -rf %s
mkdir -p %s
cd %s
cp %s %s
%s
%s
%s
""" % (logfile,logpath,account,folder_run_path,folder_run_path,folder_run_path,exe_file,folder_run_path,command_execution,mvfiles,rmfolder)
	return job_string


def GetPackRotDens(src_dir_exe, dir_output, script_dir, rotor, temperature, numbbeads):
	cmd_cp = "cp " + src_dir_exe + "*.x " + dir_output
	os.system(cmd_cp)

	os.chdir(dir_output)

	file0 = dir_output + "/rho.den000"
	if os.path.exists(file0) != True:
		print("Stop. There is no rho.dens*.")
		exit()

	cmd_exe = "./compile.x"
	os.system(cmd_exe)
	RemoveAllFiles(dir_output, rotor, temperature, numbbeads)

	os.chdir(script_dir)


def RemoveAllFiles(dir_output, rotor, temperature, numbbeads):
	filelist = glob.glob(dir_output + "/*")

	file0 = dir_output + "/rho.den"
	filelist.remove(file0)
	file1 = dir_output + "/rho.den_rho"
	filelist.remove(file1)
	file2 = dir_output + "/rho.den_eng"
	filelist.remove(file2)
	file3 = dir_output + "/rho.den_esq"
	filelist.remove(file3)
	for file in filelist:
		os.remove(file)

	subscript = (
		dir_output + "/" + rotor + "_T" + str(dropzeros(temperature)) + "t" + str(int(numbbeads))
	)
	os.rename(file1, subscript + ".rho")
	os.rename(file2, subscript + ".eng")
	os.rename(file3, subscript + ".esq")


if __name__ == "__main__":
	print("Tapas, I am here!")
