#!/usr/bin/python

import os
from subprocess import call

PATH = '/work/tapas/linear_rotors/ENT-RotDOFs-Rpt10.05Angstrom-DipoleMoment6.0Debye-beta0.2Kinv-Blocks200000-Passes200-System2HF-ParticleA1-e0vsbeads-SWAPTOUNSWAP101/results'


os.chdir(PATH)

for f in os.listdir(PATH):
	#print(os.path.splitext(f))

	if f.endswith("_sub1"):
		#print(os.path.join(f))
		file_final = f[:-5]
		file_sub1  = f[:-5]+"_sub1"
		print(file_final)
		print(file_sub1)

		filenames = [file_sub1, file_final]
		file_result = file_sub1+"_result"
		with open(file_result, 'w') as outfile:
			for fname in filenames:
				with open(fname) as infile:
					for line in infile:
						outfile.write(line)
		call(["mv", file_result, file_final])
		call(["rm", file_sub1])
