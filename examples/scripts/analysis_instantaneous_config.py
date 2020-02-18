import os
import sys
from subprocess import call
import math
from math import *
import numpy as np
from numpy import *
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
 
dir_read = sys.argv[1]
particle_index = int(sys.argv[2])
axis_plot = sys.argv[3]
axis_index = {"cost":0, "phi":1, "chi":2}

file_read = dir_read + 'results/output.xyz'
string1 =dir_read[dir_read.find("PIGS"):-1]
numb_beads = int(string1[string1.find("beads")+5:])
numb_particle = int(string1[string1.find("System")+6:string1.find("-p-H2O")])
numb_blocks = int(string1[string1.find("Blocks")+6:string1.find("-Passes")])
sim_read = string1[string1.find("PIGS"):string1.find("PIGS")+4]
dofs_read = string1[string1.find("PIGS")+5:string1.find("DOFs")]

if (sim_read == "PIGS"):
	beads_pos = int((numb_beads-1)/2)

if (dofs_read == "TransAndRot"):
	ndofs=6
else:
	ndofs=3

data_len = len(loadtxt(file_read, unpack=True, usecols=[0]))
workingNdim = int(math.log(data_len)/math.log(2))
trunc = int(data_len-2**workingNdim)
save_data = np.zeros((numb_particle,3,data_len))
for i in range(numb_particle):
	ncol1 = beads_pos+i*numb_beads
	for j in range(3):	 
		ncol = j+ncol1*ndofs
		ncol = ncol+1
		print(str(ncol)+'th column')

		save_data[i,j,:] = loadtxt(file_read, unpack=True, usecols=[ncol])

if (dofs_read == "TransAndRot"):
	distance=np.zeros((3,data_len))
	for j in range(3):	 
		distance[j,:]=save_data[1,j,:]-save_data[0,j,:]
	vec_dist=np.sqrt(np.sum(np.square(distance),axis=0))
	print('Average lattics spacing = '+str(np.mean(vec_dist))+" Angstrom")
 
#vec_plot = save_data[particle_index,axis_index[axis_plot],:]
#pyplot calling first
fig = plt.figure()
num_bins = 20
if (dofs_read == "TransAndRot"):
	plt.xlabel('Bins of lattice spacing ('+r'$\AA$'+')')
	data_plot = vec_dist[trunc:]

print(save_data[particle_index,axis_index[axis_plot],:])
if (dofs_read == "Rot"):
	data_plot = save_data[particle_index,axis_index[axis_plot],:]

	if (axis_index[axis_plot] == 0):
		plt.xlabel('Bins of '+r'$\cos(\theta)$')
		plt.xlim(-1.0,1.0)
	elif (axis_index[axis_plot] == 1):
		plt.xlabel('Bins of '+r'$\phi$')
		plt.xlim(0.0,2.0*math.pi)
	else:
		plt.xlabel('Bins of '+r'$\chi$')
		plt.xlim(0.0,2.0*math.pi)

plt.hist(data_plot, bins='auto', normed=1, color='green', alpha=0.75, edgecolor='black', label='""')
plt.ylabel('Density')
plt.show()
