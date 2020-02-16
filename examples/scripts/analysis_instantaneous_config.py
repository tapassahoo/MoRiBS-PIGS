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
from matplotlib.backends.backend_pdf import PdfPages
 
dir_read = sys.argv[1]
#sim_type = sys.argv[2]

file_read = dir_read + 'results/output.xyz'
string1 =dir_read[dir_read.find("PIGS"):-1]
numb_beads = int(string1[string1.find("beads")+5:])
numb_particle = int(string1[string1.find("System")+6:string1.find("-p-H2O")])
numb_blocks = int(string1[string1.find("Blocks")+6:string1.find("-Passes")])
sim_read = string1[string1.find("PIGS"):string1.find("PIGS")+4]

if (sim_read == "PIGS"):
	beads_pos = int((numb_beads-1)/2)

data_len = len(loadtxt(file_read, unpack=True, usecols=[0]))
save_data = np.zeros((numb_particle,3,data_len))
for i in range(numb_particle):
	ncol1 = beads_pos+i*numb_beads
	for j in range(3):	 
		ncol = j+ncol1*6
		ncol = ncol+1
		print(str(ncol)+'th column')

		save_data[i,j,:] = loadtxt(file_read, unpack=True, usecols=[ncol])


distance=np.zeros((3,data_len))
for j in range(3):	 
	distance[j,:]=save_data[1,j,:]-save_data[0,j,:]
vec_dist=np.sqrt(np.sum(np.square(distance),axis=0))
print('Average lattics spacing = '+str(np.mean(vec_dist))+" Angstrom")
 
#cost = loadtxt(file_read, unpack=True, usecols=[ncol+1])
#pyplot calling first
fig = plt.figure()
#plt.grid(True)

# Normalize density
num_bins = 20
#kwargs = dict(alpha=0.5, bins='auto', density=True, stacked=True)

#Histogram plot
plt.hist(vec_dist, bins='auto', normed=True, color='b', label='""')
#plt.hist(cost, bins='auto', normed=True, color='b', label='""')

plt.ylabel('Density')

plt.show()
