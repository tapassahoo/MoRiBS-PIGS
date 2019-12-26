import math
from math import *
import sys
import numpy as np
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
import os
 
#pyplot calling first
fig = plt.figure()
plt.grid(True)

#Input specifications
natom = int(sys.argv[1])
#nbeads = int(sys.argv[2])
nbeads = 1
naxis = int(sys.argv[2])
#NumBeads = int(sys.argv[3])
#DipoleMoment = float(sys.argv[4])

ncol1 = nbeads + (natom-1)*3 # We have considered 0, M, N-1 beads
ncol = naxis + ncol1*3 # We have considered costheta, phi, chi
ncol = ncol+1
print(str(ncol+1)+'th column')
beadposition={0:'first-bead', 1:'middle-bead', 2:'last-bead'}

#Import the data files here
#file1 = "/work/tapas/linear-rotors-PIGS/PIGS-RotDOFs-Rpt10.05Angstrom-DipoleMoment"+str(DipoleMoment)+"Debye-beta0.2Kinv-Blocks20000-Passes200-System2HF-e0vsbeads"+str(NumBeads)+"/results/output_instant.dof"
file1 = "/scratch/tapas/nonlinear-rotors/PIGS-RotDOFs-Rpt8.0Angstrom-beta0.1Kinv-Blocks20000-Passes100-System11-p-H2O-e0vsbeads51/results/output_instant.dof"

#Data processiong by Numpy 
x1 = loadtxt(file1, unpack=True, usecols=[ncol])

if (naxis != 0): 
	x1 = np.fabs(np.fmod(x1,2.0*pi))

# Normalize density
num_bins = 20
kwargs = dict(alpha=0.5, bins='auto', density=True, stacked=True)

# Lebels for plotting
#label1 = str(DipoleMoment)+' Debye'
label1 = ''
#label1 += ';  P = '+str(NumBeads)
print(label1)

#Histogram plot
plt.hist(x1, **kwargs, color='b', label='""')

#X and Y labels
if (naxis != 0): 
	plt.xlabel('bins of '+r'$\phi$')
	fname='Phi'
else:
	plt.xlabel('bins of '+r'$\cos(\theta)$')
	fname='CosTheta'

plt.ylabel('Density')

# Tweak spacing to prevent clipping of ylabel
#plt.gca().set(title=r'$\mu = $'+label1)
plt.xlim(-1,1)
plt.legend();
plt.subplots_adjust(left=0.15)
plt.legend(bbox_to_anchor=(0.25, 2.20), loc=2, borderaxespad=1., shadow=True )

#Saving the plot in a file
#outfile = "/home/tapas/ResultsOfPIGS/PIGS-RotDOFs-Rpt10.05Angstrom-DipoleMoment"+str(DipoleMoment)+"Debye-beta0.2Kinv-Blocks20000-Passes200-System2HF-e0vsbeads"+str(NumBeads)+"-histgramOf"+fname+"-for-"+str(natom)+"th-rotor-and-"+beadposition[nbeads]+".eps"
#print(outfile)
#plt.savefig(outfile, dpi = 100, format = 'eps')
plt.show()
