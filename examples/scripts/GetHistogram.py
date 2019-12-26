#import math
#from math import *
import sys
import numpy as np
import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
#from itertools import cycle
#import itertools
#from scipy.optimize import curve_fit
#from subprocess import call
#from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib.ticker import ScalarFormatter
#import os
 
#import matplotlib.mlab as mlab

fig = plt.figure()


natom = int(sys.argv[1])
nbead = int(sys.argv[2])
naxis = int(sys.argv[3])

ncol1 = nbead + (natom-1)*3 # We have considered 0, M, N-1 beads
ncol  = naxis + (ncol1-1)*3 # Here we have taken cost and phi
ncol=ncol+1
print(ncol1)
print(ncol)

exit()
num_bins = 50
file1 = "/scratch/tapas/nonlinear-rotors/PIGS-RotDOFs-Rpt6.2Angstrom-beta0.1Kinv-Blocks20000-Passes100-System11-p-H2O-e0vsbeads91/results/output_instant.dof"
#file1 = "/scratch/tapas/nonlinear-rotors/PIGS-RotDOFs-Rpt8.2Angstrom-beta0.1Kinv-Blocks20000-Passes100-System11-p-H2O-e0vsbeads91/results/output_instant.dof"
x = np.loadtxt(file1, unpack=True, usecols=[ncol])
#y = np.loadtxt(file2, unpack=True, usecols=[ncol])
#x = np.fmod(x,2.0*pi)
#y = np.fabs(np.fmod(y,2.0*pi))
plt.hist(x, bins='auto', normed=1, facecolor='blue', alpha=0.7, label = "PIMC-CL")
#plt.hist(y, bins='auto', normed=1, facecolor='green', alpha=0.7, label = "PIMC")
plt.xlabel('bins')
plt.ylabel('Probability Distribution')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.legend(bbox_to_anchor=(0.35, 0.20), loc=2, borderaxespad=1., shadow=True )

#outfile = "hist-test.eps"
#plt.savefig(outfile, dpi = 200, format = 'eps')
plt.show()
#call(["open", outfile])
