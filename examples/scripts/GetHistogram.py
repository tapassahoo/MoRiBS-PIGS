#!/usr/bin/python
#!/usr/bin/env python
import math
from math import *
import sys
import numpy as np
from numpy import *
import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
from itertools import cycle
import itertools
from scipy.optimize import curve_fit
from subprocess import call
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
import support
import os
 
import matplotlib.mlab as mlab

fig = plt.figure()


natom = int(sys.argv[1])
nbead = int(sys.argv[2])
naxis = int(sys.argv[3])

ncol1 = nbead + (natom-1)*3 # We have considered 0, M, N-1 beads
ncol  = naxis + (ncol1-1)*2 # Here we have taken cost and phi
print(ncol1)
print(ncol)

num_bins = 20
file1 = "/work/tapas/test-cluster-update/PIGS-cluster-update-correct-RotDOFs-Rpt10.05Angstrom-gFactor1.0-beta0.2Kinv-Blocks5000-Passes100-System8HF-e0vsbeads11/results/output_instant.dof"
file2 = "/work/tapas/test-cluster-update/PIGS-RotDOFs-Rpt10.05Angstrom-gFactor1.0-beta0.2Kinv-Blocks5000-Passes100-System8HF-e0vsbeads11/results/output_instant.dof"
x = loadtxt(file1, unpack=True, usecols=[ncol])
y = loadtxt(file2, unpack=True, usecols=[ncol])
y = np.fmod(y,2.0*pi)
plt.hist(x, bins='auto', normed=1, facecolor='blue', alpha=0.5, label = "PIGS-CL")
plt.hist(y, bins='auto', normed=1, facecolor='green', alpha=0.5, label = "PIGS")
plt.xlabel('bins')
plt.ylabel('Probability Distribution')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.legend(bbox_to_anchor=(0.35, 0.20), loc=2, borderaxespad=1., shadow=True )

outfile = "hist-test.eps"
plt.savefig(outfile, dpi = 200, format = 'eps')
plt.show()
#call(["open", outfile])
