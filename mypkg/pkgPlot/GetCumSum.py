#!/usr/bin/python
import os
import sys
from subprocess import call
import numpy as np
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from itertools import cycle
import itertools
import math
from math import *
 
import matplotlib.mlab as mlab

fig = plt.figure()

file1 = sys.argv[1]
print(file1)

col, x, y = loadtxt(file1, unpack=True, usecols=[0, 1, 2])
NMCumSum = np.cumsum(x)
DMCumSum = np.cumsum(y)
NormArr  = np.arange(1,int(len(NMCumSum)+1),dtype = 'float')
NMCumSum = np.divide(NMCumSum,NormArr)
DMCumSum = np.divide(DMCumSum,NormArr)
purity   = np.mean(NMCumSum)/np.mean(DMCumSum)

print("Purity is ")
print(purity)
print("")
print("Length of simulation is ")
print(len(NMCumSum))
print("")

plt.plot(col, NMCumSum, color = 'black', ls = '-', linewidth=1,  marker = 'None', markersize = 9, label = 'NM')
plt.plot(col, DMCumSum, color = 'red',   ls = '-', linewidth=1,  marker = 'None', markersize = 9, label = 'DM')
plt.ylim(0.3, 0.7)
plt.xlabel('Number of blocks')
plt.ylabel('Cumulative sum')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.legend(bbox_to_anchor=(0.35, 0.20), loc=2, borderaxespad=1., shadow=True )

#outfile = "hist-test.eps"
#plt.savefig(outfile, dpi = 200, format = 'eps')
plt.show()
