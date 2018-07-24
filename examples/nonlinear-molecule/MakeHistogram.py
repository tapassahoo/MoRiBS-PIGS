#!/usr/bin/python
#!/usr/bin/env python
import numpy as np
from numpy import *
import matplotlib
matplotlib.use('pdf')
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

num_bins = 50
x = loadtxt("output_instant.dof", unpack=True, usecols=[3])
plt.hist(x, num_bins, normed=1, facecolor='blue', alpha=0.5)
plt.xlabel('bins')
plt.ylabel('Probability')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)

outfile = "hist-test.pdf"
plt.savefig(outfile, dpi = 200, format = 'pdf')
#call(["okular", outfile])
call(["open", outfile])
