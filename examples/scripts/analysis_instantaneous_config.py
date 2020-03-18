import os
import sys
from subprocess import call
import math
from math import *
import numpy as np
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.axes
# matplotlib.use('eps')
import matplotlib.ticker as mtick
import matplotlib.ticker as mticker
#from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
from pylab import *
from scipy.optimize import curve_fit
import mypkg.pkgMoribs.support_without_parallel as support

rc('text', usetex=True)
size=24
params = {'legend.fontsize': size*0.5,
	'figure.figsize': (8,6),
	'axes.labelsize': size,
	'axes.titlesize': size,
	'xtick.labelsize': size*0.75,
	'ytick.labelsize': size*0.75,
	'axes.titlepad': size}
plt.rcParams.update(params)
matplotlib.rcParams.update({'font.size': size*0.75})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

dir_read = sys.argv[1]
numb_particle = int(sys.argv[2])
dofs_read1 = int(sys.argv[3])
#particle_index = int(sys.argv[2])
#axis_plot = sys.argv[3]
axis_index = {"cost":0, "phi":1, "chi":2}
dofs_read = {0:"TransAndRot", 1:"Rot", 2:"Trans"} 

file_read = dir_read + 'results/output.xyz'
string1 =dir_read[dir_read.find("PIGS"):-1]
numb_beads = int(string1[string1.find("beads")+5:])
numb_blocks = int(string1[string1.find("Blocks")+6:string1.find("-Passes")])
sim_read = string1[string1.find("PIGS"):string1.find("PIGS")+4]
extra_file_name=string1[string1.find("PIGS")+5:string1.find(dofs_read[dofs_read1])]
#dofs_read = string1[string1.find("moves")+6:string1.find("DOFs")]
#numb_particle = int(string1[string1.find("System")+6:string1.find("-p-H2O")])

ndofs=3
if (sim_read == "PIGS"):
	beads_pos = int((numb_beads-1)/2)
	#beads_pos = numb_beads-1

preskip = 5000
postskip = 0
data_len = len(genfromtxt(file_read, unpack=True, usecols=[0], skip_header=preskip, skip_footer=postskip))
workingNdim = int(math.log(data_len)/math.log(2))
trunc = int(data_len-2**workingNdim)
data_len=int(numb_blocks-(preskip+trunc))
save_data = np.zeros((numb_particle,ndofs,data_len))

for i in range(numb_particle):
	ncol1 = beads_pos+i*numb_beads
	for j in range(ndofs):	 
		ncol = j+ncol1*ndofs
		ncol = ncol+1
		print(str(ncol)+'th column')

		save_data[i,j,:] = genfromtxt(file_read, unpack=True, usecols=[ncol], skip_header=preskip+trunc, skip_footer=0)

if (dofs_read[dofs_read1] == "TransAndRot"):
	distance=np.zeros((ndofs,data_len),dtype=float)
	for j in range(ndofs):	 
		distance[j,:]=save_data[1,j,:]-save_data[0,j,:]
	vec_dist=np.sqrt(np.sum(np.square(distance),axis=0))
	print(vec_dist)
	rcom="{:3.6f}".format(np.mean(vec_dist))
	print('Average lattics spacing = '+str(rcom)+" Angstrom")
 
#vec_plot = save_data[particle_index,axis_index[axis_plot],:]
#pyplot calling first
'''
plt.xlabel('Number of blocks')
plt.ylabel('R ('+r'$\AA$'+')')
col_block = np.arange(numb_blocks)
data_plot = vec_dist

plt.plot(col_block, data_plot)
plt.show()
'''
if (extra_file_name == "qTIP4PF-"):
	label_str="PIGS/q-TIP4P/F"
	label_panel = "(b)"
elif (extra_file_name == "qspcfw-"):
	label_str="PIGS/q-SPC/Fw"
	label_panel = "(c)"
else:
	label_str="PIGS/TIP4P/2005"
	label_panel = "(a)"

num_bins = 20
if (dofs_read[dofs_read1] == "TransAndRot"):
	plt.xlabel(r'$\mathrm{bins \ of \ r \ (\AA)}$', labelpad=5)
	data_plot = vec_dist

plt.hist(data_plot, bins='auto', normed=1, color='green', alpha=0.75, edgecolor='black', label=label_str)
plt.ylabel(r'$\mathrm{\rho(r)}$',labelpad=10)

plt.xlim(2.4,3.4)
plt.ylim(0.0,4.5)
xmin,xmax=plt.xlim()
ymin,ymax=plt.ylim()
plt.text(xmin+0.7*(xmax-xmin),ymin+0.8*(ymax-ymin),r'$\mathrm{\langle r \rangle = \ }$'+str(rcom)+r'$\mathrm{\AA}$')
plt.text(xmin+(xmax-xmin)*0.01,ymax-(ymax-ymin)*0.08,label_panel)

plt.xlim(2.3905,3.409)
plt.ylim(0.0,4.51)
plt.subplots_adjust(top=0.99,bottom=0.13,left=0.1,right=0.99,hspace=0.0,wspace=0.0)
plt.legend(numpoints=1,loc=('upper right'))

FilePlotDensity="/home/tapas/ResultsOfPIGS"+dir_read[28:-1]+"-histogram-r-preskip"+str(preskip)+"-postskip"+str(postskip)+".eps"
print(FilePlotDensity)
plt.savefig(FilePlotDensity, dpi=50, format='eps')
plt.show()

