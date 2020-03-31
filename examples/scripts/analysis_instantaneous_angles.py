import os
import sys
from subprocess import call
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes
#matplotlib.use('eps')
from pylab import *
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

#Read some parameters from the command-line
dir_read = sys.argv[1]
mc_read = sys.argv[2]
numb_particle = int(sys.argv[3])
dofs_read1 = int(sys.argv[4])
axis_read= sys.argv[5]
numb_beads= int(sys.argv[6])
particle_index = int(sys.argv[7])

#
#particle_index = int(sys.argv[2])
axis_index = {"cost":0, "phi":1, "chi":2}
dofs_read = {0:"TransAndRot", 1:"Rot", 2:"Trans"} 

#
file_read = dir_read + 'results/output.xyz'
string1 =dir_read[dir_read.find(mc_read):-1]
numb_blocks = int(string1[string1.find("Blocks")+6:string1.find("-Passes")])
sim_read = string1[string1.find(mc_read):string1.find(mc_read)+4]
extra_file_name=string1[string1.find(mc_read)+5:string1.find(dofs_read[dofs_read1])]
gfact=string1[string1.find("gFactor")+7:string1.find("-beta")]
#numb_beads = int(string1[string1.find("beads")+5:])
#dofs_read = string1[string1.find("moves")+6:string1.find("DOFs")]
#numb_particle = int(string1[string1.find("System")+6:string1.find("-p-H2O")])
print(sim_read)

first_fragment = dir_read[:dir_read.find("gFactor")]
last_fragment = dir_read[dir_read.find("-beta"):]

ndofs=3
if (sim_read != "PIMC"):
	beads_pos = int((numb_beads-1)/2)
	#beads_pos = numb_beads-1

preskip = 0
postskip = 0

print("#Look at the below line")
print(file_read)
data_len = len(genfromtxt(file_read, unpack=True, usecols=[0], skip_header=preskip, skip_footer=postskip))
workingNdim = int(math.log(data_len)/math.log(2))
trunc = int(data_len-2**workingNdim)
data_len=int(numb_blocks-(preskip+trunc))
save_data = np.zeros((numb_particle,data_len))

label_panel = "(a)"
if (numb_particle == 4):
	gFactList=[2.0, 4.0, 8.0, 16.0]
elif (numb_particle == 16):
	gFactList=[0.5, 1.0, 1.5, 2.0]
colorList=["yellow", "green", "red", "magenta"]
j=axis_index[axis_read]
ig=0
for g in gFactList:
	file_read=first_fragment+'gFactor'+str(g)+last_fragment + 'results/output.xyz'

	for i in range(numb_particle):
		ncol1 = beads_pos+i*numb_beads
		ncol = j+ncol1*ndofs
		ncol = ncol+1
		print(str(ncol)+'th column')
		save_data[i,:] = genfromtxt(file_read, unpack=True, usecols=[ncol], skip_header=preskip+trunc, skip_footer=0)

	label_str=r'$g$='+str(g)	

	#for i in range(numb_particle):
	vec_plot = np.reshape(save_data, numb_particle*data_len)
	data_plot = vec_plot
	print(len(data_plot))
	plt.hist(data_plot, bins=50, density=True, alpha=0.5, stacked=True, color=colorList[ig], edgecolor='black', label=label_str)
	ig=ig+1

#
if (mc_read == "ENT"):
	numb_label=int(numb_particle/2)
plt.xlabel(r'$\mathrm{bins \ of \ \cos(\theta)}$', labelpad=5)
plt.ylabel(r'$\mathrm{Density}$',labelpad=10)
plt.xlim(-1.01,1.01)
xmin,xmax=plt.xlim()
ymin,ymax=plt.ylim()
plt.text(xmin+(xmax-xmin)*0.01,ymax-(ymax-ymin)*0.05,label_panel)
plt.text(xmin+0.48*(xmax-xmin),ymin+0.95*(ymax-ymin),r'$N$='+str(numb_label))
plt.subplots_adjust(top=0.99,bottom=0.13,left=0.11,right=0.98,hspace=0.0,wspace=0.0)
plt.legend(numpoints=1,loc=('upper right'))

#
index_cut=dir_read.find(mc_read)
home = os.path.expanduser("~")
final_results_path = home + "/ResultsOf" + mc_read + "/"
FilePlotDensity=final_results_path+first_fragment[index_cut:-1]+last_fragment[:-1]+"-histogram-of-"+axis_read+"-vs-gFactor-preskip"+str(preskip)+"-postskip"+str(postskip)+".pdf"
print(FilePlotDensity)

#FilePlotDensity="/home/tapas/ResultsOf"+mc_read+dir_read[28:-1]+"-histogram-r-preskip"+str(preskip)+"-postskip"+str(postskip)+".eps"
plt.savefig(FilePlotDensity, dpi=200, format='pdf')
plt.show()

#python analysis_instantaneous_angles.py /scratch/tapas/linear-rotors/ENT-RotDOFs-Rpt10.05Angstrom-gFactor10.0-beta0.2Kinv-Blocks80000-Passes100-System2HF-ParticleA1-e0vsbeads-EXTENDED_ENSMBL41/ ENT 4 1 cost 41 0 
