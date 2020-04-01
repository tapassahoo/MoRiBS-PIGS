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
	'figure.figsize': (8,12),
	'axes.labelsize': size,
	'axes.titlesize': size,
	'xtick.labelsize': size*0.75,
	'ytick.labelsize': size*0.75,
	'axes.titlepad': size}
plt.rcParams.update(params)
matplotlib.rcParams.update({'font.size': size*0.75})
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"


###  Program starts from here  ########

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
system_replaced = string1[string1.find("System"):string1.find("HF")+2]
particleA_replaced = string1[string1.find("ParticleA"):string1.find("-e0")]
gfact=string1[string1.find("gFactor"):string1.find("-beta")]
#extra_file_name=string1[string1.find(mc_read)+5:string1.find(dofs_read[dofs_read1])]
#numb_beads = int(string1[string1.find("beads")+5:])
#dofs_read = string1[string1.find("moves")+6:string1.find("DOFs")]
#numb_particle = int(string1[string1.find("System")+6:string1.find("-p-H2O")])

first_fragment = dir_read[:dir_read.find("gFactor")]
last_fragment = dir_read[dir_read.find("-beta"):]

ndofs=3
if (mc_read != "PIMC"):
	beads_pos = int((numb_beads-1)/2)
	#beads_pos = numb_beads-1

fig, ax = plt.subplots()
isubplot = 1
numb_particle_list=[4,16]
system_replaced_by = {4:"System2HF", 16:"System8HF"}
particleA_replaced_by = {4:"ParticleA1", 16:"ParticleA4"}
label_panel_list = {4:"(a)", 16:"(b)"}

print("#Look at the below line")
for numb_particle in numb_particle_list:
	plt.subplot(2, 1, isubplot)
	preskip = 0
	postskip = 0

	file1=file_read.replace(system_replaced, system_replaced_by[numb_particle])
	file2=file1.replace(particleA_replaced, particleA_replaced_by[numb_particle])
	if (numb_particle == 4):
		file2=file2.replace("-WOR-", "")
	
	data_len = len(genfromtxt(file2, unpack=True, usecols=[0], skip_header=preskip, skip_footer=postskip))
	workingNdim = int(math.log(data_len)/math.log(2))
	trunc = int(data_len-2**workingNdim)
	data_len=int(numb_blocks-(preskip+trunc))
	save_data = np.zeros((numb_particle,data_len))

	label_panel = label_panel_list[numb_particle]
	if (numb_particle == 4):
		gFactList=[2.0, 4.0, 6.0, 8.0]
	elif (numb_particle == 16):
		gFactList=[0.5, 1.0, 1.5, 2.0]
	colorList=["yellow", "green", "red", "magenta"]
	j=axis_index[axis_read]
	ig=0
	for g in gFactList:
		file3=file2.replace(gfact,"gFactor"+str(g))
		print(file3)
		for i in range(numb_particle):
			ncol1 = beads_pos+i*numb_beads
			ncol = j+ncol1*ndofs
			ncol = ncol+1
			print(str(ncol)+'th column')
			save_data[i,:] = genfromtxt(file3, unpack=True, usecols=[ncol], skip_header=preskip+trunc, skip_footer=postskip)

		label_str=r'$g$='+str(g)	

		vec_plot = np.reshape(save_data, numb_particle*data_len)
		data_plot = vec_plot
		print(len(data_plot))
		plt.hist(data_plot, bins=50, density=True, alpha=0.5, stacked=True, color=colorList[ig], edgecolor='black', label=label_str)
		ig=ig+1
	print("")

	#
	if (mc_read == "ENT"):
		numb_label=int(numb_particle/2)
	plt.ylabel(r'$\mathrm{Density}$',labelpad=10)
	plt.xlim(-1.01,1.01)
	xmin,xmax=plt.xlim()
	ymin,ymax=plt.ylim()
	plt.text(xmin+(xmax-xmin)*0.04,ymax-(ymax-ymin)*0.1,label_panel)
	plt.text(xmin+0.48*(xmax-xmin),ymin+0.9*(ymax-ymin),r'$N$='+str(numb_label))
	if (isubplot == 2):
		plt.xlabel(r'$\mathrm{bins \ of \ \cos(\theta)}$', labelpad=5)
	if (isubplot != 2):
		frame1 = plt.gca()
		frame1.axes.xaxis.set_ticklabels([])
	plt.subplots_adjust(top=0.98,bottom=0.10,left=0.11,right=0.98,hspace=0.05,wspace=0.0)
	ax.tick_params(right=True,top=True,left=True,bottom=True)
	plt.legend(numpoints=1,loc='center')
	isubplot = isubplot+1

#
index_cut=dir_read.find(mc_read)
home = os.path.expanduser("~")
final_results_path = home + "/ResultsOf" + mc_read + "/"
FilePlotDensity=final_results_path+first_fragment[index_cut:-1]+last_fragment[:-1]+"-histogram-of-"+axis_read+"-vs-gFactor-preskip"+str(preskip)+"-postskip"+str(postskip)+".pdf"
print(FilePlotDensity)

#FilePlotDensity="/home/tapas/ResultsOf"+mc_read+dir_read[28:-1]+"-histogram-r-preskip"+str(preskip)+"-postskip"+str(postskip)+".eps"
plt.savefig(FilePlotDensity, dpi=100, format='pdf')
plt.show()

#python analysis_instantaneous_angles.py /scratch/tapas/linear-rotors/ENT-RotDOFs-Rpt10.05Angstrom-gFactor10.0-beta0.2Kinv-Blocks80000-Passes100-System2HF-ParticleA1-e0vsbeads-EXTENDED_ENSMBL41/ ENT 4 1 cost 41 0 

