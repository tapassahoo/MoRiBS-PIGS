import os
import sys
from subprocess import call
import math
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib.axes
#matplotlib.use('eps')
from pylab import *
import mypkg.pkgMoribs.support_without_parallel as support

def getangle(vec1,vec2):
	abs_vec1=math.sqrt(np.dot(vec1,vec1))
	abs_vec2=math.sqrt(np.dot(vec2,vec2))
	angle = np.dot(vec1,vec2)/(abs_vec1*abs_vec2)
	return angle

def matpre(eulang):
	'''
	Formation of rotation matrix
	'''
	phi=eulang[0]
	theta=eulang[1]
	chi=eulang[2]

	cp=math.cos(phi)
	sp=math.sin(phi)
	ct=math.cos(theta)
	st=math.sin(theta)
	ck=math.cos(chi)
	sk=math.sin(chi)

	rotmat=np.zeros((3,3),dtype='float')
	rotmat[0,0]=cp*ct*ck-sp*sk
	rotmat[0,1]=-cp*ct*sk-sp*ck
	rotmat[0,2]=cp*st
	rotmat[1,0]=sp*ct*ck+cp*sk
	rotmat[1,1]=-sp*ct*sk+cp*ck
	rotmat[1,2]=sp*st
	rotmat[2,0]=-st*ck
	rotmat[2,1]=st*sk
	rotmat[2,2]=ct

	return rotmat


def rottrn(rotmat,rwf,rcom):
	'''
	It transforms body-fixed frame to space-fixed frame
	rotmat: rotation matrix
	rwf: body-fixed coordinates
	rcom: center of mass coordinates
	'''
	
	rsf=np.zeros(3,dtype='float')

	for i in range(3):
		rsf[i]=rcom[i]
		for j in range(3):
			rsf[i]=rsf[i]+rotmat[i,j]*rwf[j]

	return rsf

# main function starts here #

rc('text', usetex=True)
size=24
params = {'legend.fontsize': size*0.5,
	'figure.figsize': (16,30),
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
numb_beads= int(sys.argv[4])

#
file_read = dir_read + 'results/output.xyz'
string1 =dir_read[dir_read.find(mc_read):-1]
numb_blocks = int(string1[string1.find("Blocks")+6:string1.find("-Passes")])
system_replaced = string1[string1.find("System"):string1.find("p-H2O")+5]
rpt_read=string1[string1.find("Rpt"):string1.find("Angstrom")]

first_fragment = dir_read[:dir_read.find("Rpt")]
last_fragment = dir_read[dir_read.find("Angstrom")+8:]

ndofs=3
if (mc_read != "PIMC"):
	beads_pos = int((numb_beads-1)/2)

fig, ax = plt.subplots()
isubplot = 1
system_replaced_by = {2:"System2-p-H2O", 11:"System11-p-H2O"}
label_panel_list = {0:"(a)", 1:"(b)", 2:"(c)", 3:"(d)",4:"(e)", 5:"(f)", 6:"(g)", 7:"(h)",8:"(i)", 9:"(j)"}
#rptdict = {0:[2.5,2.6,2.7],1:[2.8,2.9,3.0],2:[3.2,3.4,3.6],3:[4.0,4.5,5.0]}
#rptdict = {0:[3.1],1:[3.2],2:[3.3],3:[3.4],4:[3.5],5:[3.6],6:[3.7],7:[3.8],8:[3.9],9:[4.0]}
#rptdict = {0:[4.1],1:[4.2],2:[4.3],3:[4.4],4:[4.5],5:[4.6],6:[4.7],7:[4.8],8:[4.9],9:[5.0]}
#rptdict = {0:[5.2],1:[5.4],2:[5.6],3:[5.8],4:[6.0],5:[6.2],6:[6.4],7:[6.6],8:[6.8],9:[7.0]}
rptdict = {0:[2.8],1:[3.2],2:[3.5],3:[10.0]}

### configurations of nonlinear molecules in the body-fixed frame

br2ang=0.52917721092

## q-TIP4P/F parameters
angHOH=107.4
dOH=0.9419
ro_wf=np.zeros(3,dtype='float')
rh1_wf=np.zeros(3,dtype='float')
rh2_wf=np.zeros(3,dtype='float')
#
ang1=(angHOH*math.pi)/180.0
zH=ro_wf[2]-math.sqrt(0.5*dOH*dOH*(1.0+math.cos(ang1)))
xH=math.sqrt(dOH*dOH-(ro_wf[2]-zH)*(ro_wf[2]-zH))
#
rh1_wf[0]=xH
rh1_wf[1]=0.0
rh1_wf[2]=zH
#
rh2_wf[0]=-rh1_wf[0]
rh2_wf[1]=rh1_wf[1]
rh2_wf[2]=rh1_wf[2]
###

for isubplot in range(4):
	label_panel = label_panel_list[isubplot]
	plt.subplot(2, 2, isubplot+1)
	preskip = 0
	postskip = 0

	rptList=rptdict[isubplot]
	print(rptList)
	colorList=["yellow", "green"]
	lsList=["dashed", "dashdot"]
	j=0
	ig=0
	for rpt in rptList:

		com1=[0,0,0]
		com2=[0,0,rpt]
		
		rpt_exact = "{:3.2f}".format(rpt)
		#file_exact1 = '/home/tapas/ResultsOfExact/ground-state-theta-distribution-lanc-2-p-H2O-jmax2-Rpt'+rpt_exact+'Angstrom-grid-24-12-niter100.txt'
		#data_exact1 = np.genfromtxt(file_exact1)
		#file_exact2 = '/home/tapas/ResultsOfExact/ground-state-theta-distribution-lanc-2-p-H2O-jmax4-Rpt'+rpt_exact+'Angstrom-grid-20-20-niter200.txt'
		#data_exact2 = np.genfromtxt(file_exact2)

		file1=file_read.replace(system_replaced, system_replaced_by[numb_particle])
		file2=file1.replace(rpt_read,"Rpt"+str(rpt))
		dir_rpt=dir_read.replace(rpt_read,"Rpt"+str(rpt))
		final_dir_in_work = dir_rpt

		file_old = final_dir_in_work+"results/output.xyz_old"
		if (os.path.isfile(file_old) == True):
			file_new = final_dir_in_work+"results/output.xyz"

			if (os.path.isfile(file_new) == True):
				print(" -- Restarted data")

				if "H2O1" in open(file_new).read():
					rmstr = int(numb_particle*numb_beads+3)
					cmd1="tail -n +"+str(rmstr)+" "+final_dir_in_work+"results/output.xyz>bb"
					os.system(cmd1)
					col_data_new = np.genfromtxt("bb")
					call(["rm", "bb"])
				else:
					col_data_new = np.genfromtxt(final_dir_in_work+"results/output.xyz")
				index = int(col_data_new[0,0])
				col_data_old = np.genfromtxt(final_dir_in_work+"results/output.xyz_old")
				marged_data  = np.concatenate((col_data_old[:index-1], col_data_new), axis=0)
				aa = col_data_new[:,0]
				final_data_set = marged_data[0:int(aa[-1]),:]
			else:
				final_data_set = np.genfromtxt(final_dir_in_work+"results/output.xyz_old", skip_header=0, skip_footer=0)
		else:
			final_data_set = np.genfromtxt(final_dir_in_work+"results/output.xyz", skip_header=0, skip_footer=0)

		data_len = len(final_data_set[:,0])
		nlen = int(len(final_data_set[:,0]))
	
		workingNdim = int(math.log(data_len)/math.log(2))
		trunc = int(data_len-2**workingNdim)
		data_len=int(data_len-(preskip+trunc))
		print(data_len)
	
		cost_data = np.zeros((numb_particle,data_len))
		phi_data = np.zeros((numb_particle,data_len))
		chi_data = np.zeros((numb_particle,data_len))

		print(file2)
		for i in range(numb_particle):
			ncol1 = beads_pos+i*numb_beads
			ncol = j+ncol1*ndofs
			ncol = ncol+1
			print(str(ncol)+'th column')
			cost_data[i,:] = final_data_set[(preskip+trunc):(nlen-postskip),ncol]
			ncol_phi=ncol+1
			phi_data[i,:] = final_data_set[(preskip+trunc):(nlen-postskip),ncol_phi]
			ncol_chi=ncol+2
			chi_data[i,:] = final_data_set[(preskip+trunc):(nlen-postskip),ncol_chi]
			

		label_str=r'$r$='+str(rpt)	

		if (numb_particle == 2):
			vec_cost1 = np.reshape(cost_data[0,:], (numb_particle-1)*data_len)
			vec_phi1 = np.reshape(phi_data[0,:], (numb_particle-1)*data_len)
			vec_chi1 = np.reshape(chi_data[0,:], (numb_particle-1)*data_len)
			vec_cost2 = np.reshape(cost_data[1,:], (numb_particle-1)*data_len)
			vec_phi2 = np.reshape(phi_data[1,:], (numb_particle-1)*data_len)
			vec_chi2 = np.reshape(chi_data[1,:], (numb_particle-1)*data_len)

		eulang1=np.zeros(3,dtype='float')
		eulang2=np.zeros(3,dtype='float')
		
		angle1_data = np.zeros(data_len)
		angle2_data = np.zeros(data_len)
		angle3_data = np.zeros(data_len)
		angle4_data = np.zeros(data_len)


		index_cut=dir_read.find(mc_read)
		home = os.path.expanduser("~")
		final_results_path = home + "/ResultsOf" + mc_read + "/"
		FileXYZ1=final_results_path+first_fragment[index_cut:-1]+last_fragment[:-1]+"-first-rotor-vs-gFactor-preskip"+str(preskip)+"-postskip"+str(postskip)+"-Rpt"+str(rpt)+"Angstrom.xyz"
		FileXYZ2=final_results_path+first_fragment[index_cut:-1]+last_fragment[:-1]+"-second-rotor-vs-gFactor-preskip"+str(preskip)+"-postskip"+str(postskip)+"-Rpt"+str(rpt)+"Angstrom.xyz"
		

		f= open(FileXYZ1,"w+")
		#g= open(FileXYZ2,"w+")

		for l in range(len(vec_cost1)):
			eulang1[0]=vec_phi1[l]
			eulang1[1]=math.acos(vec_cost1[l])
			eulang1[2]=vec_chi1[l]
			rotmat1=matpre(eulang1)
			ro_1sf=rottrn(rotmat1, ro_wf, com1)
			rh1_1sf=rottrn(rotmat1, rh1_wf, com1)
			rh2_1sf=rottrn(rotmat1, rh2_wf, com1)

			eulang2[0]=vec_phi2[l]
			eulang2[1]=math.acos(vec_cost2[l])
			eulang2[2]=vec_chi2[l]
			rotmat2=matpre(eulang2)
			ro_2sf=rottrn(rotmat2, ro_wf, com2)
			rh1_2sf=rottrn(rotmat2, rh1_wf, com2)
			rh2_2sf=rottrn(rotmat2, rh2_wf, com2)
			str_O="O"+" "+str(ro_1sf[0])+" "+str(ro_1sf[1])+" "+str(ro_1sf[2])+'\n'
			str_H1="H"+" "+str(rh1_1sf[0])+" "+str(rh1_1sf[1])+" "+str(rh1_1sf[2])+'\n'
			str_H2="H"+" "+str(rh2_1sf[0])+" "+str(rh2_1sf[1])+" "+str(rh2_1sf[2])+'\n'
			f.write('3\n')
			f.write('\n')
			f.write(str_O)
			f.write(str_H1)
			f.write(str_H2)

			str_O="O"+" "+str(ro_2sf[0])+" "+str(ro_2sf[1])+" "+str(ro_2sf[2])+'\n'
			str_H1="H"+" "+str(rh1_2sf[0])+" "+str(rh1_2sf[1])+" "+str(rh1_2sf[2])+'\n'
			str_H2="H"+" "+str(rh2_2sf[0])+" "+str(rh2_2sf[1])+" "+str(rh2_2sf[2])+'\n'
			#g.write('3\n')
			#g.write('\n')
			f.write(str_O)
			f.write(str_H1)
			f.write(str_H2)

			o_1sf_h1_1sf=np.subtract(ro_1sf,rh1_1sf)
			o_1sf_h2_1sf=np.subtract(ro_1sf,rh2_1sf)
			o_2sf_h1_1sf=np.subtract(ro_2sf,rh1_1sf)
			o_2sf_h2_1sf=np.subtract(ro_2sf,rh2_1sf)

			o_2sf_h1_2sf=np.subtract(ro_2sf,rh1_2sf)
			o_2sf_h2_2sf=np.subtract(ro_2sf,rh2_2sf)
			o_1sf_h1_2sf=np.subtract(ro_1sf,rh1_2sf)
			o_1sf_h2_2sf=np.subtract(ro_1sf,rh2_2sf)

			ang_oh1o2=getangle(o_1sf_h1_1sf,o_2sf_h1_1sf)
			ang_oh2o2=getangle(o_1sf_h2_1sf,o_2sf_h2_1sf)
			ang_oh1o1=getangle(o_2sf_h1_2sf,o_1sf_h1_2sf)
			ang_oh2o1=getangle(o_2sf_h2_2sf,o_1sf_h2_2sf)
			
			angle1_data[l]=ang_oh1o2
			angle2_data[l]=ang_oh2o2
			angle3_data[l]=ang_oh1o1
			angle4_data[l]=ang_oh2o1

		f.close()
		#g.close()
			
		#angle_oho1=np.concatenate((angle1_data))
		#angle_oho2=np.concatenate((angle3_data, angle4_data))
		angle_oho=np.concatenate((angle1_data, angle2_data, angle3_data, angle4_data))

		#plt.hist(angle2_data, bins=50, density=True, alpha=0.5, stacked=True, color=colorList[0], edgecolor='black', label='rotor1')
		#plt.hist(angle_oho2, bins=50, density=True, alpha=0.5, stacked=True, color=colorList[1], edgecolor='black', label='rotor2')
		plt.hist(angle_oho, bins=50, density=True, alpha=0.5, stacked=True, color='yellow', edgecolor='black', label=label_str+r'$\mathrm{\AA}$')

		ig=ig+1

	#
	numb_label=int(numb_particle)
	plt.xlim(-1.01,1.01)
	xmin,xmax=plt.xlim()
	ymin,ymax=plt.ylim()
	plt.text(xmin+(xmax-xmin)*0.04,ymax-(ymax-ymin)*0.1,label_panel)
	plt.text(xmin+0.45*(xmax-xmin),ymin+0.9*(ymax-ymin),r'$N$='+str(numb_label))
	#plt.text(xmin+0.4*(xmax-xmin),ymin+0.9*(ymax-ymin),r'$N$='+str(numb_label)+'; '+label_str+r'$\mathrm{\AA}$')

	plt.ylabel(r'$p(\cos(\alpha))$',labelpad=5)

	if ((isubplot == 8) or (isubplot == 9)):
		plt.xlabel(r'$\mathrm{bins \ of \ \cos(\alpha)}$', labelpad=5)
	else:
		frame1 = plt.gca()
		frame1.axes.xaxis.set_ticklabels([])

	if (numb_particle == 2):
		plt.subplots_adjust(top=0.99,bottom=0.03,left=0.06,right=0.99,hspace=0.06,wspace=0.17)
	if (numb_particle == 11):
		plt.subplots_adjust(top=0.99,bottom=0.03,left=0.06,right=0.99,hspace=0.06,wspace=0.17)
	ax.tick_params(right=True,top=True,left=True,bottom=True)
	plt.legend(numpoints=1,loc='center')
	plt.minorticks_on()
	plt.tick_params(axis="both", direction="out", which="minor", right=True, top=True, length=2)
	plt.tick_params(axis="both", direction="out", which="major", right=True, top=True, length=5)

#
index_cut=dir_read.find(mc_read)
home = os.path.expanduser("~")
final_results_path = home + "/ResultsOf" + mc_read + "/"
FilePlotDensity=final_results_path+first_fragment[index_cut:-1]+last_fragment[:-1]+"-oho-angle-vs-gFactor-preskip"+str(preskip)+"-postskip"+str(postskip)+"-Rpt"+str(rptdict[0][0])+"-"+str(rptdict[9][0])+"Angstrom.pdf"
print(FilePlotDensity)

plt.savefig(FilePlotDensity, dpi=100, format='pdf')
plt.show()
#
#python analysis_instantaneous_angles_preskip.py /Users/tsahoo/nonlinear-rotors/PIGS-qTIP4P-RotDOFs-Rpt2.5Angstrom-beta0.1Kinv-Blocks20000-Passes200-System2-p-H2O-e0vsbeads101/ PIGS 2 1 cost 101 0
