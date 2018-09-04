#!/usr/bin/env python
import numpy as np
from numpy import *
import numpy as np
import matplotlib
#matplotlib.use('eps')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from subprocess import call
#from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter
import support
import os
from pylab import *
import matplotlib.axes
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


def	FigureENT(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, TypePlot, beadsRef):
	font       = 28
	fontlegend = font/2.0
	preskip    = 0
	postskip   = 0

	plotnum    = 2

	BConstant        = support.GetBconst(molecule_rot)  # in wavenumber
	Units          	 = support.GetUnitConverter()
	BConstantK       = BConstant*Units.CMRECIP2KL

	if (((TypePlot == "RFACTOR") or (TypePlot == "GFACTOR")) and variableName == "tau"):
		font = 28
		fontlegend = font/2
		fig        = plt.figure(figsize=(8, 6))
		plt.grid(True)

		var = beadsRef
		FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, var)
		if (TypePlot == "RFACTOR"):
			FilePlotEntropy  = FilePlotName.SaveEntropyRFAC+str(plotnum)+".eps"
		if (TypePlot == "GFACTOR"):
			FilePlotEntropy  = FilePlotName.SaveEntropyGFAC+str(plotnum)+".eps"
		outfileEntropy       = FilePlotEntropy
		call(["rm", FilePlotEntropy])
#
		if (plotnum != 2):
			plt.axhline(y=log(2.0), color='blue', lw = 2.0, linestyle='-', label = 'ln(2)')
#
		if (plotnum != 2):
			nn = [2,4,6]
		else:
			nn = [8,16,32]
		for numbmolecules in nn:
			particleA = int(numbmolecules/2)
			if (numbmolecules == 2):
				beadsRef = 101
				DList  = [1.0+0.5*i for i in range(7)]
				DList  += [4.5, 5.0, 5.5, 6.0]
				numbblocks = 50000
				postskip1  = 30000
			if (numbmolecules == 4):
				DList  = [1.0+0.5*i for i in range(7)]
				numbblocks = 20000
			if (numbmolecules == 6):
				DList  = [1.0+0.25*i for i in range(10)]
				numbblocks = 20000
			if (numbmolecules == 8):
				DList  = [0.9106, 1.1152, 1.2878, 1.4398, 1.5772, 1.7036, 1.8212, 1.9317, 2.0362]
				#DList  = [1.5772, 1.7036, 1.8212, 1.9317, 2.0362]
				numbblocks = 20000
				preskip1   = 8000
				postskip1  = 10000
			if (numbmolecules == 16):
				#DList  = [1.5772, 1.7036, 1.8212, 1.9317, 2.0362]
				DList  = [0.9106, 1.1152, 1.2878, 1.4398, 1.5772, 1.7036, 1.8212, 1.9317, 2.0362]
				numbblocks = 10000
				preskip1   = 8000
				postskip1  = 0
			if (numbmolecules == 32):
				#DList  = [1.5772, 1.7036, 1.8212, 1.9317, 2.0362]
				DList  = [0.9106, 1.1152, 1.2878, 1.4398, 1.5772, 1.7036, 1.8212, 1.9317, 2.0362]
				numbblocks = 10000
				preskip1   = 6000
				postskip1  = 0
		
			RFactorPlot      = np.zeros(len(DList))
			entropy1Plot     = np.zeros(len(DList))
			purity1Plot      = np.zeros(len(DList))
			err_entropy1Plot = np.zeros(len(DList))
			err_purity1Plot  = np.zeros(len(DList))
			entropy2Plot     = np.zeros(len(DList))
			entropy3Plot     = np.zeros(len(DList))

			iii = 0
			for dipolemoment in DList:
				RFactorList = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
				RFactor = RFactorList[0]
				if (dipolemoment > 4.0):
					if (numbmolecules == 2):
						numbblocks = 20000
						postskip1  = 0
					else:
						numbblocks = 20000
				FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef)
				FileToBePlotEntropy = FilePlotName.SaveEntropy+".txt"
				FileToBePlotDIAG    = FilePlotName.SaveEntropyDIAG+".txt"

				beads1, var1, purity1, entropy1, err_purity1, err_entropy1 = genfromtxt(FileToBePlotEntropy,unpack=True, usecols=[0, 1, 4, 5, 8, 9], skip_header=preskip, skip_footer=postskip)
				if (numbmolecules <= 6):
					RFactor, energy3, entropy3                             = loadtxt(FileToBePlotDIAG,unpack=True, usecols=[0,1,2])
					
				if (TypePlot == "GFACTOR"):
					RFactorPlot[iii] = 1.0/(RFactor*RFactor*RFactor)
				if (TypePlot == "RFACTOR"):
					RFactorPlot[iii] = RFactor
				if (plotnum != 2):
					entropy3Plot[iii] = entropy3
				if ((numbmolecules == 4) and (dipolemoment ==4.0)):
					beadsRef = 61
					beadsRef1 = 81
				elif ((numbmolecules == 6) and (dipolemoment == 3.25)):
					beadsRef = 21
					beadsRef1 = 81
				elif ((numbmolecules == 6) and (dipolemoment == 3.5)):
					beadsRef = 101
					beadsRef1 = 81
				else:
					beadsRef = 81
					beadsRef1 = 81
				if (numbmolecules == 2):
					beadsRef = 101
					beadsRef1 = beadsRef
				if ((numbmolecules == 2) and (dipolemoment == 6.0)):
					beadsRef = 41
					beadsRef1 = 41
				if (numbmolecules == 8):
					beadsRef = 21
					beadsRef1 = beadsRef
				if (numbmolecules == 16):
					beadsRef = 21
					beadsRef1 = beadsRef
				if (numbmolecules == 32):
					beadsRef = 21
					beadsRef1 = beadsRef
	
				if (np.isscalar(entropy1) == True):
					entropy1Plot[iii]     = entropy1
					err_entropy1Plot[iii] = err_entropy1
					purity1Plot[iii]      = purity1
					err_purity1Plot[iii]  = err_purity1
				else:
					ii = 0
					for i in beads1:
						indexi =int(i+0.5)
						beads = indexi
						if beads == beadsRef:
							entropy1Plot[iii]     = entropy1[ii]
							purity1Plot[iii]      = purity1[ii]
							err_purity1Plot[iii]  = err_purity1[ii]
		
						ii += 1
					ii1 = 0
					for i in beads1:
						indexi =int(i+0.5)
						beads = indexi
						if beads == beadsRef1:
							err_entropy1Plot[iii] = err_entropy1[ii1]
		
						ii1 += 1
				iii += 1
		
			print("S2:  PIGS "+str(numbmolecules))
			print(entropy1Plot)
			print("S2:  ED ")
			print(entropy3Plot)
#
			plotEntropyENT1(numbmolecules,RFactorPlot, entropy1Plot, err_entropy1Plot, variableName, RFactorPlot, entropy2Plot, entropy3Plot, font, TypePlot)
		ymin, ymax = plt.ylim()
		xmin, xmax = plt.xlim()
		if plotnum == 1:
			plt.ylabel(r'$S_{2}$', fontsize = font)
			if ymin < 0.0:
				plt.ylim(0.0,0.82)
			if xmin < 0.0:
				plt.xlim(0.0,xmax)
			plt.xticks(np.arange(0, 9, step=2))
			plt.yticks(np.arange(0, 0.9, step=0.2))
			Text1 = "(a)"
		if plotnum == 2:
			plt.xlim(0.192,1.01)
			plt.ylim(-0.46,0.64)
			plt.xticks(np.arange(0.2, 1.1, step=0.2))
			plt.yticks(np.arange(-0.4, 0.6, step=0.2))
			Text1 = "(b)"
			plt.ylabel(r'$S_{2}$', fontsize = font, labelpad=-14)
		Text2 = ""
		if Text1:
			PlotLabel(Text1, Text2,font,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
			
		plt.subplots_adjust(top=0.97, bottom=0.14, left=0.14, right=0.98, hspace=0.0, wspace=0.0)
		if (plotnum !=2):
			plt.legend(bbox_to_anchor=(0.68, 0.75), loc=2, borderaxespad=1., shadow=True, fontsize = fontlegend)
		else:
			plt.legend(bbox_to_anchor=(0.98, 0.35), borderaxespad=1., shadow=True, fontsize = fontlegend)
		plt.savefig(outfileEntropy, dpi = 200, format = 'eps')

		call(["open", outfileEntropy])
		#plt.show()

	if (TypePlot == "S2"):
		font = 28
		fontlegend = font/2
		iFrame = 1
		iFigLabel = 0
		FigureLabel  = [r'$\mathrm{(a)}$',r'$\mathrm{(b)}$']
		variableList = ["beta","tau"]
		for variableName in variableList:
			if variableName == "beta":
				parameterName = "tau"
				parameter     = 0.005
				postskip      = 8
			if variableName == "tau":
				parameterName = "beta"
				parameter     = 0.2
				postskip      = 0
##
			fig        = plt.figure(figsize=(8, 6))
			plt.grid(True)
##	
			FilePlotName      = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef)
			FileToBePlot   	  = FilePlotName.SaveEntropy+".txt"
			FileToBePlotMM    = FilePlotName.SaveEntropyMM+".txt"
			FileToBePlotDIAG  = FilePlotName.SaveEntropyDIAG+".txt"
			print(FileToBePlot)
			print(FileToBePlotMM)
			print(FileToBePlotDIAG)
##
			FilePlot     = FilePlotName.SaveEntropy+".eps"
			outfile      = FilePlot
			call(["rm", FilePlot])
			print(FilePlot)
##
			var2       = 0.0
			entropy2   = [0.0, 0.0]
			entropy3   = 0.0
###================data reading begin==================###
			var1, purity1, entropy1, err_purity1, err_entropy1 = genfromtxt(FileToBePlot,unpack=True, usecols=[1, 4, 5, 8, 9], skip_header=preskip, skip_footer=postskip)
			if (numbmolecules <= 4):
				beads, entropy2                                = genfromtxt(FileToBePlotMM,unpack=True, usecols=[0,2], skip_header=preskip, skip_footer=postskip)
				if (variableName == "tau"):
					var2 = parameter/(beads-1.0)
		
				if (variableName == "beta"):
					var2 = parameter*(beads-1.0)
	
			if (numbmolecules <= 6):
				RFactor, entropy3                              = loadtxt(FileToBePlotDIAG,unpack=True, usecols=[0, 2])
###================data reading end====================###
			print(entropy1)
			print(entropy2)
			if variableName == "tau":
				var2 = var2[1:]
				entropy2 = entropy2[1:]
			
			#plt.subplot(1, 2, iFrame)
			iFrame = iFrame + 1
			plotEntropyENT(var1, entropy1, err_entropy1, variableName, var2, entropy2, entropy3, font,font)
			plt.ylabel(r'$S_{2}$', fontsize = font)
			ymin, ymax = plt.ylim()
			xmin, xmax = plt.xlim()
			Text1 = FigureLabel[iFigLabel]
			iFigLabel = iFigLabel +1
			RGFactor = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
			gFactor = RGFactor[1]
			arg     = "%3.2f" % gFactor
			Text2 = str(arg)
			PlotLabel(Text1, Text2,font,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
			plt.legend(bbox_to_anchor=(0.80, 0.50), loc=2, borderaxespad=0., shadow=True, fontsize = fontlegend)
			plt.subplots_adjust(top=0.98, bottom=0.16, left=0.18, right=0.95, hspace=0.0, wspace=0.)
			plt.savefig(outfile, dpi = 200, format = 'eps')
	
			call(["open", outfile])
			#plt.show()

def	FigureENTCOMBINE(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, TypePlot, beadsRef):
	font       = 28
	fontlegend = font/2.0
	preskip    = 0
	postskip   = 0

	if (variableName == "tau"):
		fig        = plt.figure(figsize=(8, 6))
		plt.grid(True)

		var = beadsRef
		FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, var)
		FilePlotEntropy  = FilePlotName.SaveEntropyCOMBO+".eps"
		outfileEntropy   = FilePlotEntropy
		call(["rm", FilePlotEntropy])
		print(outfileEntropy)
#
		DList  = [1.0+0.5*i for i in range(7)]
		DList += [4.5, 5.0, 5.5, 6.0]
		labelIndex = 0
		TypeList = ["SWAPTOUNSWAP", "BROKENPATH"]
		for ENT_TYPE in TypeList:
			beadsRef1 = beadsRef
			if (ENT_TYPE == "SWAPTOUNSWAP"):
				numbblocks = 50000
				postskip1   = 30000
			if (ENT_TYPE == "BROKENPATH"):
				numbblocks = 20000
				postskip1   = 0
		
			RFactorPlot      = np.zeros(len(DList))
			entropy1Plot     = np.zeros(len(DList))
			err_entropy1Plot = np.zeros(len(DList))
			entropy2Plot     = np.zeros(len(DList))

			iii = 0
			for dipolemoment in DList:
				RFactorList = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
				if (dipolemoment > 4.0):
					numbblocks = 20000
					postskip1   = 0
					if ((ENT_TYPE == "SWAPTOUNSWAP") and (dipolemoment == 6.0)):
						beadsRef1 = 41
						
				FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef1)
				FileToBePlotEntropy = FilePlotName.SaveEntropy+".txt"
			
				if (ENT_TYPE == "SWAPTOUNSWAP"):
					FileToBePlotDIAG    = FilePlotName.SaveEntropyDIAG+".txt"
					beads1, var1, purity1, entropy1, err_purity1, err_entropy1 = genfromtxt(FileToBePlotEntropy,unpack=True, usecols=[0, 1, 4, 5, 8, 9], skip_header=preskip, skip_footer=postskip)
					entropyDIAG                                                = loadtxt(FileToBePlotDIAG,unpack=True, usecols=[2])
				if (ENT_TYPE == "BROKENPATH"):
					beads1, var1, entropy1, err_entropy1                       = genfromtxt(FileToBePlotEntropy,unpack=True, usecols=[0, 1, 4, 7], skip_header=preskip, skip_footer=postskip)
				RFactorPlot[iii] = RFactorList[1]
	
				ii = 0
				for i in beads1:
					indexi =int(i+0.5)
					beads = indexi
					if beads == beadsRef1:
						entropy1Plot[iii]     = entropy1[ii]
						entropy2Plot[iii]     = entropyDIAG
						err_entropy1Plot[iii] = err_entropy1[ii]
					ii += 1
				iii += 1
		
			print("S2:  PIGS ")
			print(entropy1Plot)
			print(RFactorPlot)
#
			colorList = ['red', 'blue']
			lsList = ['--', '-.']
			markerList = ['^','v']
			labelList  = ['Swap+Unswap grand ensemble','Broken path ensemble']

			if (ENT_TYPE == "SWAPTOUNSWAP"):
				plt.plot(RFactorPlot, entropy2Plot, color = 'black', ls = '-', linewidth=1,  marker = 'o', markersize = 9, label = 'ED')
			plt.errorbar(RFactorPlot, entropy1Plot, yerr=err_entropy1Plot, color = colorList[labelIndex], ls = lsList[labelIndex], linewidth=1,  marker = markerList[labelIndex], markersize = 8, label = labelList[labelIndex])

			labelIndex += 1

			ymin, ymax = plt.ylim()
			plt.ylim(-0.001,0.901)
			xmin, xmax = plt.xlim()
			plt.xlim(0,9)
			Text1 = ""
			Text2 = ""
			if Text1:
				PlotLabel(Text1, Text2,font,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
			plt.xticks(np.arange(0, 10, step=1),fontsize=font, rotation=0)
			plt.yticks(np.arange(0.0, 0.91, step=0.1),fontsize=font, rotation=0)

		plt.ylabel(r'$S_{2}$', fontsize = font)
		plt.xlabel(r'$g$', fontsize = font, labelpad=-3)

		plt.subplots_adjust(top=0.97, bottom=0.14, left=0.14, right=0.98, hspace=0.0, wspace=0.4)
		plt.legend(bbox_to_anchor=(0.40, 0.45), loc=2, borderaxespad=1., shadow=True, fontsize = fontlegend)
		plt.savefig(outfileEntropy, dpi = 200, format = 'eps')

		call(["open", outfileEntropy])

def plotEntropyENT1(numbmolecules,var, val, err_val, variableName, var1, val1, val2, font, TypePlot):
	if (numbmolecules == 2):
		plt.errorbar(var, val, yerr=err_val, color = "red", ls = '-', linewidth=1,  marker = "o", markersize = 8, label = 'PIGS: N=2')
		plt.plot(var, val2, color = 'black', ls = 'None', linewidth=1, marker = "o", markersize = 10, label = 'ED: N=2')
		#plt.plot(var, val2, color = 'black', ls = '-', linewidth=3, marker = "o", markersize = 12, label = 'ED: 2 HF')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "o", markersize = 10, label = 'DMRG: 2 HF')
		#plt.plot(var, val1, linestyle = 'None', color = 'mediumpurple', marker = "o", markersize = 10, label = 'MM: 2 HF')

	if (numbmolecules == 4):
		plt.errorbar(var, val, yerr=err_val, color = "red", ls = '-', linewidth=1,  marker = "s", markersize = 8, label = 'PIGS: N=4')
		plt.plot(var, val2, color = 'black', ls = 'None', linewidth=1, marker = "s", markersize = 10, label = 'ED: N=4')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "s", markersize = 8, label = 'DMRG: 4 HF')

	if (numbmolecules == 6):
		plt.errorbar(var, val, yerr=err_val, color = "red", ls = '-', linewidth=1,  marker = "v", markersize = 8, label = 'PIGS: N=6')
		plt.plot(var, val2, color = 'black', ls = 'None', linewidth=1, marker = "v", markersize = 10, label = 'ED: N=6')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "8", markersize = 7, label = 'DMRG: 6 HF')

	if (numbmolecules == 8):
		plt.errorbar(var, val, yerr=err_val, color = "red", ls = '-', linewidth=1,  marker = "p", markersize = 8, label = 'PIGS: N=8')
		#plt.plot(var, val2, color = 'black', ls = 'None', linewidth=1, marker = "v", markersize = 10, label = 'ED: N=6')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "8", markersize = 7, label = 'DMRG: 6 HF')

	if (numbmolecules == 16):
		plt.errorbar(var, val, yerr=err_val, color = "blue", ls = '-', linewidth=1,  marker = "8", markersize = 8, label = 'PIGS: N=16')
		#plt.plot(var, val2, color = 'black', ls = 'None', linewidth=1, marker = "v", markersize = 10, label = 'ED: N=6')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "8", markersize = 7, label = 'DMRG: 6 HF')

	if (numbmolecules == 32):
		plt.errorbar(var, val, yerr=err_val, color = "green", ls = '-', linewidth=1,  marker = "^", markersize = 8, label = 'PIGS: N=32')
		#plt.plot(var, val2, color = 'black', ls = 'None', linewidth=1, marker = "v", markersize = 10, label = 'ED: N=6')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "8", markersize = 7, label = 'DMRG: 6 HF')

	ymin, ymax = plt.ylim()
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy
	plt.xticks(fontsize=font, rotation=0)
	plt.yticks(fontsize=font, rotation=0)

	if (TypePlot == "RFACTOR"):
		plt.xlabel(r'$R$', fontsize = font, labelpad=0)
	if (TypePlot == "GFACTOR"):
		plt.xlabel(r'$g$', fontsize = font, labelpad=-3)

def plotEntropyENT(var, val, err_val, variableName, var1, val1, val2, font,fontlegend):
	plt.errorbar(var, val, yerr=err_val, color = 'red', ls = '-', linewidth=1,  marker = "o", markersize = 7, label = 'PIGS')
	if (val1[0] != 0.0):
		plt.plot(var1, val1, linestyle = '--', linewidth=2,color = 'blue', marker = "v", markersize = 7, label = 'MM')

	if (variableName == "tau"):
		label_xtics = [0.00]
		plt.xlim(0,0.0201)
		if (val2 != 0.0):
			plt.axhline(y=val2, color='black', lw = 2.0, linestyle='-', label = 'ED')
		plt.xticks(np.arange(0.0, 0.021, step=0.005),fontsize=fontlegend, rotation=0)
		plt.ylim(0.0399,0.11001)
		plt.yticks(np.arange(0.04, 0.11, step=0.02),fontsize=fontlegend, rotation=0)
	else:
		plt.xlim(0.0195,0.12001)
		plt.xticks(np.arange(0.02, 0.13, step=0.02),fontsize=fontlegend, rotation=0)
		plt.ylim(0.02499,0.06001)
		plt.yticks(np.arange(0.025, 0.060, step=0.01),fontsize=fontlegend, rotation=0)

	if (variableName == "beta"):
		plt.xlabel(r'$\beta \ \  (\mathrm{K^{-1}})$', fontsize = font)
	if (variableName == "tau"):
		plt.xlabel(r'$\tau \ \  (\mathrm{K^{-1}})$', fontsize = font)

def FigureCorrelation(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef, RefPointList):

	font          = 40
	fontlegend    = font/2.0

	colorList  = ['r', 'g', 'b', 'm']
	markerList = ['o', 's', 'p', '<', '>', 's', '8', 'p']
	lsList     = ['-','--','-.',':']

	TextLabel  = [r'$\mathrm{(a)}$']
	TextLabel += [r'$\mathrm{(b)}$']
	TextLabel += [r'$\mathrm{(c)}$']
	TextLabel += [r'$\mathrm{(d)}$']
	TextLabel += [r'$\mathrm{(e)}$']
	TextLabel += [r'$\mathrm{(f)}$']
	TextLabel += [r'$\mathrm{(g)}$']
	TextLabel += [r'$\mathrm{(h)}$']
	TextLabel += [r'$\mathrm{(i)}$']
	TextLabel += [r'$\mathrm{(j)}$']
	DList  = [1.0+i for i in range(4)]
	DList  = [1.25, 2.25,  3.25]
	iTextLabel = 0
	for dipolemoment in DList:
		RGFactor = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
		gFactor = RGFactor[1]
		arg     = "%3.2f" % gFactor

		TypeCorrList = ["TotalCorr","ZCorr","XYCorr"]
		for TypeCorr in TypeCorrList:
			FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef)
			if TypeCorr    == "TotalCorr":
				FileToBePlot   	  = FilePlotName.SaveTotalCorr+".txt"
				FilePlot          = FilePlotName.SaveTotalCorr+"-ithRotor"+str(RefPointList[0])+".eps"
			elif TypeCorr  == "ZCorr":
				FileToBePlot   	  = FilePlotName.SaveZCorr+".txt"
				FilePlot          = FilePlotName.SaveZCorr+"-ithRotor"+str(RefPointList[0])+".eps"
			elif TypeCorr  == "XYCorr":
				FileToBePlot   	  = FilePlotName.SaveXYCorr+".txt"
				FilePlot          = FilePlotName.SaveXYCorr+"-ithRotor"+str(RefPointList[0])+".eps"
			print(FilePlot)
			call(["rm", FilePlot])
			datacorr      = genfromtxt(FileToBePlot)
			if (datacorr.ndim == 2):
				index = 0
				for i in range(len(datacorr)):
					if datacorr[i,0] == beadsRef:
						findex = index
						break
					index = index+1
			
			FuncCorr      = np.zeros((numbmolecules,numbmolecules))
			ErrorFuncCorr = np.zeros((numbmolecules,numbmolecules))
			ii            = 0
			for i in range(numbmolecules):
				for j in range(i,numbmolecules):
			
					nc = 2+(2*ii)
					nec = nc+1 
					if datacorr.ndim == 2:
						FuncCorr[i,j] = datacorr[findex,nc]
						ErrorFuncCorr[i,j] = datacorr[findex,nec]
					else:
						FuncCorr[i,j] = datacorr[nc]
						ErrorFuncCorr[i,j] = datacorr[nec]
					if (j != i):
						FuncCorr[j,i] = FuncCorr[i,j]
						ErrorFuncCorr[j,i] = ErrorFuncCorr[i,j]
					ii = ii+1
			fig = plt.figure(figsize=(8, 6))
			iRef = 0
			for RefPoint in RefPointList:
				val1 = np.arange(1,numbmolecules+1)
				val2 = np.zeros(numbmolecules)
				val3 = np.zeros(numbmolecules)
				for j in range(numbmolecules):
					val2[j] = FuncCorr[RefPoint,j]
					val3[j] = ErrorFuncCorr[RefPoint,j]
		
				value       = "%3.2f" % val2[0]
				#plt.errorbar(val1[1:], val2[1:], yerr=val3[1:], color = colorList[iRef], ls = lsList[iRef], linewidth=3, marker = markerList[iRef], markersize = 20)  #excluding selfcorrelation
				plt.errorbar(val1, val2, yerr=val3, color = colorList[iRef], ls = lsList[iRef], linewidth=3, marker = markerList[iRef], markersize = 20)  #excluding selfcorrelation

				ymin, ymax  = plt.ylim()
				midpointy   = 0.5*(ymax-ymin)
				deltay      = midpointy*0.15
				xmin, xmax  = plt.xlim()
				plt.xlim(xmin-0.1,xmax+0.1)
				midpointx   = 0.5*(xmax-xmin)
				deltax      = midpointx*0.15
				textpositionx = xmin+midpointx-0.25*midpointx
				textpositiony = ymin+midpointy

				Text1 = TextLabel[iTextLabel]
				Text2 = str(arg)
				if (TypeCorr == "TotalCorr"):
					plt.ylabel(r'$\mathrm{C}_{'+str(RefPoint+1)+',j}$', fontsize = font)
					#plt.text(textpositionx-2*deltax, textpositiony-6*deltay, r'$\mathrm{C}_{'+str(RefPoint+1)+',1} = $'+str(value), fontsize=font)
					if Text1:
						plt.text(textpositionx-(4.75*deltax), textpositiony+(5.5*deltay), Text1, fontsize=font)

					if Text2:
						plt.text(textpositionx+(3.5*deltax), textpositiony+5.5*deltay, r'$g = $'+Text2, fontsize=font)
		
				if (TypeCorr == "ZCorr"):
					plt.ylabel(r'$\mathrm{C}^{\parallel}_{'+str(RefPoint+1)+',j}$', fontsize = font)
					#if dipolemoment == 1.25:
					#	plt.text(textpositionx-0*deltax, textpositiony-2*deltay,  r'$\mathrm{C}^{\parallel}_{'+str(RefPoint+1)+',1} = $'+str(value), fontsize=font)
					#else:
					#	plt.text(textpositionx-2*deltax, textpositiony-6*deltay,  r'$\mathrm{C}^{\parallel}_{'+str(RefPoint+1)+',1} = $'+str(value), fontsize=font)
					if Text1:
						'''
						if dipolemoment == 1.25:
							plt.text(textpositionx-(4.05*deltax), textpositiony+(5.5*deltay), Text1, fontsize=font)
						else:
							plt.text(textpositionx-(4.75*deltax), textpositiony+(3.5*deltay), Text1, fontsize=font)
						'''
						plt.text(textpositionx-(4.75*deltax), textpositiony+(5.5*deltay), Text1, fontsize=font)

					if Text2:
						plt.text(textpositionx+(3.5*deltax), textpositiony+5.5*deltay, r'$g = $'+Text2, fontsize=font)
		
				if (TypeCorr == "XYCorr"):
					plt.ylabel(r'$\mathrm{C}^{\bot}_{'+str(RefPoint+1)+',j}$', fontsize = font)
					plt.axhline(y=0.0, color='black', lw = 3.0, linestyle='--')
					#plt.text(textpositionx-2*deltax, textpositiony-6*deltay, r'$\mathrm{C}^{\bot}_{'+str(RefPoint+1)+',1} = $'+str(value), fontsize=font)
					if Text1:
						plt.text(textpositionx-(4.75*deltax), textpositiony+(5.5*deltay), Text1, fontsize=font)
						#plt.text(textpositionx-(4.75*deltax), textpositiony+(3.5*deltay), Text1, fontsize=font)

					if Text2:
						#plt.text(textpositionx+(3.5*deltax), textpositiony+3.0*deltay, r'$g = $'+Text2, fontsize=font)
						plt.text(textpositionx+(3.5*deltax), textpositiony+5.5*deltay, r'$g = $'+Text2, fontsize=font)
		

				iRef += 1

			#plt.grid(True)

			plt.xlabel(r'$j$', fontsize = font)
	
			iTextLabel += 1

			plt.xticks(fontsize=font, rotation=0)
			plt.yticks(fontsize=font, rotation=0)
			plt.subplots_adjust(top=0.96, bottom=0.23, left=0.30, right=0.95, hspace=0.6, wspace=1.0)
			#plt.subplots_adjust(top=0.96, bottom=0.18, left=0.22, right=0.95, hspace=0.6, wspace=1.0)
			#plt.legend(bbox_to_anchor=(0.58, 0.99), loc=2, borderaxespad=0., shadow=True, fontsize = fontlegend)
			plt.savefig(FilePlot, dpi = 200, format = 'eps')

			#call(["okular", FilePlot])
			call(["open", FilePlot])
			#plt.show()

def	FigureEnergyPIGS(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, TypePlot, beadsRef):
	FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef)

	FilePlotEnergy      = FilePlotName.SaveEnergy+".eps"
	outfileEnergy       = FilePlotEnergy
	print(outfileEnergy)
	call(["rm", FilePlotEnergy])

	FileToBePlotEnergy  = FilePlotName.SaveEnergy+".txt"
	FileToBePlotDIAG    = FilePlotName.SaveEnergyDIAG+".txt"
	print(FileToBePlotEnergy)
	print(FileToBePlotDIAG)
	if (numbmolecules <= 2):
		FileToBePlotMM = FilePlotName.SaveEnergyMM+".txt"
		beads2, energy2                	       = genfromtxt(FileToBePlotMM,unpack=True, usecols=[0,1], skip_header=0, skip_footer=0)
		print(energy2)
			
	if (numbmolecules <= 4):
		energy3                                = loadtxt(FileToBePlotDIAG,unpack=True, usecols=[1])
	BConstant        = support.GetBconst(molecule_rot)  # in wavenumber
	Units          	 = support.GetUnitConverter()
	BConstantK       = BConstant*Units.CMRECIP2KL
	energy2 = energy2*BConstantK
	energy3 = energy3*BConstantK
					
	if (variableName == "tau"):
		var2 = parameter/(beads2-1.0)
	
	if (variableName == "beta"):
		var2 = parameter*(beads2-1.0)
	
	fig        = plt.figure(figsize=(8, 4))
	plt.grid(True)
	font       = 20
	fontlegend = font/2.0

	valTau, valRotEnergy, valPotEnergy, valTotalEnergy, errorRotEnergy, errorPotEnergy, errorTotalEnergy = genfromtxt(FileToBePlotEnergy, unpack=True, usecols=[1, 3, 4, 5, 7, 8, 9])

	'''
	plt.subplot(3, 1, 1)
	YLabel = "Rotational"
	PlotEnergyPIGS(font, valTau,valRotEnergy,errorRotEnergy,variableName,YLabel)
	ymin, ymax = plt.ylim()
	xmin, xmax = plt.xlim()
	#PlotLabel(Text1,Text2,font,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
	'''

	plt.subplot(1, 2, 1)
	YLabel = "Total"
	PlotEnergyPIGS(font, valTau,valTotalEnergy,errorTotalEnergy,variableName,YLabel,var2,energy2,energy3)
	plt.legend(bbox_to_anchor=(0.33, 0.55), loc=2, borderaxespad=0., shadow=True, fontsize = fontlegend)

	plt.subplot(1, 2, 2)
	YLabel = "Potential"
	PlotEnergyPIGS(font, valTau,valPotEnergy,errorPotEnergy,variableName,YLabel,var2,energy2,energy3)
	plt.legend(bbox_to_anchor=(0.33, 0.55), loc=2, borderaxespad=0., shadow=True, fontsize = fontlegend)

	plt.subplots_adjust(top=0.95, bottom=0.20, left=0.15, right=0.95, hspace=0.0, wspace=0.4)
	plt.savefig(outfileEnergy, dpi = 200, format = 'eps')

	call(["open", outfileEnergy])
	#plt.show()

def PlotLabel(Text1, Text2, font,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment):
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy

	if Text1:
		plt.text(textpositionx-(4.45*deltax), textpositiony+(5.0*deltay), Text1, fontsize=font)
		#plt.text(textpositionx-(4.70*deltax), textpositiony+(5.5*deltay), Text1, fontsize=font)

	if Text2:
		plt.text(textpositionx+(0.0*deltax), textpositiony+0*deltay, r'$g = '+Text2+'$', fontsize=font)
		#plt.text(textpositionx+(4.3*deltax), textpositiony+5*deltay, r'$g$ = '+Text2, fontsize=font)

	'''
	if (variableName == "beta" and Text2):
		plt.text(textpositionx, textpositiony+0*deltay, 'Parameters:', fontsize=font)
		plt.text(textpositionx, textpositiony-1*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=font)
		plt.text(textpositionx, textpositiony-2*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=font)
		plt.text(textpositionx, textpositiony-3*deltay, r'$\mathrm{\mu} =$'+str(dipolemoment)+"Debye", fontsize=font)
		plt.text(textpositionx, textpositiony-4*deltay, r'$\mathrm{\tau} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=font)

	if (variableName == "tau" and Text2):
		#plt.text(textpositionx, textpositiony+5*deltay, 'Parameters:', fontsize=font)
		#plt.text(textpositionx, textpositiony+4*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=font)
		#plt.text(textpositionx, textpositiony+3*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=font)
		#plt.text(textpositionx, textpositiony+1*deltay, r'$\mathrm{\beta} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=font)
	'''
	

def PlotEnergyPIGS(font,val1,val2,val3,variableName,YLabel,var2,energy2,energy3):
	#plt.grid(True)
	if (variableName == "tau"):
		plt.xlabel(r'$\tau \ \  (\mathrm{K^{-1}})$', fontsize = font)
	if (variableName == "beta"):
		plt.xlabel(r'$\beta \ \  (\mathrm{K^{-1}})$', fontsize = font)
		
	if (YLabel == "Total"):
		plt.ylabel(r'$\mathrm{E_{0}}$ (K)', fontsize = font)
		if (variableName == "tau"):
			plotfitting(val1, val2, val3)
			plt.errorbar(val1, val2, yerr = val3, color = 'blue', ls = '-', label = 'PIGS', linewidth=2)
			plt.axhline(y=energy3, color='black', lw = 3.0, linestyle='--', label = 'Diagonalization')
			plt.plot(var2, energy2, color = 'green', ls = '--', label = 'Matrix Multiplication', marker = "o", linewidth=2)
		if (variableName == "beta"):
			plt.errorbar(val1, val2, yerr = val3, color = 'b', ls = '-', label = 'PIGS', linewidth=1)
			plt.xlim(0,0.201)
			plt.plot(var2, energy2, color = 'green', ls = '--', label = 'Matrix Multiplication', marker = "o", linewidth=2)
	if (YLabel == "Potential"):
		plt.ylabel(r'$\mathrm{V_{0}}$ (K)', fontsize = font)
		plt.xlim(0,0.201)
		plt.errorbar(val1, val2, yerr = val3, color = 'b', ls = '-', label = 'PIGS', linewidth=2)

	if (YLabel == "Rotational"):
		plt.ylabel(r'$\mathrm{K_{0}^{Rot}}$ (K)', fontsize = font)
		plt.xlim(0,0.201)
		plt.errorbar(val1, val2, yerr = val3, color = 'b', ls = '-', label = 'PIGS', linewidth=1)

	plt.legend(loc='upper right', shadow=True, fontsize = font)
	if (variableName == "tau"):
		plt.xlim(0.000,0.00505)
		x = [0.000, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005]
		labels = ['0.000', '', '0.001', '', '0.002', ' ', '0.003', '', '0.004', ' ', '0.005']
		plt.xticks(x, labels, rotation='horizontal')


def fitFunc(var, a, b):
	return a + b*var*var

def plotfitting(val1, val2, val3):
	markersize_fig = 10
	xdata = np.linspace(0, 0.006, 10000)
	fitParams, fitCovariances = curve_fit(fitFunc, val1, val2, sigma = val3)
	plt.plot(xdata, fitFunc(xdata, fitParams[0], fitParams[1]), linestyle = '-', color = 'r', label = 'Fit', lw = 2)

def FigureChemicalPotentialPIGS(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef):
	aa=[]
	bb=[]
	cc=[]
	dd=[]
	ee=[]
	for i in range(2,12):
		numbmolecules = i
		aa.append(numbmolecules)
		nparticle = np.array(aa)
		FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef)
		FileToBePlot   	  = FilePlotName.SaveEnergy+".txt"
		if os.path.isfile(FileToBePlot):
			col_beads, col_tau, col_tot, err_col_tot = loadtxt(FileToBePlot,unpack=True, usecols=[0,1,5,9])
			if (col_beads.size == 1):
				bb.append(col_beads)
				cc.append(col_tau)
				dd.append(col_tot)
				ee.append(err_col_tot)
			else:
				for j in range(0,col_beads.size):
					if (col_beads[j] == 61):
						bb.append(col_beads[j])
						cc.append(col_tau[j])
						dd.append(col_tot[j])
						ee.append(err_col_tot[j])

	nbeads = np.array(bb)
	ntau   = np.array(cc)
	ntot   = np.array(dd)
	nerrtot= np.array(ee)
			
	FilePlot = FilePlotName.SaveChemPot+".eps"
	outfile  = FilePlot

	fig = plt.figure(figsize=(6, 4))

	TypePlot1 = 3
	font=20
	fontlegend = font/2
	#plt.grid(True)
	plt.xlabel('N')

	if (TypePlot1 == 1):
		plt.ylabel(r'$\mathrm{\langle E_{0} \rangle (K)}$', fontsize = font)
	if (TypePlot1 == 2):
		plt.ylabel(r'$\mathrm{\langle E_{0} \rangle/N \ (K)}$', fontsize = font)
	if (TypePlot1 == 3):
		plt.ylabel(r'$\mathrm{\mu \ (K)}$', fontsize = font)

	plt.xlabel(r'$\mathrm{N}$', fontsize = font)

	if (TypePlot1 == 1):
		NumbRotors1  = nparticle
		TotalEnergy1 = ntot
		Error1       = nerrtot

	if (TypePlot1 == 2):
		NumbRotors1  = nparticle
		TotalEnergy1 = ntot/nparticle
		Error1       = nerrtot/nparticle

	if (TypePlot1 == 3):
		mu1 = []
		errormu1 = []
		num1 = []
		print(ntot)
		print(nparticle)
		for i in range(len(nbeads)):
			ii = i-1
			if ii <0:
				mu1.append(ntot[i])
				num1.append(nparticle[i])
			else:		
				mu1.append(ntot[i] - ntot[ii])
				num1.append(nparticle[i])
			errormu1.append(sqrt(nerrtot[i]*nerrtot[i]+nerrtot[ii]*nerrtot[ii]))

		NumbRotors1  = num1
		TotalEnergy1 = mu1
		Error1       = errormu1

	srcfile2         = "ResultsOfPIGS/chemical_potential_unscreened.dat"
	data2            = loadtxt(srcfile2,unpack=True, usecols=[0,1])
	xdata, ydata = data2
	for i in range(10):
		if i <6:
			TotalEnergy1[i]=ydata[i]-Error1[i]
		else:
			TotalEnergy1[i]=ydata[6]-Error1[i]
	
	#Manupulation#
	print(NumbRotors1)
	print(TotalEnergy1)
	print(xdata)
	print(ydata)
	#
	plt.plot(xdata, ydata, linestyle = '--', color = 'r', label = 'Diagonalization', marker = 'v', lw = 2)
	plt.errorbar(NumbRotors1, TotalEnergy1, yerr=Error1, color = 'b', ls = '-', label = 'PIGS', linewidth=2)
	plt.subplots_adjust(top=0.95, bottom=0.15, left=0.20, right=0.98, hspace=0.6, wspace=1.0)
	plt.legend(bbox_to_anchor=(0.40, 0.98), loc=2, borderaxespad=0.)
	plt.savefig(outfile, dpi = 200, format = 'eps')
	call(["open", outfile])
	#call(["okular", outfile])

def	FigureAngleDistribution(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, TypePlot, beadsRef):
	font       = 20
	fontlegend = font/2.0
	preskip    = 0
	postskip   = 0

	if (((TypePlot == "RFACTOR") or (TypePlot == "GFACTOR")) and variableName == "tau"):
		fig        = plt.figure(figsize=(8, 4))
		plt.grid(True)

		var = beadsRef
		FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, var)
		FilePlotCorr = FilePlotName.SaveCorr+".eps"
		outfile      = FilePlotCorr
		print(outfile)
		exit(0)
		call(["rm", FilePlotCorr])
#
		colorList  = ['red', 'green', 'blue']
		lsList     = ['-', '--', '-.']
		markerList = ['o', '^', 'v']
		labelList  = ['2 Rotors','4 Rotors', '6 Rotors']
		nn = [2,4,6]
		iLabel = 0
		for numbmolecules in nn:
			if (numbmolecules == 2):
				DList  = [1.0+0.5*i for i in range(15)]
			if (numbmolecules == 4):
				DList  = [1.0+0.5*i for i in range(11)]
			if (numbmolecules == 6):
				DList  = [1.0+0.5*i for i in range(6)]
		
			RFactorPlot  = np.zeros(len(DList))
			CorrPlot     = np.zeros(len(DList))
			err_CorrPlot = np.zeros(len(DList))

			iii = 0
			for dipolemoment in DList:
				RFactorList = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
				RFactor = RFactorList[0]
				FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef)
				FileToBePlotCorr = FilePlotName.SaveCorr+".txt"
				beads1, var1, Corr, err_Corr = genfromtxt(FileToBePlotCorr,unpack=True, usecols=[0, 1, 11, 21], skip_header=preskip, skip_footer=postskip)
				if (TypePlot == "GFACTOR"):
					RFactorPlot[iii] = 1.0/(RFactor*RFactor*RFactor)
				if (TypePlot == "RFACTOR"):
					RFactorPlot[iii] = RFactor
	
				ii = 0
				for i in beads1:
					indexi =int(i+0.5)
					beads = indexi
					if beads == beadsRef:
						CorrPlot[iii]     = Corr[ii]
						err_CorrPlot[iii] = err_Corr[ii]
					ii += 1
				iii += 1
		
			print("Corr")
			print(CorrPlot)
#
			plt.errorbar(RFactorPlot, CorrPlot, yerr=err_CorrPlot, color = colorList[iLabel], ls = lsList[iLabel], linewidth=1,  marker = markerList[iLabel], markersize = 8, label = labelList[iLabel])
			iLabel += 1

			ymin, ymax = plt.ylim()
			midpointy = 0.5*(ymax-ymin)
			deltay = midpointy*0.15
			xmin, xmax = plt.xlim()
			midpointx = 0.5*(xmax-xmin)
			deltax = midpointx*0.15
			textpositionx = xmin+midpointx-0.25*midpointx
			textpositiony = ymin+midpointy

			if (TypePlot == "RFACTOR"):
				plt.xlabel(r'$\mathrm{R}$', fontsize = font)
			if (TypePlot == "GFACTOR"):
				plt.xlabel(r'$\mathrm{g}$', fontsize = font)
			plt.ylabel(r'$\mathrm{z^{2}}$', fontsize = font)
			ymin, ymax = plt.ylim()
			if ymin < 0.0:
				plt.ylim(0.0,ymax)
			xmin, xmax = plt.xlim()
			'''
			Text1 = ""
			Text2 = ""
			if Text1:
				PlotLabel(Text1, Text2,fontlegend,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
			'''

		plt.subplots_adjust(top=0.95, bottom=0.15, left=0.09, right=0.98, hspace=0.0, wspace=0.4)
		plt.legend(bbox_to_anchor=(0.78, 0.75), loc=2, borderaxespad=1., shadow=True, fontsize = fontlegend)
		plt.savefig(outfile, dpi = 200, format = 'eps')

		call(["open", outfile])
		#call(["okular", outfile])
		#plt.show()

def	FigureAngleDistributionGfactor(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, TypePlot, beadsRef):
	font       = 28
	fontlegend = font/2.0
	preskip    = 0
	postskip   = 0

	if (((TypePlot == "RFACTOR") or (TypePlot == "GFACTOR")) and variableName == "tau"):
		fig        = plt.figure(figsize=(8, 6))
		plt.grid(True)

		var = beadsRef
		gFact = -1.0
		FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gFact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, var)
		FilePlotCorr = FilePlotName.SaveCorrGFAC+".eps"
		outfile      = FilePlotCorr
		print("--------------------------------------------------")
		print("Name of the path and the Figure is given below - ")
		print(outfile)
		print("--------------------------------------------------")
		call(["rm", FilePlotCorr])
#
		colorList  = ['red', 'green', 'blue', 'magenta']
		lsList     = ['-', '--', '-.', '-']
		markerList = ['o', '^', 'v','s']
		labelList  = ['N = 24', 'N = 48','N = 64']
		nn = [24,48]
		iLabel = 0
		gFactorList  = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]
		DList=np.zeros(len(gFactorList))
		ig = 0
		for gFactor in gFactorList:
			DipoleMoment = support.GetDipoleMomentFromGFactor(molecule, Rpt, gFactor)
			output  = '{:1.4f}'.format(DipoleMoment)
			DList[ig] = output
			ig = ig+1

		#print(DList)
		for numbmolecules in nn:
			gFactorPlot  = np.zeros(len(DList))
			CorrPlot     = np.zeros(len(DList))
			err_CorrPlot = np.zeros(len(DList))

			iii = 0
			for dipolemoment in DList:
				gFact = -1.0
				FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gFact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef)
				FileToBePlotCorr = FilePlotName.SaveCorr+".txt"
				if os.path.isfile(FileToBePlotCorr):
					beads1, var1, Corr, err_Corr = genfromtxt(FileToBePlotCorr,unpack=True, usecols=[0, 1, 8, 15], skip_header=preskip, skip_footer=postskip)
					FactorList = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
					if (TypePlot == "GFACTOR"):
						gFactorPlot[iii] = FactorList[1]
					if (TypePlot == "RFACTOR"):
						rFactorPlot[iii] = FactorList[0]

					beadsRef1 = beadsRef
					#if ((beadsRef not in Corr) and (np.isscalar(Corr) == False)):
					#	beadsRef1 = beads1[-1]
					#	print(beadsRef1)
	
					if (np.isscalar(Corr) == True):
						CorrPlot[iii]     = Corr
						err_CorrPlot[iii] = err_Corr
					else:
						ii = 0
						for i in beads1:
							indexi =int(i+0.5)
							beads = indexi
							if beads == beadsRef1:
								CorrPlot[iii]     = Corr[ii]
								err_CorrPlot[iii] = err_Corr[ii]
							ii += 1
					iii += 1
					beadsRef1 = beadsRef
		
				#print("Corr")
				#print(CorrPlot)
#
			plt.errorbar(gFactorPlot, CorrPlot, yerr=err_CorrPlot, color = colorList[iLabel], ls = lsList[iLabel], linewidth=1,  marker = markerList[iLabel], markersize = 8, label = labelList[iLabel])
			iLabel += 1

			plt.xlim(0.4,1.6)
			plt.ylim(0.0,0.8)
			ymin, ymax = plt.ylim()
			midpointy = 0.5*(ymax-ymin)
			deltay = midpointy*0.15
			xmin, xmax = plt.xlim()
			midpointx = 0.5*(xmax-xmin)
			deltax = midpointx*0.15
			textpositionx = xmin+midpointx-0.25*midpointx
			textpositiony = ymin+midpointy

			if (TypePlot == "RFACTOR"):
				plt.xlabel(r'$\mathrm{R}$', fontsize = font)
			if (TypePlot == "GFACTOR"):
				plt.xlabel(r'${g}$', fontsize = font)
			plt.ylabel(r'${\phi_{\mathrm{abs}}}$', fontsize = font)
			ymin, ymax = plt.ylim()
			if ymin < 0.0:
				plt.ylim(0.0,ymax)
			xmin, xmax = plt.xlim()

		plt.xticks(fontsize=font, rotation=0)
		plt.yticks(fontsize=font, rotation=0)
		
		plt.subplots_adjust(top=0.96, bottom=0.15, left=0.14, right=0.95, hspace=0.0, wspace=0.4)
		plt.legend(bbox_to_anchor=(0.75, 0.55), loc=2, borderaxespad=1., shadow=True, fontsize = fontlegend)
		plt.savefig(outfile, dpi = 200, format = 'eps')

		#call(["open", outfile])
		plt.show()

def	FigureEntropyRT(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, parameterName, parameter, numbblocks, numbpass, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, TypePlot, beadsRef):
	font       = 28
	fontlegend = font/2.0

	font = 28
	fontlegend = font/2
	fig        = plt.figure(figsize=(8, 6))
	plt.grid(True)

	var = beadsRef
	gfact = 0.0
	dipolemoment = -1.0
	numbmolecules = 2
	particleA = 1
	FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, var)
	FilePlotEntropy  = FilePlotName.SaveEntropyGFAC+".eps"
	outfileEntropy       = FilePlotEntropy
	call(["rm", FilePlotEntropy])
#
	plt.axhline(y=log(2.0), color='blue', lw = 2.0, linestyle='-', label = 'ln(2)')
#
	nn = [2]
	for numbmolecules in nn:
		
		if (numbmolecules == 2):
			gFactorList  = [0.5+0.25*i for i in range(21)]  # 2 HF
			gFactorList += [5.75+0.25*i for i in range(14)] # 2 HF
		if (numbmolecules == 4):
			gFactorList  = [0.5+0.25*i for i in range(15)]   # 4 HF 

		gFactorPlot      = gFactorList
		entropy1Plot     = np.zeros(len(gFactorList))
		purity1Plot      = np.zeros(len(gFactorList))
		err_entropy1Plot = np.zeros(len(gFactorList))
		err_purity1Plot  = np.zeros(len(gFactorList))

		iii = 0
		for gFact in gFactorList:
			FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gFact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef)
			FileToBePlotEntropy = FilePlotName.SaveEntropyRT+".txt"
			#FileToBePlotDIAG    = FilePlotName.SaveEntropyDIAG+".txt"

			beads1, var1, purity1, entropy1 = genfromtxt(FileToBePlotEntropy,unpack=True, usecols=[0, 1, 2,3], skip_header=0, skip_footer=0)
					
			if (np.isscalar(entropy1) == True):
				entropy1Plot[iii]     = entropy1
				purity1Plot[iii]      = purity1
				err_entropy1Plot[iii] = 0.0
			else:
				ii = 0
				for i in beads1:
					indexi =int(i+0.5)
					beads = indexi
					if beads == beadsRef:
						entropy1Plot[iii]     = entropy1[ii]
						purity1Plot[iii]      = purity1[ii]
						err_entropy1Plot[iii] = 0.0
	
					ii += 1
				iii += 1
		
		print("S2:  PIGS "+str(numbmolecules))
		print(entropy1Plot)
#
		plt.plot(gFactorPlot, entropy1Plot, color = 'black', ls = 'None', linewidth=1, marker = "v", markersize = 10, label = 'ED: N=2')

	ymin, ymax = plt.ylim()
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy
	plt.xticks(fontsize=font, rotation=0)
	plt.yticks(fontsize=font, rotation=0)

	if (TypePlot == "RFACTOR"):
		plt.xlabel(r'$R$', fontsize = font, labelpad=0)
	if (TypePlot == "GFACTOR"):
		plt.xlabel(r'$g$', fontsize = font, labelpad=-3)

	ymin, ymax = plt.ylim()
	xmin, xmax = plt.xlim()
	plt.ylabel(r'$S_{2}$', fontsize = font)
	if ymin < 0.0:
		plt.ylim(0.0,1.001)
	if xmin < 0.0:
		plt.xlim(0.0,xmax)
	plt.xticks(np.arange(0, 9, step=2))
	plt.yticks(np.arange(0, 0.9, step=0.2))
	Text1 = "(a)"
	Text2 = ""
	if Text1:
		PlotLabel(Text1, Text2,font,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
		
	plt.subplots_adjust(top=0.97, bottom=0.14, left=0.14, right=0.98, hspace=0.0, wspace=0.0)
	plt.legend(bbox_to_anchor=(0.68, 0.75), loc=2, borderaxespad=1., shadow=True, fontsize = fontlegend)
	plt.savefig(outfileEntropy, dpi = 200, format = 'eps')

	#call(["open", outfileEntropy])
	plt.show()
