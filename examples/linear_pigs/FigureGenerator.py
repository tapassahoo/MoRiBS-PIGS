#!/usr/bin/env python
import numpy as np
from numpy import *
import numpy as np
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
from pylab import *

class Specification:
	def __init__(self):
		self.color1               = itertools.cycle(('b', 'r', 'g', 'm'))
		self.marker1              = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))
		self.line1                = itertools.cycle(('--','-','-.',':'))
		self.color2               = itertools.cycle(('b', 'r', 'g', 'm'))
		self.marker2              = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))
		self.line2                = itertools.cycle(('--','-','-.',':'))

def	FigureENT(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, TypePlot, beadsRef):

	fig        = plt.figure(figsize=(8, 4), dpi=200)
	plt.grid(True)

	font       = 40
	fontlegend = font/2.0
	preskip    = 0
	postskip   = 0

	if (((TypePlot == "RFACTOR") or (TypePlot == "GFACTOR")) and variableName == "tau"):
		var = beadsRef
		FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, var)
		if (TypePlot == "RFACTOR"):
			FilePlotEntropy  = FilePlotName.SaveEntropyRFAC+".pdf"
			FilePlotEnergy   = FilePlotName.SaveEnergyRFAC+".pdf"
		if (TypePlot == "GFACTOR"):
			FilePlotEntropy  = FilePlotName.SaveEntropyGFAC+".pdf"
			FilePlotEnergy   = FilePlotName.SaveEnergyGFAC+".pdf"
		outfileEntropy       = FilePlotEntropy
		outfileEnergy        = FilePlotEnergy
		print(outfileEntropy)
		print(outfileEnergy)
		call(["rm", FilePlotEntropy])
		call(["rm", FilePlotEnergy])
#
		DMRGdata                        	       = genfromtxt("ResultsOfPIGSENT/large_g_S2.txt")
		DMRGdatagFac = DMRGdata[:,0]
#
		nn = [2,4,6]
		iDMRGdata = 1
		for numbmolecules in nn:
			DMRGdataPlot = DMRGdata[:,iDMRGdata]
			iDMRGdata += 1

			particleA = numbmolecules/2
			if (numbmolecules == 2):
				DList  = [1.0+0.5*i for i in range(7)]
				DList  += [4.5, 5.0, 6.0]
				numbblocks = 50000
			if (numbmolecules == 4):
				DList  = [1.0+0.5*i for i in range(7)]
				numbblocks = 20000
			if (numbmolecules == 6):
				DList  = [1.0+0.5*i for i in range(6)]
				numbblocks = 20000
		
			RFactorPlot      = np.zeros(len(DList))
			entropy1Plot     = np.zeros(len(DList))
			purity1Plot      = np.zeros(len(DList))
			err_entropy1Plot = np.zeros(len(DList))
			err_purity1Plot  = np.zeros(len(DList))
			entropy2Plot     = np.zeros(len(DList))
			entropy3Plot     = np.zeros(len(DList))

			energy1Plot      = np.zeros(len(DList))
			err_energy1Plot  = np.zeros(len(DList))
			energy2Plot      = np.zeros(len(DList))
			err_energy2Plot  = np.zeros(len(DList))
			energy3Plot      = np.zeros(len(DList))
			err_energy3Plot  = np.zeros(len(DList))
			BConstant        = support.GetBconst(molecule_rot)  # in wavenumber
			Units          	 = support.GetUnitConverter()
			BConstantK       = BConstant*Units.CMRECIP2KL

			iii = 0
			for dipolemoment in DList:
				RFactor = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
				if (dipolemoment > 4.0):
					numbblocks = 20000
				FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef)
				FileToBePlotEntropy = FilePlotName.SaveEntropy+".txt"
				FileToBePlotEnergy  = FilePlotName.SaveEnergy+".txt"
				FileToBePlotDIAG    = FilePlotName.SaveEntropyDIAG+".txt"
				if (numbmolecules <= 2):
					FileToBePlotMM    = FilePlotName.SaveEntropyMM+".txt"
					beads2, energy2, entropy2                      	       = genfromtxt(FileToBePlotMM,unpack=True, usecols=[0,1,2], skip_header=preskip, skip_footer=postskip)
			
				beads1, var1, purity1, entropy1, err_purity1, err_entropy1 = genfromtxt(FileToBePlotEntropy,unpack=True, usecols=[0, 1, 4, 5, 8, 9], skip_header=preskip, skip_footer=postskip)
				beads1, var1, energy1, err_energy1                         = genfromtxt(FileToBePlotEnergy,unpack=True, usecols=[0, 1, 5, 9], skip_header=preskip, skip_footer=postskip)
				if (numbmolecules <= 6):
					RFactor, energy3, entropy3                             = loadtxt(FileToBePlotDIAG,unpack=True, usecols=[0,1,2])
					
				print(beads1)
				if (TypePlot == "GFACTOR"):
					RFactorPlot[iii] = 1.0/(RFactor*RFactor*RFactor)
				if (TypePlot == "RFACTOR"):
					RFactorPlot[iii] = RFactor
				entropy3Plot[iii] = entropy3
				energy3Plot[iii]  = energy3*BConstantK
	
				if (numbmolecules <=2 ):
					ii = 0
					for i in beads2:
						indexi =int(i+0.5)
						beads = indexi
						if beads == beadsRef:
							entropy2Plot[iii] = entropy2[ii]
							energy2Plot[iii]  = energy2[ii]*BConstantK
						ii += 1
				if ((numbmolecules == 4) and (dipolemoment ==3.5)):
					beadsRef = 41
				elif ((numbmolecules == 6) and (dipolemoment == 3.5)):
					beadsRef = 21
				else:
					beadsRef = 61
	
				if ((numbmolecules == 8) and (dipolemoment > 2.5)):
					entropy1Plot[iii]     = entropy1
					err_entropy1Plot[iii] = err_entropy1
					purity1Plot[iii]      = purity1
					err_purity1Plot[iii]  = err_purity1
	
					energy1Plot[iii]     = energy1
					err_energy1Plot[iii] = err_energy1
				else:
					ii = 0
					for i in beads1:
						indexi =int(i+0.5)
						beads = indexi
						if beads == beadsRef:
							entropy1Plot[iii]     = entropy1[ii]
							err_entropy1Plot[iii] = err_entropy1[ii]
							purity1Plot[iii]      = purity1[ii]
							err_purity1Plot[iii]  = err_purity1[ii]
		
							energy1Plot[iii]     = energy1[ii]
							err_energy1Plot[iii] = err_energy1[ii]
						ii += 1
				iii += 1
		
			'''
			print("S2:  PIGS ")
			print(entropy1Plot)
			print("S2:  MM ")
			print(entropy2Plot)
			print("S2:  DIAG ")
			print(entropy3Plot)
#
			print("")
			print("")
			print("")
			print("E0:  PIGS ")
			print(energy1Plot)
			print("E0:  MM ")
			print(energy2Plot)
			print("E0:  DIAG ")
			print(energy3Plot)
#
			fig       = plt.figure(figsize=(4, 4), dpi=200)
			plt.grid(True)
			#rc('axes', linewidth=1)
	
			plotEntropyENT1(numbmolecules,RFactorPlot, energy1Plot, err_energy1Plot, variableName, RFactorPlot, energy2Plot, energy3Plot, font, TypePlot)
			plt.ylabel(r'$\mathrm{\langle E_{0} \rangle \ \  (Kelvin^{-1}})$', fontsize = font)
			ymin, ymax = plt.ylim()
			xmin, xmax = plt.xlim()
			Text1 = ""
			Text2 = ""
			if Text1:
				PlotLabel(Text1, Text2,fontlegend,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
			plt.subplots_adjust(top=0.95, bottom=0.20, left=0.1, right=0.95, hspace=0.0, wspace=0.4)
			plt.legend(bbox_to_anchor=(0.73, 0.95), loc=2, borderaxespad=0., shadow=True, fontsize = fontlegend)
			plt.savefig(outfileEnergy, dpi = 200, format = 'pdf')
			call(["open", outfileEnergy])
			'''
#
			plotEntropyENT1(numbmolecules,RFactorPlot, entropy1Plot, err_entropy1Plot, variableName, RFactorPlot, entropy2Plot, entropy3Plot,DMRGdatagFac,DMRGdataPlot, font, TypePlot)
			plt.ylabel(r'$\mathrm{S_{2}}$', fontsize = font)
			ymin, ymax = plt.ylim()
			if ymin < 0.0:
				plt.ylim(0.0,ymax)
			xmin, xmax = plt.xlim()
			Text1 = ""
			Text2 = ""
			if Text1:
				PlotLabel(Text1, Text2,fontlegend,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
			
		
			'''
			plt.subplot(1, 2, 2)
			plt.grid(True)
			purity2Plot = np.exp(-entropy2Plot)
			purity3Plot = np.exp(-entropy3Plot)
			plotEntropyENT1(numbmolecules,RFactorPlot, purity1Plot, err_purity1Plot, variableName, RFactorPlot, purity2Plot, purity3Plot, font, TypePlot)
			plt.ylabel(r'$\mathrm{\langle P_{A} \rangle}$', fontsize = font)
			ymin, ymax = plt.ylim()
			xmin, xmax = plt.xlim()
			Text1 = "(b)"
			Text2 = ""
			if Text1:
				PlotLabel(Text1, Text2,fontlegend,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
			plt.legend(bbox_to_anchor=(0.48, 0.75), loc=2, borderaxespad=0., shadow=True, fontsize = fontlegend)
			plt.subplots_adjust(top=0.95, bottom=0.20, left=0.1, right=0.95, hspace=0.0, wspace=0.4)
			'''

		plt.subplots_adjust(top=0.95, bottom=0.20, left=0.1, right=0.95, hspace=0.0, wspace=0.4)
		plt.legend(bbox_to_anchor=(0.78, 0.75), loc=2, borderaxespad=1., shadow=True, fontsize = fontlegend)
		plt.savefig(outfileEntropy, dpi = 200, format = 'pdf')

		call(["open", outfileEntropy])
		#call(["okular", outfileEntropy])
		#plt.show()

	if (TypePlot == "S2"):
		iFrame = 1
		iFigLabel = 0
		FigureLabel  = ["(a)","(b)"]
		variableList = ["beta","tau"]
		for variableName in variableList:
			if variableName == "beta":
				parameterName = "tau"
				parameter     = 0.005
			if variableName == "tau":
				parameterName = "beta"
				parameter     = 0.2
##
			fig        = plt.figure(figsize=(8, 6), dpi=200)
			plt.grid(True)
##	
			FilePlotName      = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip1, postskip1, extra_file_name, src_dir, particleA, beadsRef)
			FileToBePlot   	  = FilePlotName.SaveEntropy+".txt"
			FileToBePlotMM    = FilePlotName.SaveEntropyMM+".txt"
			FileToBePlotDIAG  = FilePlotName.SaveEntropyDIAG+".txt"
			print(FileToBePlotMM)
			print(FileToBePlotDIAG)
##
			FilePlot     = FilePlotName.SaveEntropy+".pdf"
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
	
			#plt.subplot(1, 2, iFrame)
			iFrame = iFrame + 1
			plotEntropyENT(var1, entropy1, err_entropy1, variableName, var2, entropy2, entropy3, font,fontlegend)
			plt.ylabel(r'$\mathrm{S_{2}}$', fontsize = font)
			ymin, ymax = plt.ylim()
			xmin, xmax = plt.xlim()
			Text1 = FigureLabel[iFigLabel]
			iFigLabel = iFigLabel +1
			RGFactor = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
			gFactor = RGFactor[1]
			arg     = "%3.2f" % gFactor
			Text2 = arg
			PlotLabel(Text1, Text2,fontlegend,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
			plt.legend(bbox_to_anchor=(0.65, 0.75), loc=2, borderaxespad=0., shadow=True, fontsize = fontlegend)
			plt.subplots_adjust(top=0.95, bottom=0.17, left=0.20, right=0.95, hspace=0.0, wspace=0.2)
	
			plt.savefig(outfile, dpi = 200, format = 'pdf')
	
			call(["open", outfile])
			#call(["okular", outfile])
			#plt.show()


def plotEntropyENT1(numbmolecules,var, val, err_val, variableName, var1, val1, val2, DMRGdatagFac, DMRGdataPlot, font, TypePlot):
	#plt.xlim(0,0.201)
	x = Specification()
	if (numbmolecules == 2):
		plt.errorbar(var, val, yerr=err_val, color = "red", ls = '-', linewidth=1,  marker = "o", markersize = 10, label = 'PIGS: 2 HF')
		plt.plot(var, val2, color = 'black', ls = 'None', linewidth=1, marker = "o", markersize = 8, label = 'DIAG: 2 HF')
		#plt.plot(var, val2, color = 'black', ls = '-', linewidth=3, marker = "o", markersize = 12, label = 'DIAG: 2 HF')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "o", markersize = 10, label = 'DMRG: 2 HF')
		#plt.plot(var, val1, linestyle = 'None', color = 'mediumpurple', marker = "o", markersize = 10, label = 'MM: 2 HF')

	if (numbmolecules == 4):
		plt.errorbar(var, val, yerr=err_val, color = "red", ls = '-', linewidth=1,  marker = "s", markersize = 10, label = 'PIGS: 4 HF')
		plt.plot(var, val2, color = 'black', ls = 'None', linewidth=1, marker = "s", markersize = 8, label = 'DIAG: 4 HF')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "s", markersize = 8, label = 'DMRG: 4 HF')

	if (numbmolecules == 6):
		plt.errorbar(var, val, yerr=err_val, color = "red", ls = '-', linewidth=1,  marker = "v", markersize = 10, label = 'PIGS: 6 HF')
		plt.plot(var, val2, color = 'black', ls = 'None', linewidth=1, marker = "v", markersize = 8, label = 'DIAG: 6 HF')
		#plt.plot(DMRGdatagFac, DMRGdataPlot, linestyle = 'None', color = 'darkcyan', marker = "8", markersize = 7, label = 'DMRG: 6 HF')

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

def plotEntropyENT(var, val, err_val, variableName, var1, val1, val2, font,fontlegend):
	if (variableName == "beta"):
		plt.xlim(0,0.201)
	plt.errorbar(var, val, yerr=err_val, color = 'red', ls = '-', linewidth=2,  marker = "o", markersize = 10, label = 'PIGS')
	if (val1[0] != 0.0):
		plt.plot(var1, val1, linestyle = '--', linewidth=2,color = 'mediumpurple', marker = "v", markersize = 10, label = 'MM')

	if (variableName == "tau"):
		plt.xlim(0,0.051)
		if (val2 != 0.0):
			plt.axhline(y=val2, color='black', lw = 2.0, linestyle='-', label = 'DIAG')

	ymin, ymax = plt.ylim()
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy
	print(textpositionx)
	print(textpositiony)
	plt.xticks(fontsize=fontlegend, rotation=0)
	plt.yticks(fontsize=fontlegend, rotation=0)


	if (variableName == "beta"):
		plt.xlabel(r'$\mathrm{\beta \ \  (Kelvin^{-1}})$', fontsize = font)
	if (variableName == "tau"):
		plt.xlabel(r'$\mathrm{\tau \ \  (Kelvin^{-1}})$', fontsize = font)

def FigureCorrelation(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef, RefPoint):

	font          = 56
	fontlegend    = font/2.0

	TextLabel = ["(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)"]
	DList  = [1.0+i for i in range(4)]
	iTextLabel = 0
	for dipolemoment in DList:
		RGFactor = support.GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
		gFactor = RGFactor[1]
		arg     = "%7.2f" % gFactor

		TypeCorrList = ["TotalCorr","ZCorr","XYCorr"]
		for TypeCorr in TypeCorrList:
			FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef)
			if TypeCorr    == "TotalCorr":
				FileToBePlot   	  = FilePlotName.SaveTotalCorr+".txt"
				FilePlot          = FilePlotName.SaveTotalCorr+".pdf"
			elif TypeCorr  == "ZCorr":
				FileToBePlot   	  = FilePlotName.SaveZCorr+".txt"
				FilePlot          = FilePlotName.SaveZCorr+".pdf"
			elif TypeCorr  == "XYCorr":
				FileToBePlot   	  = FilePlotName.SaveXYCorr+".txt"
				FilePlot          = FilePlotName.SaveXYCorr+".pdf"
			print(FilePlot)
			call(["rm", FilePlot])
			datacorr      = genfromtxt(FileToBePlot)
			print(datacorr[1,:])
			index = 0
			for i in range(len(datacorr)):
				if datacorr[i,0] == 41:
					findex = index
					break
				index = index+1
			
			FuncCorr      = np.zeros((numbmolecules,numbmolecules))
			ErrorFuncCorr = np.zeros((numbmolecules,numbmolecules))
			ii            = 0
			for i in range(numbmolecules-1):
				for j in range((i+1),numbmolecules):
			
					nc = 2+(2*ii)
					nec = nc+1 
					FuncCorr[i,j] = datacorr[1,nc]
					FuncCorr[j,i] = FuncCorr[i,j]
					ErrorFuncCorr[i,j] = datacorr[1,nec]
					ErrorFuncCorr[j,i] = ErrorFuncCorr[i,j]
					ii = ii+1
			val1 = np.arange(numbmolecules+1)
			val2 = np.zeros(numbmolecules)
			val3 = np.zeros(numbmolecules)
			for j in range(numbmolecules):
				if (j != RefPoint):
					val2[j] = FuncCorr[RefPoint:(RefPoint+1),j]
					val3[j] = ErrorFuncCorr[RefPoint:(RefPoint+1),j]
		
			fig = plt.figure(figsize=(8, 6), dpi=200)
			plt.errorbar(val1[RefPoint+2:], val2[RefPoint+1:], yerr=val3[RefPoint+1:], color = 'red', ls = '-', linewidth=3, marker = "o", markersize = 20)
			print(val1[RefPoint+2:])

			plt.grid(True)

			ymin, ymax = plt.ylim()
			midpointy = 0.5*(ymax-ymin)
			deltay = midpointy*0.15
			xmin, xmax = plt.xlim()
			plt.xlim(xmin-0.1,xmax+0.1)
			midpointx = 0.5*(xmax-xmin)
			deltax = midpointx*0.15
			textpositionx = xmin+midpointx-0.25*midpointx
			textpositiony = ymin+midpointy

			plt.xlabel(r'$\mathrm{j}^{th}$'+"rotor", fontsize = font)
	
			if (TypeCorr == "TotalCorr"):
				plt.ylabel(r'$\mathrm{C_{1j}}$', fontsize = font)
			if (TypeCorr == "XCorr"):
				plt.ylabel(r'$\mathrm{C^{X}_{1j}}$', fontsize = font)
			if (TypeCorr == "YCorr"):
					plt.ylabel(r'$\mathrm{C^{Y}_{1j}}$', fontsize = font)
			if (TypeCorr == "ZCorr"):
				plt.ylabel(r'$\mathrm{C^{Z}_{1j}}$', fontsize = font)
			if (TypeCorr == "XYCorr"):
				plt.ylabel(r'$\mathrm{C^{X,Y}_{1j}}$', fontsize = font)

			Text1 = TextLabel[iTextLabel]
			Text2 = str(arg)
			if Text1:
				PlotLabel(Text1, Text2,fontlegend,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment)
		
			iTextLabel += 1

			plt.xticks(fontsize=fontlegend, rotation=0)
			plt.yticks(fontsize=fontlegend, rotation=0)
			plt.subplots_adjust(top=0.95, bottom=0.25, left=0.35, right=0.95, hspace=0.6, wspace=1.0)
			plt.legend(bbox_to_anchor=(0.70, 0.70), loc=2, borderaxespad=0., shadow=True, fontsize = font)
			plt.savefig(FilePlot, dpi = 200, format = 'pdf')

			#call(["okular", FilePlot])
			call(["open", FilePlot])
			#plt.show()

def	FigureEnergyPIGS(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, TypePlot, beadsRef):
	FilePlotName = support.GetFileNamePlot(TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, preskip, postskip, extra_file_name, src_dir, particleA, beadsRef)

	FilePlotEnergy      = FilePlotName.SaveEnergy+".pdf"
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
	
	fig        = plt.figure(figsize=(8, 4), dpi=200)
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
	plt.savefig(outfileEnergy, dpi = 200, format = 'pdf')

	call(["open", outfileEnergy])
	#call(["okular", outfileEnergy])
	#plt.show()

def PlotLabel(Text1, Text2, font,xmin,xmax,ymin,ymax,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment):
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy

	if Text1:
		plt.text(textpositionx-(4.75*deltax), textpositiony+(5.5*deltay), Text1, fontsize=font)

	if Text2:
		plt.text(textpositionx-(1.0*deltax), textpositiony-1*deltay, r'$\mathrm{g} =$'+Text2, fontsize=font)

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
		plt.xlabel(r'$\mathrm{\tau \ \  (Kelvin^{-1}})$', fontsize = font)
	if (variableName == "beta"):
		plt.xlabel(r'$\mathrm{\beta \ \  (Kelvin^{-1}})$', fontsize = font)
		
	if (YLabel == "Total"):
		plt.ylabel(r'$\mathrm{E_{0}}$ (Kelvin)', fontsize = font)
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
		plt.ylabel(r'$\mathrm{V_{0}}$ (Kelvin)', fontsize = font)
		plt.xlim(0,0.201)
		plt.errorbar(val1, val2, yerr = val3, color = 'b', ls = '-', label = 'PIGS', linewidth=2)

	if (YLabel == "Rotational"):
		plt.ylabel(r'$\mathrm{K_{0}^{Rot}}$ (Kelvin)', fontsize = font)
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
			
	FilePlot = FilePlotName.SaveChemPot+".pdf"
	outfile  = FilePlot

	fig = plt.figure(figsize=(6, 4), dpi=200)

	TypePlot1 = 3
	font=20
	fontlegend = font/2
	#plt.grid(True)
	plt.xlabel('N')

	if (TypePlot1 == 1):
		plt.ylabel(r'$\mathrm{\langle E_{0} \rangle (Kelvin)}$', fontsize = font)
	if (TypePlot1 == 2):
		plt.ylabel(r'$\mathrm{\langle E_{0} \rangle/N \ (Kelvin)}$', fontsize = font)
	if (TypePlot1 == 3):
		plt.ylabel(r'$\mathrm{\mu \ (Kelvin)}$', fontsize = font)

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
	plt.savefig(outfile, dpi = 200, format = 'pdf')
	call(["open", outfile])
	#call(["okular", outfile])
