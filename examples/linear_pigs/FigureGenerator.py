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

class Specification:
	def __init__(self):
		self.color1               = itertools.cycle(('b', 'r', 'g', 'm'))
		self.marker1              = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))
		self.line1                = itertools.cycle(('--','-','-.',':'))
		self.color2               = itertools.cycle(('b', 'r', 'g', 'm'))
		self.marker2              = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))
		self.line2                = itertools.cycle(('--','-','-.',':'))

def	FigureENT(FileToBePlot,FilePlot,TypeCal,variableName,parameter,ExactValue,numbmolecules,molecule,Rpt,dipolemoment):
	outfile = FilePlot
	fig = plt.figure(figsize=(6, 4), dpi=200)

	font=18
	#plt.grid(True)
	#plt.xlim(0,0.25)

	#x      = Specification()

	data1            = loadtxt(FileToBePlot,unpack=True, usecols=[1,5,9])
	beta1, entropy1, error1 = data1

	plt.errorbar(beta1, entropy1, yerr=error1, color = 'b', ls = '-', label = 'PIGS', linewidth=1)
	plt.axhline(y=ExactValue, color='black', lw = 2.0, linestyle='--', label = 'Exact')

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


	if (variableName == "beta"):
		plt.xlabel(r'$\mathrm{\beta \ \  Kelvin^{-1}}$', fontsize = font)
		plt.text(textpositionx, textpositiony+0*deltay, 'Parameters:', fontsize=10)
		plt.text(textpositionx, textpositiony-1*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=10)
		plt.text(textpositionx, textpositiony-2*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=10)
		plt.text(textpositionx, textpositiony-3*deltay, r'$\mathrm{\mu} =$'+str(dipolemoment)+"Debye", fontsize=10)
		plt.text(textpositionx, textpositiony-4*deltay, r'$\mathrm{\tau} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=10)
	if (variableName == "tau"):
		plt.xlabel(r'$\mathrm{\tau \ \  Kelvin^{-1}}$', fontsize = font)
		plt.text(textpositionx, textpositiony+5*deltay, 'Parameters:', fontsize=10)
		plt.text(textpositionx, textpositiony+4*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=10)
		plt.text(textpositionx, textpositiony+3*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=10)
		plt.text(textpositionx, textpositiony+2*deltay, r'$\mathrm{\mu} =$'+str(dipolemoment)+"Debye", fontsize=10)
		plt.text(textpositionx, textpositiony+1*deltay, r'$\mathrm{\beta} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=10)
	
	plt.ylabel(r'$\mathrm{S_{2}}$', fontsize = font)


	plt.subplots_adjust(top=0.95, bottom=0.20, left=0.15, right=0.95, hspace=0.6, wspace=1.0)
	plt.legend(bbox_to_anchor=(0.70, 0.70), loc=2, borderaxespad=0.)
	plt.savefig(outfile, dpi = 200, format = 'pdf')

	call(["open", outfile])
	#plt.show()

def	FigurePIGS(FileToBePlot,FilePlot,TypeCorr,variableName,parameter,numbmolecules,molecule,Rpt,dipolemoment):
	outfile = FilePlot
	fig = plt.figure(figsize=(6, 4), dpi=200)

	font=18
	val1, val2, val3 = genfromtxt(FileToBePlot, unpack=True, usecols=[0, 1, 2], skip_footer = 10)

	plt.errorbar(val1, val2, yerr=val3, color = 'b', ls = '-', label = 'PIGS', linewidth=1)

	plt.grid(True)
	plt.xlim(1.98, 6.02)

	ymin, ymax = plt.ylim()
	midpointy = 0.5*(ymax-ymin)
	deltay = midpointy*0.15
	xmin, xmax = plt.xlim()
	midpointx = 0.5*(xmax-xmin)
	deltax = midpointx*0.15
	textpositionx = xmin+midpointx-0.25*midpointx
	textpositiony = ymin+midpointy

	if (variableName == "beta"):
		plt.xlabel(r'$\mathrm{j}$', fontsize = font)
		plt.text(textpositionx, textpositiony+0*deltay, 'Parameters:', fontsize=10)
		plt.text(textpositionx, textpositiony-1*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=10)
		plt.text(textpositionx, textpositiony-2*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=10)
		plt.text(textpositionx, textpositiony-3*deltay, r'$\mathrm{\mu} =$'+str(dipolemoment)+"Debye", fontsize=10)
		plt.text(textpositionx, textpositiony-4*deltay, r'$\mathrm{\tau} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=10)
		plt.text(textpositionx, textpositiony-5*deltay, r'$\mathrm{i} = 1$', fontsize=10)

	if (variableName == "tau"):
		plt.xlabel(r'$\mathrm{\tau \ \  Kelvin^{-1}}$', fontsize = font)
		plt.text(textpositionx, textpositiony+5*deltay, 'Parameters:', fontsize=10)
		plt.text(textpositionx, textpositiony+4*deltay, r'$\mathrm{N} =$'+str(numbmolecules)+" "+molecule, fontsize=10)
		plt.text(textpositionx, textpositiony+3*deltay, r'$\mathrm{R} =$'+str(Rpt)+r'$\mathrm{\AA}$', fontsize=10)
		plt.text(textpositionx, textpositiony+2*deltay, r'$\mathrm{\mu} =$'+str(dipolemoment)+"Debye", fontsize=10)
		plt.text(textpositionx, textpositiony+1*deltay, r'$\mathrm{\beta} =$' +str(parameter) +' '+"K"+r'$^{-1}$', fontsize=10)
		plt.text(textpositionx, textpositiony+0*deltay, r'$\mathrm{i} = 1$', fontsize=10)
	
	if (TypeCorr == "Total"):
		plt.ylabel(r'$\mathrm{C_{ij}}$', fontsize = font)
	if (TypeCorr == "XCorr"):
		plt.ylabel(r'$\mathrm{C^{X}_{ij}}$', fontsize = font)
	if (TypeCorr == "YCorr"):
		plt.ylabel(r'$\mathrm{C^{Y}_{ij}}$', fontsize = font)
	if (TypeCorr == "ZCorr"):
		plt.ylabel(r'$\mathrm{C^{Z}_{ij}}$', fontsize = font)
	if (TypeCorr == "XYCorr"):
		plt.ylabel(r'$\mathrm{C^{X,Y}_{ij}}$', fontsize = font)

	plt.subplots_adjust(top=0.95, bottom=0.20, left=0.17, right=0.95, hspace=0.6, wspace=1.0)
	plt.legend(bbox_to_anchor=(0.70, 0.70), loc=2, borderaxespad=0.)
	plt.savefig(outfile, dpi = 200, format = 'pdf')

	call(["open", outfile])
	#plt.show()
