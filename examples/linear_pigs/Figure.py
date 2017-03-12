import numpy as np
from numpy import *
import matplotlib.pyplot as plt
from itertools import cycle
import itertools
from scipy.optimize import curve_fit

#===============================================================================
#                                                                              |
#   Some parameters for submission of jobs and analysis outputs.               |
#   Change the parameters as you requied.                                      |
#                                                                              |
#===============================================================================
#molecule            = "HF-C60"                                                     #change param1
molecule            = "HF"                                                     #change param1
#molecule            = "H2"                                                         #change param1
molecule_rot        = "HF"                                                         #change param2

numbblocks          = 20000                                                        #change param3
numbmolecules       = 1                                                            #change param4
nrange              = 1                                                            #change param5

Rpt                 = 3.5                                                          #change param6
#values              = [0.128, 0.256, 0.512]                                        #change param8
values              = [0.001, 0.002, 0.004, 0.01]                                  #change param7
values              = [0.004]                                  #change param7
values              = [0.128]
Avg_total_energy    = -11.823                                                      #change param9
Avg_potential_energy= -21.1273                                                     #change param10
var2                = "beta"                                                       #change param11
var1                = "tau"                                                        #change param12
num1                = 2
trunc 				= 7

#======================================================================================================================

class Specification:
	def __init__(self):
		self.color1               = itertools.cycle(('r', 'b', 'g', 'm'))
		self.marker1              = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))
		self.line1                = itertools.cycle(('-','--','-.',':'))
		self.color2               = itertools.cycle(('r', 'b', 'g', 'm'))
		self.marker2              = itertools.cycle(('o', 'v', '^', '<', '>', 's', '8', 'p'))
		self.line2                = itertools.cycle(('-','--','-.',':'))

def plotenergy(func,var1,var2,var3,numbmolecules,molecule,numbblocks,x,trunc):
	markersize_fig = 1
	srcfile     = "Energy-vs-"+var1+"-"+str(numbmolecules)+"-"+str(molecule)+"-fixed-"+var2+str(var3)+"-blocks"+str(numbblocks)+".txt"
	label_fig   = r'$\tau$ = '+str(var3)+' K $^{-1}$'
	data        = loadtxt(srcfile,unpack=True, usecols=[1,2,3,5,6])
	var, pot, tot, err_pot, err_tot = data
	var = var[0:trunc]
	pot = pot[0:trunc]
	tot = tot[0:trunc]
	err_pot = err_pot[0:trunc]
	err_tot = err_tot[0:trunc]

	if func == 'Energy':
		val     = tot
		err_val = err_tot
	if func == 'Potential':
		val     = pot
		err_val = err_pot
		
	plt.plot(var, val, linestyle = x.line1.next(), color = x.color1.next(), label = label_fig, lw = 2)
	plt.errorbar(var, val, yerr = err_val, ecolor = x.color2.next(), fmt = x.marker2.next(), markersize = markersize_fig)

def fitFunc(var, a, b):
	return a + b*var*var

def plotfitting(func,Avg_energy,b,var1,var2,var3,numbmolecules,molecule,numbblocks,x,trunc):
	markersize_fig = 10
	srcfile     = "Energy-vs-"+var1+"-"+str(numbmolecules)+"-"+str(molecule)+"-fixed-"+var2+str(var3)+"-blocks"+str(numbblocks)+".txt"
	label_fig   = r'$\tau$ = '+str(var3)+' K $^{-1}$'
	data        = loadtxt(srcfile,unpack=True, usecols=[1,2,3,5,6])
	var, pot, tot, err_pot, err_tot = data
	var = var[0:trunc]
	pot = pot[0:trunc]
	tot = tot[0:trunc]
	err_pot = err_pot[0:trunc]
	err_tot = err_tot[0:trunc]

	if func == 'Energy':
		val     = tot
		err_val = err_tot
	if func == 'Potential':
		val     = pot
		err_val = err_pot
		
	fitParams, fitCovariances = curve_fit(fitFunc, var, val, sigma = err_val)

	plt.plot(var, fitFunc(var, fitParams[0], fitParams[1]), linestyle = x.line1.next(), color = 'g', label = label_fig, lw = 2)


def plotenergylebel(num1,num2,Avg_energy):
	if num1 == num2:
		#plt.text(0.1, -2.0, r'$\tau$ = 0.001 K $^{-1}$')
		plt.ylabel('V'r'$ [K^{-1}]$')
		plt.xlabel(r'$\beta$' r'$ [K^{-1}]$')

		#plt.grid(True)
		#plt.xlim((0,0.201))
		#plt.xticks(np.linspace(-0, 0.2, 5, endpoint=True))
		plt.axhline(y=Avg_energy, color='black', lw = 3.0, linestyle='--')
		plt.legend(loc='upper right')
	if num1 != num2:
		plt.title(r'$\beta$' ' convergence for 1 HF')
		#plt.text(0.1, -2.0, r'$\tau$ = 0.001 K $^{-1}$')
		plt.ylabel('E'r'$_{0}  [K^{-1}]$')
		#plt.grid(True)
		#plt.xlim((0,0.205))
		#plt.xticks(np.linspace(-0.1, 1.4, 5, endpoint=True))
		#plt.yticks(np.linspace(-4.5, -0.5, 5, endpoint=True))
		plt.axhline(y=Avg_energy, color='black', lw = 3.0, linestyle='--')
		plt.legend(loc='upper right')
		#plt.ylim((-3.5,-3.2))


plt.figure(figsize=(4, 8), dpi=100)

num2 = 1
plt.subplot(num1, 1, num2)
x = Specification()
func = "Energy"
for i in range(nrange):
	var3     = values[i]
	plotenergy(func,var1,var2,var3,numbmolecules,molecule,numbblocks,x,trunc)
	plotenergylebel(num1,num2,Avg_total_energy)

# fitting done here
	if var1 == 'tau':
		b = -100.0
		plotfitting(func,Avg_total_energy,b,var1,var2,var3,numbmolecules,molecule,numbblocks,x,trunc)

num2 = 2
plt.subplot(num1, 1, num2)
x = Specification()
func = "Potential"
for i in range(nrange):
	var3     = values[i]
	plotenergy(func,var1,var2,var3,numbmolecules,molecule,numbblocks,x,trunc)
	plotenergylebel(num1,num2,Avg_potential_energy)

# fitting done here
	if var1 == 'tau':
		b = -100.0
		plotfitting(func,Avg_potential_energy,b,var1,var2,var3,numbmolecules,molecule,numbblocks,x,trunc)

plt.show()
