"""
Simple demo with multiple subplots.
"""
import numpy as np
from numpy import *
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6), dpi=80)
col_beta1, col_pot1, col_tot1, col_err_pot1, col_err_tot1 = loadtxt('Energy-vs-beta-2-HF-fixed-tau0.001-blocks80000.txt',unpack=True, usecols=[1,2,3,5,6])
col_beta2, col_pot2, col_tot2, col_err_pot2, col_err_tot2 = loadtxt('Energy-vs-beta-2-HF-fixed-tau0.002-blocks80000.txt',unpack=True, usecols=[1,2,3,5,6])
plt.subplot(2, 1, 1)
plt.plot(col_beta1, col_tot1, 'o-', label= r'$\tau$ = 0.001 K $^{-1}$')
plt.errorbar(col_beta1, col_tot1, yerr=col_err_tot1, fmt='o')
plt.plot(col_beta2, col_tot2, 'o-', label= r'$\tau$ = 0.002 K $^{-1}$')
plt.errorbar(col_beta2, col_tot2, yerr=col_err_tot2, fmt='o')
plt.title(r'$\beta$' ' convergence')
#plt.text(0.1, -2.0, r'$\tau$ = 0.001 K $^{-1}$')
plt.ylabel('E'r'$_{0}  [K^{-1}]$')
#plt.grid(True)
plt.xlim((0,0.201))
plt.xticks(np.linspace(-0, 0.2, 5, endpoint=True))
plt.yticks(np.linspace(-4.5, -0.5, 5, endpoint=True))
plt.axhline(y=-3.3248566, color='black', lw = 1.0, linestyle='--')
plt.legend(loc='upper right')
#plt.ylim((-3.5,-3.2))

#---------------
#
#second pot ----
#
#---------------
plt.subplot(2, 1, 2)

plt.plot(col_beta1, col_pot1, 'o-', label = r'$\tau$ = 0.001 K $^{-1}$')
plt.errorbar(col_beta1, col_pot1, yerr=col_err_pot1, fmt='o')
plt.plot(col_beta2, col_pot2, 'o-', label = r'$\tau$ = 0.002 K $^{-1}$')
plt.errorbar(col_beta2, col_pot2, yerr=col_err_pot2, fmt='o')

#plt.text(0.1, -2.0, r'$\tau$ = 0.001 K $^{-1}$')
plt.ylabel('V'r'$ [K^{-1}]$')
plt.xlabel(r'$\beta$' r'$ [K^{-1}]$')

#plt.grid(True)
plt.xlim((0,0.201))
plt.xticks(np.linspace(-0, 0.2, 5, endpoint=True))
plt.axhline(y=-3.3248566, color='black', lw = 1.0, linestyle='--')
plt.legend(loc='upper right')

plt.show()

