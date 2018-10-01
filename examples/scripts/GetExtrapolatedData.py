import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math

def fitFunc(t, a, b, c):
	return a + b*t*t + c*t*t*t

fileRead = "/home/tapas/ResultsOfENT/ENT-RotDOFs-Rpt10.05Angstrom-gFactor2.0-Entropy-vs-tau-fixed-beta0.2Kinv-Blocks20000-Passes100-System2HF-preskip10000-postskip0-SWAPTOUNSWAP.txt"
data      = np.loadtxt(fileRead)
xdata     = data[1:,1]
ydata     = data[1:,3]
ydata_err = data[1:,5]

plt.axhline(y=0.2627921404306721, color='black', lw = 2.0, linestyle='-', label = 'DMRG')
#plt.axhline(y=0.0780713110258388, color='black', lw = 2.0, linestyle='-', label = 'DMRG')
print(data[0,0])
plt.ylabel('Purity', fontsize = 16)
plt.xlabel('Imaginary time (K'+r'$^{-1}$'+')', fontsize = 16)
plt.xlim(-0.0001,0.021)
print(xdata)
print(ydata)
print(ydata_err)
0.196381214201262

# plot the data as red circles with errorbars in the vertical direction
popt, pcov = curve_fit(fitFunc, xdata, ydata, sigma=ydata_err)
print(popt)
perr       = np.sqrt(np.diag(pcov))
print(perr)
popt_up    = popt + perr
popt_dw    = popt - perr

t          = np.linspace(0, np.max(xdata), 10)
fit        = fitFunc(t, *popt)
fit_up     = fitFunc(t, *popt_up)
fit_dw     = fitFunc(t, *popt_dw)
fit_er     = fit_up - fit
print(fit[0])
print(fit_er[0])
# now plot the best fit curve and also +- 3 sigma curves
# the square root of the diagonal covariance matrix element 
# is the uncertianty on the corresponding fit parameter.

plt.errorbar(t, fit, yerr=fit_er, fmt = 'bv')
plt.errorbar(xdata, ydata, yerr=ydata_err, fmt = 'ro')
#plt.ylim(0,1)
plt.show()
