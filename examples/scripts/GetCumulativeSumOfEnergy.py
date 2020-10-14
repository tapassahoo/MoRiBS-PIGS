import os
import sys
from subprocess import call
import math
from math import *
import numpy as np
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.mlab as mlab
from itertools import cycle
import itertools

fig=plt.figure()

file1=sys.argv[1]
print(file1)

if (os.path.isfile(file1+"_old")==True):
	print(file1+"_old")
	col_data_new=np.genfromtxt(file1)
	index=int(col_data_new[0,0])
	col_data_old=np.genfromtxt(file1+"_old")
	merged_data=np.concatenate((col_data_old[:index-1], col_data_new), axis=0)
else:
	merged_data=np.genfromtxt(file1)

col=merged_data[:,0]
x=merged_data[:,3]
engCumSum=np.cumsum(x)
Norm=np.arange(1,int(len(engCumSum)+1),dtype = 'float')
engCumSum=np.divide(engCumSum,Norm)


print("Length of simulation is ")
print(len(engCumSum))
print("")

plt.xlim(-1,len(col)+500)
plt.plot(col,engCumSum,color='black',ls='-',linewidth=1,marker='None',markersize=9,label='Total Energy')
plt.xlabel('Number of blocks')
plt.ylabel('Cumulative sum')

plt.minorticks_on()
plt.tick_params(axis="both", direction="in", which="minor", right=True, top=False, length=2)
plt.tick_params(axis="both", direction="in", which="major", right=True, top=True, length=6)

plt.legend(numpoints=1,loc=('center right'))
plt.subplots_adjust(top=0.99,bottom=0.12,left=0.15,right=0.98,hspace=0.0,wspace=0.0)

plt.show()
