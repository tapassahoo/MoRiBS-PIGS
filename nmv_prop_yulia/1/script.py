import numpy as np

import subprocess
from os import system
for jmax in range(35,100):
#        print(jmax)
        command="./asymrho.x 1. 64 -1 10 10 0.6666525 0.2306476 0.1769383 "+str(jmax)+"|grep \"CLASSICAL: Z\"" 
        system(command)
      

