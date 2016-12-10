#!/usr/bin/python
# Example PBS cluster job submission in Python
 
from popen2 import popen2
import time
import os

file="yw*"
cmd_rm="rm "+file
os.system(cmd_rm)
 
    # Open a pipe to the qsub command.
output, input = popen2('qsub')
     
    # Customize your options here
job_name = "test-py"
command = "time ./pimc"
 
job_string = """#!/bin/bash
             #PBS -N %s
             cd $PBS_O_WORKDIR
             %s""" % (job_name, command)
     
    # Send job_string to qsub
input.write(job_string)
input.close()
     
    # Print your job and the system response to the screen as it's submitted
print job_string
print output.read()
     
time.sleep(0.1)
