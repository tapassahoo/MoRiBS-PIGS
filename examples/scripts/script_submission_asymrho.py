# Copyright 2019 Dr. Tapas Sahoo
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import decimal
import os
import sys
import time
from os import system
from subprocess import call

import numpy as np

import support_asymrho as asym
import inputFile

parser = argparse.ArgumentParser(
    description="It is a script file, written in Python, used to submit jobs in a queue for generating the rotational density matrix of an asymmetric rotor. Note: Module support_asymrho.py consists of many functions and it is not permitted to modify without consulting the developer - Dr. Tapas Sahoo. User can easily modify module inputFile.py to generate lists of beads (see Getbeads function), step lengths for rotational and translational motions, and levels for Bisection move (see class GetStepAndLevel) as needed."
)
parser.add_argument(
    "--Partition",
    help="It allows user to submit jobs in a specific CPU.",
    type=str,
    default="ntapas",
)

parser.add_argument(
    "--Rotor",
    help="Name of a rotor. E.g. - H2O. It is needed to save rotational density matrix.",
    type=str,
)

parser.add_argument(
    "--Temperature",
    help="Temperature in Kelvin. It is required to generate rotational density matrix.",
    type=float,
    default=0.0,
)

parser.add_argument(
    "--Beads", help="Number of beads for Path Integral Monte Carlo simulation.", type=int, default=4
)

parser.add_argument(
    "--spin",
    help="Spin coupling. It is required to generate rotational density matrix.",
    type=int,
    default=-1,
)

parser.add_argument("--Jmax", help="Maximum J quantum number.", type=int, default=0)
parser.add_argument(
    "--TypeCal", help="It is either submission for analysis.", type=str, default="A"
)
args = parser.parse_args()
#Argparser ends here

myhost = os.uname()[1]
if ((myhost == 'gra-login1') or (myhost == 'gra-login2') or (myhost == 'gra-login3')):
	NameOfServer = "graham"
else:
	NameOfServer = "nlogn"
NameOfPartition = args.Partition

rotor = args.Rotor
temperature = args.Temperature
numbbeads = args.Beads
jmax = args.Jmax
iodevn = args.spin
TypeCal = args.TypeCal


# RUNDIR              = "work"
RUNDIR = "scratch"
extra_file_name = ""
script_dir = os.getcwd()
user_name = os.getlogin()
home = os.path.expanduser('~')

if user_name == "tsahoo":
    dir_user = "/Users/tsahoo/mount_feynman/"
else:
    dir_user = home+"/"
moribs_dir = "MoRiBS-PIGS/"
asymrho_dir = "nmv_prop/"
dir_store = "rot-dens-asymmetric-top/"

src_dir_exe = dir_user + moribs_dir + asymrho_dir
if TypeCal == "S":
    asym.MakeExecutable(script_dir, src_dir_exe)

dir_name = asym.GetDirNameSubmission(rotor, temperature, numbbeads, iodevn, jmax)

if (RUNDIR == "scratch"):
    dir_job = "/scratch/" + user_name + "/"
if (NameOfServer == "graham"):
    dir_job = "/scratch/" + user_name +"/"+dir_store+dir_name

if NameOfServer == "graham":
    dir_input = "/scratch/" + user_name + "/" + dir_store + dir_name
    dir_output = "/scratch/" + user_name + "/" + dir_store + dir_name
else:
    dir_input = "/work/" + user_name + "/" + dir_store + dir_name
    dir_output = "/work/" + user_name + "/" + dir_store + dir_name

if os.path.isdir(dir_input) == False:
    call(["rm", "-rf", dir_input])
    call(["mkdir", "-p", dir_input])

execution_file = src_dir_exe + "asymrho.x"
call(["cp", execution_file, dir_input])

if TypeCal == "S":
    asym.Submission(
        RUNDIR,
        dir_job,
        script_dir,
        execution_file,
        numbbeads,
        temperature,
        rotor,
        dir_input,
        dir_output,
        NameOfPartition,
        user_name,
        dir_store,
        dir_name,
        iodevn,
        jmax,
        NameOfServer,
    )
else:
    asym.GetPackRotDens(src_dir_exe, dir_output, script_dir,rotor,temperature,numbbeads)
