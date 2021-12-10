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
from subprocess import call
import getpass

import numpy as np

import mypkg.pkgAsymrho.support_asymrho as asym

parser = argparse.ArgumentParser(
    description="It is used to submit jobs in a queue for generating the rotational density matrix elements of an asymmetric rotor. Note: Module support_asymrho.py consists of many functions and it is not permitted to modify without consulting the developer - Dr. Tapas Sahoo." 
)

parser.add_argument(
    "-group",
    "--Partition",
    help="It allows user to submit jobs in a specific CPU. Note that user does not need flag --Partition or -group.",
    type=str,
    default="ntapas",
    choices=["ntapas", "tapas"],
)

parser.add_argument(
    "-P",
    "--Beads",
    help="Number of beads for Path Integral Monte Carlo simulation.",
    type=int,
    default=4,
)

parser.add_argument(
    "-iodevn",
    "--spin",
    help="Spin coupling. It is required to generate rotational density matrix for a particular spin isomer - -ortho, -para, or spinless. For more informations read the README file placed in 'nmv_prop' directory.",
    type=int,
    default=-1,
    choices=[-1, 0, 1],
)

parser.add_argument("-J", "--Jmax", help="Maximum J quantum number.", type=int, default=0)

parser.add_argument(
    "-job",
    "--TypeCal",
    help="It is either submission for analysis. A and S atand for submission and analysis, respectively.",
    type=str,
    default="A",
    choices=["S", "A"],
)

parser.add_argument(
    "Rotor",
    help="Name of a rotor. E.g. - H2O. It is needed to save rotational density matrix.",
    type=str,
)

parser.add_argument(
    "param", help="Name of parameter: ether beta or tau.", type=str, choices=["tau", "beta"]
)

parser.add_argument(
    "value",
    help="Value of the parameter. It is the value of either beta or tau.",
    type=float,
    default=1.0,
)
args = parser.parse_args()
# Argparser ends here

myhost = os.uname()[1]
myhost = myhost[0:3]
if ((myhost == "gra") or (myhost == "ced")):
    NameOfServer = "graham"
else:
    NameOfServer = "nlogn"

'''
if ((myhost == "gra-login1") or (myhost == "gra-login2") or (myhost == "gra-login3") or (myhost == "gra768") or (myhost == "gra796")):
    NameOfServer = "graham"
elif ((myhost == "cedar1.cedar.computecanada.ca") or (myhost == "cedar2.cedar.computecanada.ca") or (myhost == "cedar3.cedar.computecanada.ca") or (myhost == "cedar4.cedar.computecanada.ca") or (myhost == "cedar5.cedar.computecanada.ca")):
    NameOfServer = "graham"
else:
    NameOfServer = "nlogn"
'''

NameOfPartition = args.Partition

rotor = args.Rotor
param = args.param
numbbeads = args.Beads
jmax = args.Jmax
iodevn = args.spin
TypeCal = args.TypeCal

if param == "tau":
    beta = args.value * numbbeads
else:
    beta = args.value

temperature = 1.0 / beta
temperature = "%8.6f" % temperature

extra_file_name = ""
script_dir = os.getcwd()
user_name = getpass.getuser()
home = os.path.expanduser("~")

if user_name == "tsahoo":
    dir_user = "/Users/tsahoo/mount_feynman/"
else:
    dir_user = home + "/"
moribs_dir = "MoRiBS-PIGS/"
asymrho_dir = "nmv_prop/"
dir_store = "rot-dens-asymmetric-top/"

src_dir_exe = dir_user + moribs_dir + asymrho_dir
if TypeCal == "S":
    asym.MakeExecutable(script_dir, src_dir_exe)

dir_name = asym.GetDirNameSubmission(rotor, temperature, numbbeads, iodevn, jmax)

if NameOfServer == "graham":
    dir_job = "/scratch/" + user_name + "/" + dir_store + dir_name
else:
    dir_job = "/work/" + user_name + "/" + dir_store + dir_name

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
    asym.Submission(dir_job, script_dir, execution_file, numbbeads, temperature, rotor, dir_input, dir_output, NameOfPartition, user_name, dir_store, dir_name, iodevn, jmax, NameOfServer)
else:
    asym.GetPackRotDens(src_dir_exe, dir_output, script_dir, rotor, temperature, numbbeads)
