#!/usr/bin/python
 
import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import support
import sys
import argparse

dipolemoment=1.86
Rpt=10.05
molecule="HF"
support.GetrAndgFactor(molecule, Rpt, dipolemoment)
support.GetPreFactDDPot(molecule, Rpt, dipolemoment)
