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

dipolemoment=1.8548
Rpt=6.0
molecule="HF"
support.GetrAndgFactor(molecule, Rpt, dipolemoment)
