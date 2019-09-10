#!/usr/bin/python

import os
from os import system
import support

N=6
offset = 0
for i in range(0,N-1):
	for j in range(i+1,N):
		jj = j-i-1+offset
		print("i ", i, " j ",j," offset ",offset, jj)
	offset += (N-i-1)
