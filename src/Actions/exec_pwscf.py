#!/usr/bin/python

from pwscf import *
import sys

infile_pre = 'pwscf.0.in'
infile_post = 'pwscf.1.in'

if(len(sys.argv) != 3):
  print "Usage: >./exec_pwscf.py <coordfile.dat> <energy.out>"
  sys.exit()
infile_coords = sys.argv[1]
energy_out = sys.argv[2]
execute(infile_pre, infile_post, infile_coords, energy_out)
