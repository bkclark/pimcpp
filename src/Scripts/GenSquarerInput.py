import sys
from math import *
from numpy import *
import subprocess

# Input Parameters

## General Parameters
prefix = 'e-e' # filename prefix
unit1 = 'r' # energy unit
unit2 = 'a' # length unit
type1 = 'e' # particle type
type2 = 'e' # particle type
lam1 = 0.0625 # hbar^2/2m
lam2 = 0.0625 # hbar^2/2m

## Potgen
Z1Z2 = 0.5 # e^2
Npoints = 450 # number of grid points
r0 = 0.01 # first grid point
L = 5.17051953218 # length of box
kCut = 4.86077676958 # k cutoff for ewald
rCut = 2.51215075962 # r cutoff for ewald
breakup = 1 # 1 - Optimized breakup, 0 - Classical Ewald breakup
diagonal = 1 # 1 - Do FT along diagonal, 0 - Do FT in x-direction

## Squarer
tlow = 0.0913541505403 # lowest temperature at which the density matrix is generated
ntemp = 8 # number of temperatures at which to compute the density matrix
ndim = 3 # dimension
norder = 3 # highest order of polynomial fit
nl = 30 # number of partial waves
nsquare = 14 # total number of squarings to reach lowest temperature

# Generate Potential
subprocess.call(['ewald',str(L),str(kCut),str(r0),str(rCut),str(Z1Z2),str(Npoints),str(breakup),str(diagonal)])
kData = loadtxt('kData.txt')
Vk0 = kData[0][1]

# Write .dm file for squarer
f = open(prefix+'.dm','w')
f.write(' UNITS '+unit1+' '+unit2)
f.write('\n TYPE '+type1+' '+str(lam1))
f.write('\n TYPE '+type2+' '+str(lam2))
f.write('\n GRID '+str(Npoints)+' LINEAR '+str(r0)+' '+str(L))
f.write('\n SQUARER '+str(tlow)+' '+str(ntemp)+' '+str(ndim)+' '+str(norder)+' '+str(nl)+' '+str(nsquare))
f.write('\n POT COULOPT '+str(rCut)+' '+str(kCut)+' '+str(Z1Z2)+' '+str(ndim)+' 0.D0 1.0 '+str(L)+' '+str(L)+' '+str(L))
f.write('\n VIMAGE '+str(Vk0))
f.write('\n POTTAIL 0.0')
f.write('\n RANK 2 '+str(Npoints)+' 1')
f.write('\n GRID '+str(Npoints)+' LINEAR '+str(r0)+' '+str(L))
f.write('\n LABEL 1 r')
f.write('\n BEGIN potential 0')
rData = loadtxt('rData.txt')
count = 0
f.write('\n  ')
for [r,VShort,VLong] in rData:
  f.write('%.10E'%(VShort)+'  ')
  count += 1
  if count % 5 == 0:
    f.write('\n  ')
f.close()
