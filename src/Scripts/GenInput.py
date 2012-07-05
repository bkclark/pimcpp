import sys
from math import *
from numpy import *
import subprocess

def GetUnique(a):
  seen = set()
  return [x for x in a if str(x) not in seen and not seen.add(str(x))]

# Input Parameters

## General Parameters
prefix = 'e-e' # filename prefix
unit1 = 'H' # energy unit
unit2 = 'A' # length unit
type1 = 'e' # particle type
type2 = 'e' # particle type
lam1 = 1.0 # hbar^2/2m
lam2 = 1.0 # hbar^2/2m

## Potgen
Z1Z2 = 1.0 # e^2
Npoints = 450 # number of grid points
r0 = 0.01 # first grid point
L = 16.0 # length of box
kCut = 5.0 # k cutoff for ewald
rCut = 5.0 # r cutoff for ewald
breakup = 1 # 1 - Optimized breakup, 0 - Classical Ewald breakup
diagonal = 1 # 1 - Do FT along diagonal, 0 - Do FT in x-direction

## PIMC Simulation
T = 0.025 # desired temperature of PIMC simulation
tau = 0.125 # desired timestep of PIMC simulation
M = int((1/T)/tau) # number of time slices
tau = (1/T)/M # used tau
print 'M', M
print 'tau', tau

## Squarer
ntemp = 1 # number of temperatures at which to compute the density matrix
tlow = 1/tau # lowest temperature at which the density matrix is generated
ndim = 3 # dimension
norder = 3 # highest order of polynomial fit
nl = 30 # number of partial waves
nsquare = 14 # total number of squarings to reach lowest temperature

# Generate Potential
print '**** Performing Breakup ****'
subprocess.call(['ewald',str(L),str(kCut),str(r0),str(rCut),str(Z1Z2),str(Npoints),str(breakup),str(diagonal)])

# Write .yk file
print '**** Creating Square Inputs ****'
kData = loadtxt('kData.txt')
ks = kData[:,0]
Vks = kData[:,1]
kData = GetUnique(kData.tolist())
kData.sort()
Vk0 = kData[0][1]
print 'Vk0', Vk0
f = open(prefix+'.yk','w')
for kDatum in kData:
  f.write('  %.10E'%kDatum[0]+'       %.10E'%kDatum[1]+'\n')
f.close()

# Write .dm and .in file for squarer
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

# Squarer
print '**** Performing Squaring Procedure ****'
subprocess.call(['squarer',prefix])

# Density Matrix Parser
print '**** Parsing Density Matrix ****'
subprocess.call(['python','dmparse.py',prefix+'.dm'])
