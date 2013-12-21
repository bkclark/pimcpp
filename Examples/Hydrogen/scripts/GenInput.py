import sys
from math import *
from numpy import *
import subprocess

def GetUnique(a):
  seen = set()
  return [x for x in a if str(x) not in seen and not seen.add(str(x))]

# Input Parameters
params = [['e','e',0.5,0.5,1.0],['e','p',0.5,0.0002723089072243553,-1.0],['p','p',0.0002723089072243553,0.0002723089072243553,1.0]]

for [type1,type2,lam1,lam2,Z1Z2] in params:
  print type1, lam1, type2, lam2, 'Z1Z2', Z1Z2

  ## General Parameters
  prefix = type1+'-'+type2 # filename prefix
  unit1 = 'H' # energy unit
  unit2 = 'A' # length unit

  ## Potgen
  Npoints = 450 # number of grid points
  r0 = 0.1 # first grid point
  L = 20.0 # length of box
  kCut = 3.0 # k cutoff for ewald
  rCut = 8.0 # r cutoff for ewald
  breakup = 2 # 2 - Short-ranged only, 1 - Optimized breakup, 0 - Classical Ewald breakup
  gridType = "LINEAR" # LOG/LINEAR
  diagonal = 1 # 1 - Do FT along diagonal, 0 - Do FT in x-direction

  ## PIMC Simulation
  T = 0.025 # desired temperature of PIMC simulation
  tau = 0.125 # desired timestep of PIMC simulation
  M = int((1/T)/tau) # number of time slices
  tau = (1/T)/M # used tau
  print 'M', M
  print 'tau', tau
  print 'beta', M*tau

  ## Squarer
  ntemp = 1 # number of temperatures at which to compute the density matrix
  tlow = (1./tau)/(1<<(ntemp-1)) # lowest temperature at which the density matrix is generated
  ndim = 3 # dimension
  norder = 3 # highest order of polynomial fit
  nl = 30 # number of partial waves
  nsquare = 14 + ntemp-1 # total number of squarings to reach lowest temperature

  # Generate Potential
  print '**** Performing Breakup ****'
  if breakup!=2:
    subprocess.call(['ewald',str(L),str(kCut),str(r0),str(rCut),str(Z1Z2),str(Npoints),str(breakup),str(diagonal)])
    rMax = 0.75*sqrt(3)*L
  else:
    f = open('rData.txt','w')
    if gridType=="LINEAR":
      rs = linspace(r0,L,num=Npoints,endpoint=True)
    elif gridType=="LOG":
      rs = logspace(log10(r0),log10(L),num=Npoints,endpoint=True)
    else:
      print 'Unrecognized grid'
    for r in rs:
      f.write(str(r)+' '+str(Z1Z2/r)+' 0.0\n')
    f.close()
    rMax = L

  print '**** Creating Squarer Inputs ****'
  # Write .yk file
  Vk0 = 0.0
  if breakup!=2:
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
  f.write('\n GRID '+str(Npoints)+' '+gridType+' '+str(r0)+' '+str(rMax))
  f.write('\n SQUARER '+str(tlow)+' '+str(ntemp)+' '+str(ndim)+' '+str(norder)+' '+str(nl)+' '+str(nsquare))
  f.write('\n POT COULOPT '+str(rCut)+' '+str(kCut)+' '+str(Z1Z2)+' '+str(ndim)+' 0.D0 1.0 '+str(L)+' '+str(L)+' '+str(L))
  f.write('\n VIMAGE '+str(Vk0))
  f.write('\n POTTAIL 0.0')
  f.write('\n RANK 2 '+str(Npoints)+' 1')
  f.write('\n GRID 1 '+gridType+' '+str(r0)+' '+str(rMax))
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
  f.write('\n')
  f.close()

  # Squarer
  print '**** Performing Squaring Procedure ****'
  subprocess.call(['squarer',prefix])

  # Density Matrix Parser
  print '**** Parsing Density Matrix ****'
  subprocess.call(['python','../scripts/dmparse.py',prefix+'.dm'])

  # Density Matrix Parser
  print '**** Creating PairAction File ****'
  f = open(prefix+'.PairAction','w')
  f.write('Section (Fits)')
  f.write('\n{')
  f.write('\n  string Type="DavidFit";')
  f.write('\n  int NumOffDiagonalTerms = 3;')
  f.write('\n  Section (Particle1)')
  f.write('\n  {')
  f.write('\n    string Name = "'+type1+'";')
  f.write('\n    double lambda ='+str(lam1)+';')
  f.write('\n    int Ndim = 3;')
  f.write('\n  }')
  f.write('\n  Section (Particle2)')
  f.write('\n  {')
  f.write('\n    string Name = "'+type2+'";')
  f.write('\n    double lambda ='+str(lam2)+';')
  f.write('\n    int Ndim = 3;')
  f.write('\n  }')
  f.write('\n  string Daviddmfile = "inputs/'+prefix+'.h5";')
  f.write('\n}')
  f.close()
