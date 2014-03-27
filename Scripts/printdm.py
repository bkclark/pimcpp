import sys
import os
from numpy import *
import h5py as h5
from scipy import interpolate

f = h5.File(sys.argv[1])

def printitout(name):
  print name
  data = f[name+'/Data']
  data = data[:].flatten()
  print len(data)

  start = f[name+'/Grid/Start'][0]
  end = f[name+'/Grid/End'][0]
  N = f[name+'/Grid/NumGridPoints'][0]
  rs = linspace(start,end,num=N,endpoint=True)

  spline = interpolate.splrep(rs,data,s=0)
  rsnew = arange(0.1,3.01,0.1)
  datanew = interpolate.splev(rsnew,spline,der=0)

  for i in range(0,len(rsnew)):
    print rsnew[i], datanew[i]

printitout('Ukj0')
#printitout('Ukj1')
#printitout('Ukj2')
#printitout('Ukj3')
printitout('dUkjdBeta0')
#printitout('dUkjdBeta1')
#printitout('dUkjdBeta2')
#printitout('dUkjdBeta3')
f.close()
