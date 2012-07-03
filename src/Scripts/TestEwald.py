import sys
import subprocess
from math import *
from numpy import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os.path

L = 10.0
Params = [[5.0,3.0,1],[5.0,3.0,0],[5.0,4.0,1],[5.0,4.0,0],[5.0,3.5,1],[5.0,3.5,0]]
Z1Z2 = 1.0
NumPoints = 1000
diagonal = 1
rerun = 1
r0 = 0.001

for [kCut,rCut,breakup] in Params:
  filename='FTData-'+str(L)+'-'+str(kCut)+'-'+str(r0)+'-'+str(rCut)+'-'+str(Z1Z2)+'-'+str(NumPoints)+'-'+str(breakup)+'-'+str(diagonal)
  if (not os.path.isfile(filename)) or rerun:
    subprocess.call(['ewald',str(L),str(kCut),str(r0),str(rCut),str(Z1Z2),str(NumPoints),str(breakup),str(diagonal)])
    subprocess.call(['cp','FTData.txt',filename])
  a=loadtxt(filename,skiprows=3)
  if breakup:
    plt.plot(a[:,0],a[:,2]+a[:,1]-1./a[:,0],'--',label=filename)
  else:
    plt.plot(a[:,0],a[:,2]+a[:,1]-1./a[:,0],'-',label=filename)
#plt.xlim([0,L/2])
#plt.ylim([-1,1])
plt.legend()
plt.savefig('comparison.png')
plt.clf()




