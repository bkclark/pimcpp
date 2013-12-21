from numpy import *
import string
import re
from math import *

def C(g,t,mean,sigma2,N):
  if (sigma2 == 0):
    return 0
  total = 0.0
  diff = g - mean
  total = sum(diff[0:N-t]*diff[t:N])
  return (1.0/(sigma2*(N-t)))*total

def Kappa(g,mean,sigma2,N):
  total = 0.0
  for i in range(1,len(g)):
    c = C(g,i,mean,sigma2,N)
    if (c <= 0):
      break
    else:
      total += c
  return 1.0 + 2.0*total

def Mean(g):
  return sum(g)/len(g)

def Mean2(g):
  return sum(g**2)/len(g)

def Var(g):
  mean = Mean(g)
  mean2 = Mean(g**2)
  return mean2 - mean**2

def Sigma(g):
  return sqrt(Var(g))

def Error(g):
  kappa = Kappa(g,Mean(g),Var(g),len(g))
  return Sigma(g)/sqrt(len(g)/kappa)

def StdError(g):
  return Sigma(g)/sqrt(len(g))

def UnweightedAvg(stats):
    means = [s[0] for s in stats]
    stddevs = [s[1] for s in stats]
    kappas = [s[2] for s in stats]
    errors =[s[3] for s in stats]
    mean=sum(array(means))/(len(means)+0.0)
    stddev=sqrt(sum(array(stddevs)**2))/len(stddevs)
    kappa=sqrt(sum(array(kappas)**2))/len(kappas)
    error=sqrt(sum(array(errors)**2))/len(errors)
    return [mean,stddev,kappa,error]

def WeightedAvg(means, errors):
    zeroErrors = False
    for i in errors:
        if i == 0.0:
            zeroErrors = True
    if (not zeroErrors):
        weights = map (lambda x: 1.0/(x*x), errors)
        norm = 1.0/sum(weights)
        weights = map(lambda x: x*norm, weights)
        avg = 0.0
        error2 = 0.0
        for i in range (0,len(means)):
            avg = avg + means[i]*weights[i]
            error2 = error2 + weights[i]**2*errors[i]*errors[i]
        return (avg, math.sqrt(error2))
    else:
        return (sum(means)/len(means), 0.0)

def stats(g):
  N = len(g)
  mean = sum(g)/N
  mean2 = sum(g**2)/N
  var = mean2 - mean**2
  if abs(var) < 1e-12:
    #print 'Warning: var < 1e-12! Setting var = 0'
    var = 0.
  try:
    sigma = sqrt(var)
  except:
    print 'Warning: var < 0! Setting sigma = 0'
    var = 0.
    sigma = 0.
  kappa = Kappa(g,mean,var,N)
  Neff = N/kappa
  return [mean,sigma,kappa,sigma/sqrt(Neff)]

def getAndOutputStats(data):
  [myMean,myStdDev,myKappa,myStdErr] = stats(data)
  print 'Mean     = %8.5f' % (myMean) 
  print 'Standard Deviation     = %8.5f' % (myStdDev) 
  print 'Kappa     = %8.5f' % (myKappa) 
  print 'Standard Error    = %8.5f' % (myStdErr)
  return [myMean,myStdDev,myKappa,myStdErr]

import sys
def main():
  try:
    filename = str(sys.argv[1])
  except:
    print 'Need a file name.'
    abort()
  try:
    print 'Reading data from ', filename
    data = transpose(loadtxt(filename))
  except:
    print 'Data file not properly formed.'
    abort()
  try:
    startCut = int(sys.argv[2])
  except:
    startCut = 0
  if len(data.shape) == 1:
    assert startCut < len(data)
    getAndOutputStats(data[startCut:-1])
  else:
    assert startCut < data.shape[1]
    for datum in data:
      getAndOutputStats(datum[startCut:-1])

if __name__ == "__main__":
  main()
