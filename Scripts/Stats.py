from numpy import *
import string
import re
from math import *

def C(g,t,mean,var,N):
  if (var == 0):
    return 0
  total = 0.0
  diff = g - mean
  total = sum(diff[0:N-t]*diff[t:N])
  return (1.0/(var*(N-t)))*total

def Kappa(g,mean,var):
  total = 0.0
  N = len(g)
  for i in range(1,N):
    c = C(g,i,mean,var,N)
    if (c <= 0):
      break
    else:
      total += c
  return 1.0 + 2.0*total

def Mean(g):
  return sum(array(g))/len(g)

def Mean2(g):
  return sum(array(g)**2)/len(g)

def RMS(g):
  return sqrt(Mean2(g))

def MSR(g):
  return sqrt(sum(array(g)**2))/len(g)

# Sample variance
def SampVar(g):
  return Mean2(g) - Mean(g)**2

# Unbiased estimate of population variance
def Var(g):
  N = len(g)
  if N > 1:
    return ((N+0.)/(N-1.)) * SampVar(g)
  else:
    return 0.

def Sigma(g):
  return sqrt(Var(g))

def Error(g):
  kappa = Kappa(g,Mean(g),Var(g))
  return Sigma(g)/sqrt(len(g)/kappa)

def StdError(g):
  return Sigma(g)/sqrt(len(g))

def UnweightedReAvg(stats):
  totMean,totMean2,totError,totN = 0.,0.,0.,0.
  for s in stats:
    [mean,error,N] = s[:3]
    mean2 = error*error*(N-1.) + mean*mean
    totMean += mean*N
    totMean2 += mean2*N
    totN += N
  if totN != 0:
    totMean /= totN
    totMean2 /= totN
    var = totMean2 - totMean*totMean
    if totN > 1 and var > 0:
      totError = sqrt(var/(totN-1.))
    else:
      totError = 0.
  return [totMean,totError,totN]

def UnweightedAvg(stats):
  if len(stats) > 0:
    mean = Mean([s[0] for s in stats])
    var = Mean2([s[1] for s in stats])
    stddev = sqrt(var)
    if len(stats[0]) == 4:
      kappa = sqrt(Mean2([s[2] for s in stats]))
      error = MSR([s[3] for s in stats])
      return [mean,stddev,kappa,error]
    else:
      error = stddev/sqrt(len([s[1] for s in stats]))
      return [mean,error]
  else:
    return [0.,0.,0.,0.]

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

def chunkIt(seq, num):
  avg = len(seq) / float(num)
  out = []
  last = 0.0
  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg
  return out

def BlockAnalysis(g):
  N = len(g)
  for n in range(1,int(N/2.)):
    h = chunkIt(g,n)
    stats = []
    for i in h:
      mean = sum(i)/len(i)
      mean2 = sum(i**2)/len(i)
      var = mean2 - mean**2
      kappa = Kappa(i,mean,var)
      Neff = N/kappa
      stats.append([mean,sqrt(var),kappa,sqrt(var/Neff)])
    #print n, UnweightedAvg(stats)

def stats(g):
  N = len(g)
  mean = Mean(g)
  mean2 = Mean2(g)
  if N > 1:
    var = ((N+0.)/(N-1.)) * (mean2-mean*mean)
  else:
    var = 0. #float('infinity')
  try:
    stddev = sqrt(var)
  except:
    print 'Warning: var < 0! Setting stddev = 0'
    var = 0.
    stddev = 0.
  kappa = Kappa(g,mean,var)
  Neff = N/kappa
  return [mean,stddev,kappa,stddev/sqrt(Neff)]

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
