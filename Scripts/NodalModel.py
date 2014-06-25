import sys
import numpy as np
import h5py as h5
import Stats
import matplotlib as mpl
mpl.use('Agg')
#import matplotlib.pyplot as plt

nType = int(sys.argv[1])
files = sys.argv[2:]
StartCut = 10
stats = {}
for iType in range(nType):
  for fname in files:
    print fname
    f = h5.File(fname,'r')

    fails = True
    while fails:
      try:
        if nType > 1:
          Ts = np.array(f['Observables/NodalModelTime_'+str(iType)+'/Time'])
        else:
          Ts = np.array(f['Observables/NodalModelTime/Time'])
        print Ts
        fails = False

        Ts = Ts.T
        NumModels = len(Ts)
        for i in range(0,NumModels):
          try:
              stats[i].append(Stats.stats(Ts[i][StartCut:]))
          except:
              stats[i] = [Stats.stats(Ts[i][StartCut:])]
      except:
        fails = False

  if nType > 1:
    f = open('NodalModelTime_'+str(iType)+'.dat','w')
  else:
    f = open('NodalModelTime.dat','w')
  for i in stats.iterkeys():
    tmpStats = Stats.UnweightedAvg(stats[i])
    for stat in tmpStats:
      f.write(str(stat)+' ')
    f.write('\n')
  f.close()

#plt.xlabel("Model")
#plt.ylabel("Relative Time")
#plt.errorbar(range(0,5),means,yerr=errs,fmt='o')
#plt.savefig('ftest.png')
#plt.clf()

#  AcceptArray = np.zeros((NumModels,NumModels))
#  AttemptArray = np.zeros((NumModels,NumModels))
#  AcceptRatioArray = np.zeros((NumModels,NumModels))
#  for fname in sys.argv[2:]:
#    f = h5.File(fname,'r')
#
#    if nType > 1:
#      AcceptArray += sum(1.0*np.array(f['Moves/NodalModel_'+str(iType)+'/Stage/AcceptRatios'])[StartCut:])
#      AttemptArray += sum(1.0*np.array(f['Moves/NodalModel_'+str(iType)+'/Stage/AttemptRatios'][StartCut:]))
#    else:
#      AcceptArray += sum(1.0*np.array(f['Moves/NodalModel/Stage/AcceptRatios'])[StartCut:])
#      AttemptArray += sum(1.0*np.array(f['Moves/NodalModel/Stage/AttemptRatios'][StartCut:]))
#
#  for i in range(NumModels):
#      for j in range(NumModels):
#          if AttemptArray[i][j] != 0:
#              AcceptRatioArray[i][j] += float(AcceptArray[i][j])/float(AttemptArray[i][j])
#  if nType > 1:
#    np.savetxt('NodalAcceptArray_'+str(iType)+'.dat',AcceptArray,fmt='%d')
#    np.savetxt('NodalAttemptArray_'+str(iType)+'.dat',AttemptArray,fmt='%d')
#    np.savetxt('NodalAcceptRatio_'+str(iType)+'.dat',AcceptRatioArray,fmt='%f')
#  else:
#    np.savetxt('NodalAcceptArray.dat',AcceptArray,fmt='%d')
#    np.savetxt('NodalAttemptArray.dat',AttemptArray,fmt='%d')
#    np.savetxt('NodalAcceptRatio.dat',AcceptRatioArray,fmt='%f')


#  norm = []
#  tot = 0
#  for i in range(0,5):
#    norm.append(sum(AttemptArray[i]))
#    tot += sum(AttemptArray[i])
#    AcceptRatioArray[i] = AcceptArray[i]/norm[i]
#  print AcceptRatioArray
#  print np.array(norm)/tot
#
#
#  TotColRatios,TotRowRatios = [],[]
#  for i in range(0,5):
#    TotColRatios.append(sum([x[i] for x in AcceptRatioArray]))
#    TotRowRatios.append(sum(AcceptRatioArray[i]))
#  NormColRatios = np.array(TotColRatios)
#  NormRowRatios = np.array(TotRowRatios)
#
#  print NormColRatios, NormRowRatios
#  print (NormColRatios/NormRowRatios)/sum(NormColRatios/NormRowRatios)

