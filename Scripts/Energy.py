import sys
import numpy as np
import h5py as h5
import Stats

StartCut = int(sys.argv[1])
suffix = sys.argv[2]

ENames = ['Kinetic','Node','dUExt','dULong','dUShort','ShortRange','VShort','DavidLongRange','VLong','Total']
EStats = {}
for fname in sys.argv[3:]:
  f = h5.File(fname,'r')

  for EName in ENames:
    try:
      Es = np.array(f['Observables/Energy/'+EName][StartCut:])
      try:
        EStats[EName].append(Stats.stats(Es))
      except:
        EStats[EName] = [Stats.stats(Es)]
    except:
      pass

f = open('Energy-'+suffix+'.dat','w')
for EName in ENames:
  try:
    TotStats = Stats.UnweightedAvg(EStats[EName])
    print EName, TotStats
    f.write('%s %f %f %f %f\n'%(EName, TotStats[0], TotStats[1], TotStats[2], TotStats[3]))
  except:
    pass
f.close()

