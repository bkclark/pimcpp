import sys
import numpy as np
import h5py as h5
import Stats

StartCut = 1

ENames = ['Kinetic','Node','Total','dUExt','dULong','dUShort']
EStats = {}
for fname in sys.argv[1:]:
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


for EName in ENames:
  try:
    TotStats = Stats.UnweightedAvg(EStats[EName])
    print EName, TotStats
  except:
    pass

