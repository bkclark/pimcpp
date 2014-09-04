import sys
import numpy as np
import h5py as h5
import Stats

StartCut = int(sys.argv[1])
suffix = sys.argv[2]

ENames = ['Kinetic','Node','dUExt','dULong','dUShort','ShortRange','VShort','IlkkaLongRange','DavidLongRange','VLong','Total']
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
  f.close()


g = open('Energy-'+suffix+'.dat','w')
for EName in ENames:
  try:
    TotStats = Stats.UnweightedAvg(EStats[EName])
    print EName, TotStats
    g.write('%s %f %f %f %f\n'%(EName, TotStats[0], TotStats[1], TotStats[2], TotStats[3]))
  except:
    pass

tailNames = ['duLong_k0','duLong_r0','vLong_k0','vLong_r0']
fname = sys.argv[3:][-1]
f = h5.File(fname,'r')
for tailName in tailNames:
  try:
    tails = np.array(f['Observables/Energy/'+tailName])
    print tailName, tails[0]
    g.write('%s '%(tailName))
    for tail in tails[0]:
        g.write('%f '%(tail))
    g.write('\n')
  except:
    pass
f.close()

g.close()
