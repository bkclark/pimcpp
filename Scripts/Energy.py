import sys
import numpy as np
import h5py as h5
import Stats

StartCut = int(sys.argv[1])
N = int(sys.argv[2])
suffix = sys.argv[3]

EStats = {}
for fname in sys.argv[4:]:
    f = h5.File(fname,'r')
    ENames = f['Observables/Energy'].keys()
    for EName in ENames:
        try:
            Es = f['Observables/Energy/'+EName][StartCut:]
            if EName in EStats:
                EStats[EName].append(Stats.stats(Es))
            else:
                EStats[EName] = [Stats.stats(Es)]
        except:
            pass

    f.flush()
    f.close()

g = open('Energy-'+suffix+'.dat','w')
perParticleE, perParticleEErr = 0., 0.
for EName in EStats.keys():
  try:
    TotStats = Stats.UnweightedAvg(EStats[EName])
    print EName, TotStats
    g.write('%s %f %f %f %f\n'%(EName, TotStats[0], TotStats[1], TotStats[2], TotStats[3]))
    if EName == 'Total':
        perParticleE += TotStats[0]/N
        perParticleEErr += TotStats[3]/N
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

    if tailName == 'duLong_k0' or tailName == 'duLong_r0':
        for tail in tails[0]:
            perParticleE += tail/N
  except:
    pass
f.close()

print 'PerParticleE', perParticleE, perParticleEErr
g.write('PerParticleE %f %f'%(perParticleE,perParticleEErr))

g.close()
