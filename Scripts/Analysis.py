import sys
import getopt
sys.path.reverse()
sys.path.append('/g/g15/brown281/lib/python2.6/site-packages')
sys.path.reverse()
from numpy import *
from math import *
import os
import string
import h5py as h5
import Stats


def choose(n, k):
  if 0 <= k <= n:
    ntok = 1
    ktok = 1
    for t in xrange(1, min(k, n - k) + 1):
      ntok *= n
      ktok *= t
      n -= 1
    return ntok // ktok
  else:
    return 0


# Grabs data from given Section
def FetchData(files,Section,StartCut,Col):
  AllData = []
  count = 0
  for f in files:
    try:
      #f = h5.File(filename,'r')
      data = f[Section]
      if len(data.shape) == 1:
        if len(data) > StartCut:
          AllData.append(data[StartCut:])
        elif len(data) == 1:
          AllData.append(data[:])
        else:
          print 'Warning: StartCut > len(data)!'
      else:
        data = transpose(data)
        if len(data[Col]) > StartCut:
          AllData.append(data[Col][StartCut:])
        elif len(data[Col]) == 1:
          AllData.append(data[Col][:])
        else:
          print 'Warning: StartCut > len(data)!'
      #f.close()
    except:
      continue
    count += 1
  return AllData


# Grabs data from given Section
def FetchDataList(files,Base,Sections,StartCut,Col):
  AllData = []
  for i in range(0,len(Sections)):
    AllData.append([])
  count = 0
  for f in files:
    i = 0
    for Section in Sections:
      try:
        data = f[Base+'/'+Section]
        if len(data.shape) == 1:
          if len(data) > StartCut:
            AllData[i].append(data[StartCut:])
          elif len(data) == 1:
            AllData[i].append(data[:])
          else:
            print 'Warning: StartCut > len(data)!'
        else:
          data = transpose(data)
          if len(data[Col]) > StartCut:
            AllData[i].append(data[Col][StartCut:])
          elif len(data[Col]) == 1:
            AllData[i].append(data[Col][:])
          else:
            print 'Warning: StartCut > len(data)!'
      except:
        AllData[i].append(array([0.]))
        print 'Error in', Section
      i += 1
    count += 1
    #f.close()
  return AllData,Sections


# Grabs data from given Section
def FetchDataArray(files,Section,StartCut):
  AllData = []
  count = 0
  for f in files:
    try:
      #f = h5.File(filename,'r')
      data = f[Section]
      data = transpose(data)
      nCol = len(data)
      for i in range(0,nCol):
        if (len(AllData) != nCol):
          AllData.append([])
        if len(data[i]) > StartCut:
          AllData[i].append(data[i][StartCut:])
        elif len(data[i]) == 1:
          AllData[i].append(data[i][:])
        #else:
        #  print 'Warning: StartCut > len(data)!', filename, Section
      #f.close()
    except:
      continue
    count += 1
  return AllData


# Grabs data from given Section
def FetchDataArray2(files,Section,StartCut):
  AllData = []
  count = 0
  for f in files:
    try:
      data = f[Section]
      nCol = len(data)
      AllData.append(array(data))
    except:
      continue
    count += 1
  return AllData


def FetchStrings(files,SectionBase):
  #f = h5.File(files[0],'r')
  f = files[0]
  strings = f[SectionBase]
  return strings[:]


def RunThroughSections(files,basename,prefix,outfile,Sections,StartCut):
  data,Sections = FetchDataList(files,prefix,Sections,StartCut,0)
  g = open(basename+'/'+outfile,'w')
  i = 0
  for Section in Sections:
    stats = []
    count = 0
    for datum in data[i]:
      [mean,stddev,kappa,stderr] = Stats.stats(datum)
      if stderr != 0.:
        neff = (stddev/stderr)**2.
      else:
        neff = 1.
      stat = [mean,stderr,neff]
      #stat = [mean,stddev,kappa,stderr]
      if stat[0] < 1e50:
        stats.append(stat)
      count += 1
    [Mean,StdErr,Neff] = Stats.UnweightedReAvg(stats)
    #[Mean,StdDev,Kappa,StdErr] = Stats.UnweightedAvg(stats)
    g.write('%s %e %e %e\n' % (Section.replace(' ',''), Mean, StdErr, Neff))
    #g.write('%s %e %e\n' % (Section.replace(' ',''), Mean, StdErr))
    i += 1
  g.close()


def GetSubSections(files,SectionName):
  #f = h5.File(files[0],'r')
  f = files[0]
  Section = f[SectionName]
  Sections = Section.keys()
  if 'Description' in Sections:
    Sections.remove('Description')
  if 'Type' in Sections:
    Sections.remove('Type')
  return Sections


def ProcessMoves(files,basename,StartCut,clone):
  print '* Processing Moves'
  Sections = [x+'/AcceptRatio' for x in GetSubSections(files,'Moves')]
  RunThroughSections(files,basename,'Moves','AvgMoves'+clone+'.dat',Sections,StartCut)


def ProcessEnergy(files,basename,StartCut,clone):
  print '* Processing Energies'
  #Sections = ['Kinetic','Node','Residual Energy','Total','VLong','VShort','VTail Long Range','VTail Short Range','dULong','dUShort','dUNonlocal']
  Sections = GetSubSections(files,'Observables/Energy')
  RunThroughSections(files,basename,'Observables/Energy','AvgEnergies'+clone+'.dat',Sections,StartCut)


def ProcessSign(files,basename,StartCut,clone):
  print '* Processing Signs'
  Sections = ['Total']
  RunThroughSections(files,basename,'Observables/Sign','AvgSigns'+clone+'.dat',Sections,StartCut)


def ProcessImportanceWeight(files,basename,StartCut,clone):
  print '* Processing Importance Weight'
  Sections = ['Total']
  RunThroughSections(files,basename,'Observables/ImportanceWeight','AvgImportanceWeight'+clone+'.dat',Sections,StartCut)

BlockSize = 100

def ProcessSectorCounts(files,basename,StartCut,clone):
  Sections = [['CycleCount/SectorCount','AvgSectorCounts']]
  for [Section,OutPrefix] in Sections:
    print '* Processing', Section
    try:
      data = FetchDataArray2(files,'Observables/'+Section,StartCut)
      Data = {}
      nFiles = len(data)
      nSamples = 0
      for i in range(0,nFiles):
        for j in range(StartCut*BlockSize,len(data[i])):
          label = int(data[i][j][0])
          val = data[i][j][1]
          nSamples += val
          if label in Data:
            Data[label].append(val)
          else:
            Data[label] = [val]
      for label in Data.iterkeys():
        Data[label] = sum(array(Data[label]))/nSamples
      g = open(basename+'/'+OutPrefix+clone+'.dat','w')
      for label in sorted(Data.iterkeys()):
        [mean,sigma,kappa,err] = [Data[label],0.,1.,0.]
        #print label, mean, sigma, kappa ,err
        g.write('%i %e %e %i\n' % (label, mean, err, nSamples))
      g.close()
    except (KeyboardInterrupt, SystemExit):
      raise
    except:
      print 'Error in', Section


def ProcessPermEnergies(files,basename,StartCut,clone):
  Sections = [['Energy/PermEnergy','AvgPermEnergies']]
  for [Section,OutPrefix] in Sections:
    print '* Processing', Section
    try:
      data = FetchDataArray2(files,'Observables/'+Section,StartCut)
      Data = {}
      nFiles = len(data)
      for i in range(0,nFiles):
        nSamples = len(data[i])
        for j in range(StartCut*BlockSize,nSamples):
          label = int(data[i][j][0])
          val = data[i][j][1]
          N = data[i][j][2]
          err = data[i][j][3]
          if label in Data:
            Data[label].append([val,err,N])
          else:
            Data[label] = [[val,err,N]]
      g = open(basename+'/'+OutPrefix+clone+'.dat','w')
      for label in sorted(Data.iterkeys()):
        [mean,err,N] = Stats.UnweightedReAvg(Data[label])
        g.write('%i %e %e %i\n' % (label, mean, err, N))
      g.close()
    except:
      print 'Error in', Section


def ProcessCycleCounts(files,basename,StartCut,clone):
  Sections = [['CycleCount/y','AvgCycleCounts']]
  for [Section,OutPrefix] in Sections:
    print '* Processing', Section
    try:
      data = FetchDataArray(files,'Observables/'+Section,StartCut)
      stats = []
      for i in range(0,len(data)):
        stats.append([])
        for f in range(0,len(data[0])):
          stats[i].append(Stats.stats(data[i][f]))
        stats[i] = Stats.UnweightedAvg(stats[i])
      Means = [s[0] for s in stats]
      StdDevs = [s[1] for s in stats]
      Kappas = [s[2] for s in stats]
      StdErrs = [s[3] for s in stats]
      g = open(basename+'/'+OutPrefix+clone+'.dat','w')
      for i in range(0,len(stats)):
        #print i, Means[i], StdDevs[i], Kappas[i], StdErrs[i]
        g.write('%i %e %e\n' % (i, Means[i], StdErrs[i]))
      g.close()
    except (KeyboardInterrupt, SystemExit):
      raise
    except:
      print 'Error in', Section


def ProcessPerms(files,basename,StartCut,clone):
  ProcessCycleCounts(files,basename,StartCut,clone)
  ProcessPermEnergies(files,basename,StartCut,clone)
  ProcessSectorCounts(files,basename,StartCut,clone)


def ProcessCentroid(files,basename,StartCut,clone):
  print '* Processing Centroid'
  ProcessCentroidDist(files,basename,StartCut,clone)
  ProcessCentroidSpread(files,basename,StartCut,clone)


def ProcessCentroidDist(files,basename,StartCut,clone):
  Section = 'Centroid'
  try:
    rs = FetchData(files,'Observables/'+Section+'/x',0,0)
    rs = rs[0]
    grs = FetchDataArray(files,'Observables/'+Section+'/Spread',StartCut)
    stats = []
    for i in range(0,len(grs)):
      stats.append([])
      for f in range(0,len(grs[0])):
        stats[i].append(Stats.stats(grs[i][f]))
      stats[i] = Stats.UnweightedAvg(stats[i])
    Means = [s[0] for s in stats]
    StdDevs = [s[1] for s in stats]
    Kappas = [s[2] for s in stats]
    StdErrs = [s[3] for s in stats]
    g = open(basename+'/'+Section+'SpreadDist'+clone+'.dat','w')
    for i in range(0,len(rs)):
      g.write('%e %e %e\n' % (rs[i], Means[i], StdErrs[i]))
    g.close()
  except:
    print 'Error in', Section


def ProcessCentroidSpread(files,basename,StartCut,clone):
  Section = 'Centroid'
  try:
    vals = FetchDataArray2(files,'Observables/'+Section+'/SpreadVals',StartCut)
    vecs = FetchDataArray2(files,'Observables/'+Section+'/SpreadVecs',StartCut)
    N = len(vals[0][0]) # Particles
    D = len(vals[0][0][0]) # Dimensions
    ValStats,VecStats = [],[]
    for i in range(0,N): # particle
      ValStats.append([])
      VecStats.append([])
      for d in range(0,D): # dimension
        ValStats[i].append([])
        for f in range(0,len(vals)): # file
          for m in range(0,len(vals[f])): # measurement
            ValStats[i][d].append(vals[f][m][i][d])
        ValStats[i][d] = Stats.stats(ValStats[i][d])
        VecStats[i].append([])
        for d2 in range(0,D):
          VecStats[i][d].append([])
          for f in range(0,len(vecs)): # file
            for m in range(0,len(vecs[f])): # measurement
              VecStats[i][d][d2].append(vecs[f][m][i][d][d2])
          VecStats[i][d][d2] = Stats.stats(VecStats[i][d][d2])
    g = open(basename+'/'+Section+'Spread'+clone+'.dat','w')
    for i in range(0,N):
      g.write('%i ' % (i))
      for d in range(0,D):
        g.write('%e %e ' % (ValStats[i][d][0], ValStats[i][d][3]))
      for d1 in range(0,D):
        for d2 in range(0,D):
          g.write('%e %e ' % (VecStats[i][d1][d2][0], VecStats[i][d1][d2][3]))
      g.write('\n')
    g.close()
  except:
    print 'Error in', Section


def AdjustValuesByImportanceWeight2(basename,Sections):
  print '* Adjusting Values By Importance Weight'
  try:
    f = open(basename+'/AvgImportanceWeight.dat','r')
  except:
    f = open(basename+'/AvgSigns.dat','r')
  for line in f:
    sgn = float(line.split()[-2])
    dsgn = float(line.split()[-1])
  f.close()
  for Section in Sections:
    try:
      g = open(basename+'/'+Section+'.dat','r')
      h = open(basename+'/Adj'+Section+'.dat','w')
      for line in g:
        line = line.split()
        label = line[0]
        x = float(line[1])
        dx = float(line[2])

        if (sgn != 0) and (label != 'VTailLongRange') and (label != 'VTailShortRange'):
          x = x/sgn
          if (x != 0) and (dx != float('inf')) and (dsgn != float('inf')):
            dx = abs(x/sgn)*sqrt((dx/x)**2 + (dsgn/sgn)**2)
        #print label, x, dx
        h.write('%s %e %e\n' % (label, x, dx))
      g.close()
      h.close()
    except:
      print 'Error in', Section


def AdjustValuesByImportanceWeight(basename,Sections,clone):
  print '* Adjusting Values By Importance Weight'
  try:
    f = open(basename+'/AvgImportanceWeight'+clone+'.dat','r')
  except:
    f = open(basename+'/AvgSigns'+clone+'.dat','r')
  for line in f:
    sgn = float(line.split()[1])
    dsgn = float(line.split()[2])
  f.close()
  for Section in Sections:
    try:
      g = open(basename+'/'+Section+clone+'.dat','r')
      h = open(basename+'/Adj'+Section+clone+'.dat','w')
      for line in g:
        line = line.split()
        label = line[0]
        x = float(line[1])
        dx = float(line[2])
        if (sgn != 0) and (label != 'VTailLongRange') and (label != 'VTailShortRange'):
          x = x/sgn
          if (x != 0) and (dx != float('inf')) and (dsgn != float('inf')):
            dx = abs(x/sgn)*sqrt((dx/x)**2 + (dsgn/sgn)**2)
        #print label, x, dx
        h.write('%s %e %e\n' % (label, x, dx))
      g.close()
      h.close()
    except:
      print 'Error in', Section


def ProcessTimeAnalysis(files,basename,StartCut,clone):
  print '* Processing Time Analysis'
  Sections = ['Action','Move','Observable']
  g = open(basename+'/TimeAnalysis'+clone+'.dat','w')
  for Section in Sections:
    try:
      names = FetchStrings(files,'Observables/TimeAnalysis/'+Section+'Names')
      data = FetchDataArray(files,'Observables/TimeAnalysis/'+Section+'Time',StartCut)
      for i in range(0,len(names)):
        stats = []
        for datum in data[i]:
          stats.append(Stats.stats(datum))
        [Mean,StdDev,Kappa,StdErr] = Stats.UnweightedAvg(stats)
        #print names[i], Mean, StdDev, Kappa, StdErr
        g.write('%s %e %e\n' % (names[i], Mean, StdErr))
    except (KeyboardInterrupt, SystemExit):
      raise
    except:
      print 'Error in Time Anaylsis,', Section
  g.close()


def ProcessStructureFactor(files,basename,StartCut,Section,clone):
  print '* Processing', Section
  try:
    xs = FetchData(files,'Observables/'+Section+'/x',0,0)
    xs = xs[0].round(decimals=8)
    Species1 = FetchStrings(files,'Observables/'+Section+'/Species1')[0]
    Species2 = FetchStrings(files,'Observables/'+Section+'/Species2')[0]
    ys = FetchDataArray(files,'Observables/'+Section+'/y',StartCut)
    ks = unique(xs) # Pick out unique k values
    Sk = []
    for i in range(0,len(ks)):
      Sk.append([])
      for f in range(0,len(ys[0])):
        Sk[i].append(array([]))
    for i in range(0,len(ks)):
      for j in range(0,len(xs)):
        if ks[i] == xs[j]:
          for f in range(0,len(ys[j])):
            try:
              Sk[i][f] = append(Sk[i][f],ys[j][f]) # Combine data for unique k values
            except:
              print i, j, f
              exit()
    stats = []
    for i in range(0,len(Sk)):
      stats.append([])
      for f in range(0,len(ys[0])):
        stats[i].append(Stats.stats(Sk[i][f]))
      stats[i] = Stats.UnweightedAvg(stats[i])
    Means = [s[0] for s in stats]
    StdDevs = [s[1] for s in stats]
    Kappas = [s[2] for s in stats]
    StdErrs = [s[3] for s in stats]
    #StdErrs /= Means[argmax(ks)]
    #Means /= Means[argmax(ks)]
    Everything = zip(ks,Means,StdDevs,Kappas,StdErrs)
    Everything.sort()
    ks,Means,StdDevs,Kappas,StdErrs = zip(*Everything)
    g = open(basename+'/'+Section+clone+'.dat','w')
    for i in range(0,len(ks)):
      #print ks[i], Means[i], StdDevs[i], Kappas[i], StdErrs[i]
      g.write('%e %e %e\n' % (ks[i], Means[i], StdErrs[i]))
    g.close()
  except:
    print 'Error in', Section


def ProcessPairCorrelation(files,basename,StartCut,Section,clone,y='y'):
  print '* Processing', Section
  try:
    rs = FetchData(files,'Observables/'+Section+'/x',0,0)
    rs = rs[0]
    Species1 = FetchStrings(files,'Observables/'+Section+'/Species1')[0]
    Species2 = FetchStrings(files,'Observables/'+Section+'/Species2')[0]
    grs = FetchDataArray(files,'Observables/'+Section+'/'+y,StartCut)
    Box = FetchData(files,'System/Box',0,0)
    stats = []
    for i in range(0,len(grs)):
      stats.append([])
      for f in range(0,len(grs[0])):
        stats[i].append(Stats.stats(grs[i][f]))
      stats[i] = Stats.UnweightedAvg(stats[i])
    Means = [s[0] for s in stats]
    StdDevs = [s[1] for s in stats]
    Kappas = [s[2] for s in stats]
    StdErrs = [s[3] for s in stats]
    if y=='y':
      g = open(basename+'/'+Section+clone+'.dat','w')
    else:
      g = open(basename+'/'+Section+'Neg'+clone+'.dat','w')
    for i in range(0,len(rs)):
      if rs[i] <= Box[0][0]/2:
        #print rs[i], Means[i], StdDevs[i], Kappas[i], StdErrs[i]
        g.write('%e %e %e\n' % (rs[i], Means[i], StdErrs[i]))
    g.close()
  except:
    print 'Error in', Section


def ProcessNofr(files,basename,StartCut,Section,clone):
  print '* Processing', Section
  try:
    rs = FetchData(files,'Observables/'+Section+'/x',0,0)
    rs = rs[0]
    nrs = FetchDataArray(files,'Observables/'+Section+'/y',StartCut)
    Box = FetchData(files,'System/Box',0,0)
    stats = []
    for i in range(0,len(nrs)):
      stats.append([])
      for f in range(0,len(nrs[0])):
        stats[i].append(Stats.stats(nrs[i][f]))
      stats[i] = Stats.UnweightedAvg(stats[i])
    Means = [s[0] for s in stats]
    StdDevs = [s[1] for s in stats]
    Kappas = [s[2] for s in stats]
    StdErrs = [s[3] for s in stats]
    g = open(basename+'/'+Section+clone+'.dat','w')
    for i in range(0,len(rs)):
      #print rs[i], Means[i], StdDevs[i], Kappas[i], StdErrs[i]
      g.write('%e %e %e\n' % (rs[i], Means[i], StdErrs[i]))
    g.close()
  except:
    print 'Error in', Section


def CheckAcceptance(files,baseame,Section):
  Moves = [Section+'/AcceptRatio']
  data,Moves = FetchDataList(files,'Moves',Moves,0,0)
  i = 0
  for Move in Moves:
    try:
      stats = []
      for datum in data[i]:
        stat = Stats.stats(datum)
        if stat[0] < 1e100:
          stats.append(stat)
      [Mean,StdDev,Kappa,StdErr] = Stats.UnweightedAvg(stats)
    except:
      print 'Error in CheckAcceptance'
      return 0
    if Mean == 0.0:
      return 0
    i += 1
  return 1


def getnumber(s):
  if s.isdigit():
    return int(s)
  else:
    try:
      s = float(s)
      return s
    except ValueError:
      return s


def CombineData(basename,Sections,clones):
  for Section in Sections:
    hs = []
    for clone in clones:
      filename = basename+'/'+Section+clone+'.dat'
      try:
        hs.append(open(filename,'r'))
      except:
        print 'ERROR: DID NOT FIND', filename
    data = {}
    Ns = []
    for h in hs:
      i = 0
      N = 1
      for line in h:
        line = line.split()
        label = getnumber(line[0])
        y = float(line[1])
        dy = float(line[2])
        if len(line) == 4:
          N = float(line[3])
          if label in data:
            data[label].append([y,dy,N,i])
          else:
            data[label] = [[y,dy,N,i]]
        else:
          if label in data:
            data[label].append([y,dy])
          else:
            data[label] = [[y,dy]]
      i += 1
      Ns.append(N)
    totN = sum(Ns)
    for h in hs:
      h.close()
    g = open(basename+'/'+Section+'.dat','w')
    for label in sorted(data.iterkeys()):
      if len(data[label][0]) == 4: # Since we have N, we can do something more precise
        if Section == 'AvgSectorCounts':
          totX = 0.
          totX2 = 0.
          for x in data[label]:
            totX += x[0]
            totX2 += x[0]*x[0]
          N = len(clones)
          avgX = totX/N
          avgX2 = totX2/N
          dX = sqrt(abs(avgX2-avgX*avgX))/N
          y,dy,N = avgX,dX,N
        else:
          y,dy,N = Stats.UnweightedReAvg(data[label])
      else:
        y,dy = Stats.UnweightedAvg(data[label])
        N = len(data[label])
      g.write('%s %e %e %i\n' % (label, y, dy, N))
    g.close()


from optparse import OptionParser

def cb(option, opt_str, value, parser):
  args=[]
  for arg in parser.rargs:
    if arg[0] != "-":
      args.append(arg)
    else:
      del parser.rargs[:len(args)]
      break
  if getattr(parser.values, option.dest):
    args.extend(getattr(parser.values, option.dest))
  setattr(parser.values, option.dest, args)


def main():

  # Options
  usage = "usage: python %prog [options] ../path/to/fileprefix"
  parser = OptionParser(usage=usage)
  parser.add_option("-f", "--filelist", action="store", type="string", dest="basename", help="Use a file list")
  parser.add_option("-x", "--startcut", action="callback", callback=cb, dest="StartCuts", default=[], help="Starting data point [default: %default]")
  parser.add_option("-m", "--maxfiles", action="store", type="int", dest="maxFiles", default=1000000000, help="Maximum number of files [default: %default]")
  parser.add_option("-s", "--tracksign", action="store_true", dest="TrackSign", default=False, help="Adjust by average value of the sign [default: %default]")
  parser.add_option("-c", "--clones", action="callback", callback=cb, dest="clones", default=[], help="Use data from previous runs (give clone numbers) [default: %default]")
  parser.add_option("-d", "--checkdead", action="store_true", dest="CheckDead", default=False, help="Preemptively check for runs with zero acceptance ratio [default: %default]")
  parser.add_option("-n", "--noanalysis", action="store_false", dest="DoAnalysis", default=True, help="Don't rerun the analysis [default: %default]")
  parser.add_option("-u", "--usedclones", action="store_true", dest="UsedClones", default=False, help="Use UsedClones.dat file [default: %default]")
  parser.add_option("-g", "--genusedclones", action="store_true", dest="GenUsedClones", default=False, help="Generate UsedClones.dat file [default: %default]")
  parser.add_option("-a", "--appendusedclones", action="store_true", dest="AppendUsedClones", default=False, help="Append UsedClones.dat file [default: %default]")
  parser.add_option("-e", "--everyclone", action="store_true", dest="EveryClone", default=False, help="Do analysis for every clone [default: %default]")
  parser.add_option("-r", "--excluderun", action="callback", callback=cb, dest="ExcludedRuns", default=[], help="Exclude specific runs (e.g. the first one, 0) [default: %default]")
  parser.add_option("-q", "--startingrun", action="store", type="int", dest="StartingRun", default=0, help="Start at specific run (e.g. the first one, 0) [default: %default]")
  (opts, args) = parser.parse_args()
  print opts,args


  # Get Base File Name
  UseFileList = False
  SingleFile = False
  if opts.basename != None:
    listFile = file(basename)
    filesList = listFile.readlines()
    for counter in range(0,len(filesList)):
      filesList[counter] = filesList[counter][:-1]
    print "Using file list:",filesList
    basename=basename+"d"
    UseFileList = True
  if not UseFileList:
    fullpath = args[0].split('/')
    if len(fullpath) > 2:
      basename = fullpath[2]
      basedir = fullpath[0]+'/'+fullpath[1]+'/'
    else:
      basename = fullpath[-1]
      basedir = '../rawdata/'
    splitBaseName=string.split(basename,'.')
    if splitBaseName[-1]=='h5':
      basename=string.join(splitBaseName[0:-2],'.')
      SingleFile = True

  #basedir = '../rawdata/'
  #opts.GenUsedClones = True
  #opts.AppendUsedClones = True

  i = 0
  while i < opts.StartingRun:
    opts.ExcludedRuns.append(str(i))
    i += 1
  print 'Excluded Runs:', opts.ExcludedRuns

  # Use all clones
  if opts.EveryClone:
    clone = 0
    run = 0
    while str(run) in opts.ExcludedRuns:
      run += 1
    while 1:
      filename = basedir+basename+'.'+str(run)+'.'+str(clone)+'.h5'
      if not os.path.exists(filename):
        break
      opts.clones.append(str(clone))
      clone += 1

  # Adjust clones and StartCut lists
  if opts.UsedClones:
    try:
      f = open(basename+'/UsedClones.dat','r')
      g = []
      for line in f:
        g.append(line.split())
      opts.clones = g[0]
      if len(g) > 1:
        opts.StartCuts = g[1]
      else:
        opts.StartCuts = []
    except:
      print basename, ' has no UsedClones.dat !'
      sys.exit()
  if len(opts.clones) == 0:
    opts.clones = ['-1']
  if len(opts.StartCuts) == 0:
    opts.StartCuts = ['0' for clone in opts.clones]
  elif len(opts.StartCuts) == 1:
    opts.StartCuts = [opts.StartCuts[0] for clone in opts.clones]
  if len(opts.clones) == len(opts.StartCuts) - 1:
    opts.clones.append('-1')
  StartCuts = [int(StartCut) for StartCut in opts.StartCuts]
  print 'Clones:', opts.clones
  print 'StartCuts:', StartCuts
  #opts.clones = [opts.clones[1]]
  #StartCuts = [StartCuts[1]]

  # Make Files List
  runs = []
  if opts.DoAnalysis:
    if SingleFile:
      files = [[string.join(fullpath,'/')]]
      cloneCount = 1
      clones = ['-1']
    else:
      files = [[]]
      cloneStr = ''
      cloneCount = 0
      clones = []
      for clone in opts.clones:
        if clone == '-1':
          cloneStr = ''
        else:
          cloneStr = '.'+clone
        run = 0
        badCount = 0
        addRun = True
        while 1:
          if not str(run) in opts.ExcludedRuns:
            filename = basedir+basename+'.'+str(run)+cloneStr+'.h5'
            if not os.path.exists(filename) or run >= opts.maxFiles:
              badCount += 1
              if badCount > 2:
                badCount = 0
                break
            else:
              badCount = 0
              files[cloneCount].append(filename)
              if addRun:
                runs.append(run)
                addRun = False
          run += 1
          if run > opts.maxFiles:
            break
        if run == 0:
          break
        files.append([])
        cloneCount += 1
        clones.append(clone)
  else:
    clones = opts.clones
    cloneCount = len(clones)
    run = 0
    while 1:
      if not str(run) in opts.ExcludedRuns:
        break
      else:
        run += 1
    files = [basedir+basename+'.'+str(run)+'.'+clones[0]+'.h5']
    print files

  # Create directory
  if not os.path.exists(basename):
    print 'Creating New Directory', basename
    os.makedirs(basename)


  # Run Through Clones
  AllSections = []
  AdjSections = []
  cloneStrs = ['_'+clone for clone in clones]

  if not opts.DoAnalysis:
    # Only need 1 file to get Sections
    try:
      Sections = GetSubSections([h5.File(files[0],'r')],'Observables')
    except:
      try:
        f = open(basename+'/Sections.dat','r')
        Sections = f.readline().split()
        f.close()
      except:
        dirfiles = os.listdir(basename)
        for f in dirfiles:
          if "TimeAnalysis" in f:
            Sections.append("TimeAnalysis")
          if "AvgEnergies" in f:
            Sections.append("Energy")
          if "AvgMoves" in f:
            Sections.append("Moves")
          if "AvgCycleCounts" in f:
            Sections.append("CycleCount")
          if "Sign" in f:
            Sections.append("Sign")
          if "ImportanceWeight" in f:
            Sections.append("ImportanceWeight")

  for i in range(0,cloneCount):
    print '*** Working on clone', clones[i], '***'

    h5files = []
    if opts.DoAnalysis:
      Sections = []

      # Open h5 files
      for f in files[i]:
        try:
          h5files.append(h5.File(f,'r'))
        except:
          print 'Error with:',f
          #sys.exit()

      # Check for dead runs
      if opts.CheckDead:
        print '* Checking for dead runs *'
        count = 0
        toremove = []
        Sections = GetSubSections(h5files,'Moves')
        for Section in Sections:
          if "BisectionBlock" in Section:
            break
        for f in h5files:
          if not CheckAcceptance([f],basename,Section):
            print f, 'is bad'
            toremove.append(count)
          count += 1
        toremove.reverse()
        map(h5files.pop,toremove)

      # File List Info
      nrun = len(h5files)
      print 'Found', nrun, 'files for clone', clones[i]
      if nrun == 0:
        print 'No Files!'
        continue
        #exit()
      else:
        Sections = GetSubSections(h5files,'Observables')


    # Process Moves
    if opts.DoAnalysis:
      ProcessMoves(h5files,basename,StartCuts[i],cloneStrs[i])
    AllSections.append('AvgMoves')

    # Process Observables
    f = open(basename+'/Sections.dat','w')
    for Section in Sections:
      f.write(Section+' ')
      if "TimeAnalysis" in Section:
        if opts.DoAnalysis:
          ProcessTimeAnalysis(h5files,basename,StartCuts[i],cloneStrs[i])
        AllSections.append(Section)
      elif "Energy" in Section:
        if opts.DoAnalysis:
          ProcessEnergy(h5files,basename,StartCuts[i],cloneStrs[i])
        AdjSections.append('AvgEnergies')
        AllSections.append('AvgEnergies')
      elif "Centroid" in Section:
        if opts.DoAnalysis:
          ProcessCentroid(h5files,basename,StartCuts[i],cloneStrs[i])
        AllSections.append('CentroidSpread')
        #AllSections.append('CentroidSpreadVals')
      elif "CycleCount" in Section:
        if opts.DoAnalysis:
          ProcessPerms(h5files,basename,StartCuts[i],cloneStrs[i])
        AllSections.append('AvgCycleCounts')
        AllSections.append('AvgSectorCounts')
        AllSections.append('AvgPermEnergies')
      elif "Sign" in Section:
        if opts.DoAnalysis:
          ProcessSign(h5files,basename,StartCuts[i],cloneStrs[i])
        AllSections.append('AvgSigns')
      elif "ImportanceWeight" in Section:
        if opts.DoAnalysis:
          ProcessImportanceWeight(h5files,basename,StartCuts[i],cloneStrs[i])
        AllSections.append('AvgImportanceWeight')
      elif "PairCorrelation" in Section:
        if opts.DoAnalysis:
          ProcessPairCorrelation(h5files,basename,StartCuts[i],Section,cloneStrs[i])
        AdjSections.append(Section)
        AllSections.append(Section)
      elif "nofr" in Section:
        if opts.DoAnalysis:
          ProcessNofr(h5files,basename,StartCuts[i],Section,cloneStrs[i])
        AdjSections.append(Section)
        AllSections.append(Section)
      elif "StructureFactor" in Section:
        if opts.DoAnalysis:
          ProcessStructureFactor(h5files,basename,StartCuts[i],Section,cloneStrs[i])
        AdjSections.append(Section)
        AllSections.append(Section)
    f.close()

    # Adjust by sign
    if opts.TrackSign:
      AdjustValuesByImportanceWeight(basename,AdjSections,cloneStrs[i])

    for file in h5files:
      file.close()

  # Combine Clones
  print '*** Combining Clone Data ***'
  if opts.AppendUsedClones:
    oldClones = []
    oldStartCuts = []
    try:
      f = open(basename+'/UsedClones.dat','r')
      g = []
      for line in f:
        g.append(line.split())
      oldClones = g[0]
      if len(g) > 1:
        oldStartCuts = g[1]
      else:
        oldStartCuts = []
    except:
      print basename, ' has no UsedClones.dat !'
      sys.exit()
    oldCloneStrs = ['_'+clone for clone in oldClones]
    cloneStrs = cloneStrs + oldCloneStrs
    clones = clones + oldClones
    StartCuts = StartCuts + oldStartCuts
    tmpClones,tmpCloneStrs,tmpStartCuts = [],[],[]
    for i in range(0,len(clones)):
      if not clones[i] in tmpClones:
        tmpClones.append(clones[i])
        tmpCloneStrs.append(cloneStrs[i])
        tmpStartCuts.append(StartCuts[i])
    clones = tmpClones
    StartCuts = tmpStartCuts
    cloneStrs = tmpCloneStrs
  print 'clones:',clones
  print 'StartCuts:',StartCuts
  CombineData(basename,list(set(AllSections)),cloneStrs)

  # Adjust by sign
  if opts.TrackSign:
    AdjustValuesByImportanceWeight2(basename,AdjSections)

  # Generate UsedClones file
  if opts.GenUsedClones:
    g = open(basename+'/UsedClones.dat','w')
    for clone in clones:
      g.write(clone+' ')
    g.write('\n')
    for StartCut in StartCuts:
      g.write(str(StartCut)+' ')
    g.write('\n')
    for run in runs:
      g.write(str(run)+' ')
    g.close()


if __name__ == "__main__":
  sys.exit(main())
