import sys
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

# Grab StartCut (NEquilibrate)
def getStartCut(f):
  if StartCutOverride:
    return 0
  else:
    StartCut = f['System/NEquilibrate'][0]
  return StartCut

# Grabs data from given Section
def FetchData(files,Section,StartCut,Col):
  AllData = []
  for filename in files:
    try:
      f = h5.File(filename,'r')
      data = f[Section]
      if StartCut < 0:
        StartCut = getStartCut(f)
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
      f.close()
    except:
      continue
  return AllData

# Grabs data from given Section
def FetchDataList(files,Base,Sections,StartCut,Col):
  AllData = []
  for i in range(0,len(Sections)):
    AllData.append([])
  for filename in files:
      try:
        f = h5.File(filename,'r')
      except:
        print filename, 'is not found'
        return AllData,Sections
      if StartCut < 0:
        StartCut = getStartCut(f)
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
      f.close()
  return AllData,Sections

# Grabs data from given Section
def FetchDataArray(files,Section,StartCut,nCol):
  AllData = []
  for i in range(0,nCol):
    AllData.append([])
  for filename in files:
    try:
      f = h5.File(filename,'r')
      if StartCut < 0:
        StartCut = getStartCut(f)
      data = f[Section]
      data = transpose(data)
      for i in range(0,nCol):
        if len(data[i]) > StartCut:
          AllData[i].append(data[i][StartCut:])
        elif len(data[i]) == 1:
          AllData[i].append(data[i][:])
        #else:
        #  print 'Warning: StartCut > len(data)!', filename, Section
      f.close()
    except:
      continue
  return AllData

def FetchStrings(files,SectionBase):
  f = h5.File(files[0],'r')
  strings = f[SectionBase]
  return strings[:]

def RunThroughSections(files,basename,prefix,outfile,Sections):
  data,Sections = FetchDataList(files,prefix,Sections,-1,0)
  g = open(basename+'/'+outfile,'w')
  i = 0
  for Section in Sections:
    try:
      stats = []
      for datum in data[i]:
        stat = Stats.stats(datum)
        if stat[0] < 1e100:
          stats.append(stat)
      [Mean,StdDev,Kappa,StdErr] = Stats.UnweightedAvg(stats)
      print Section, Mean, StdDev, Kappa, StdErr
      g.write('%s %f %f\n' % (Section, Mean, StdErr))
    except:
      print 'Error in Section: ', Section
      g.write('%s %f %f\n' % (Section, 0.0, 0.0))
    i += 1
  g.close()

# Moves
def ProcessMoves(files,basename):
  print '*** Processing Moves ***'
  Sections = ['BisectionBlock/AcceptRatio','RefSlice/AcceptRatio','ShiftMove/AcceptRatio']
  RunThroughSections(files,basename,'Moves','AvgMoves.dat',Sections)

# Energy
def ProcessEnergy(files,basename):
  print '*** Processing Energies ***'
  Sections = ['Energy Vals','Kinetic','Node','Residual Energy','Total','VLong','VShort','VTail Long Range','VTail Short Range','dULong','dUShort','dUNonlocal']
  RunThroughSections(files,basename,'Observables/Energy','AvgEnergies.dat',Sections)

# EnergySign
def ProcessEnergySign(files,basename):
  print '*** Processing Energy Signs ***'
  Sections = ['Kinetic','Node','Total','VLong','VShort','dULong','dUShort','dUNonlocal']
  RunThroughSections(files,basename,'Observables/EnergySign','AvgEnergySigns.dat',Sections)

# Sign
def ProcessSign(files,basename):
  print '*** Processing Signs ***'
  Sections = ['Total']
  RunThroughSections(files,basename,'Observables/Sign','AvgSigns.dat',Sections)

# Adjust Energy Values by Sign
def AdjustValuesBySign(files,basename):
  print '*** Adjusted Values ***'
  g = open(basename+'/AvgEnergySigns.dat','r')
  h = open(basename+'/AvgSigns.dat','r')
  i = open(basename+'/AvgAdjEnergySigns.dat','w')
  for line in h:
    sgn = float(line.split()[1])
    dsgn = float(line.split()[2])
  for line in g:
    line = line.split()
    x = float(line[1])
    dx = float(line[2])
    if (x != 0) and (sgn != 0):
      dx = abs(x/sgn)*sqrt((dx/x)**2 + (dsgn/sgn)**2)
      x = x/sgn
    print line[0], x, dx
    i.write('%s %f %f\n' % (line[0], x, dx))

# Time Analysis
def ProcessTimeAnalysis(files,basename):
  print '*** Processing Time Analysis ***'
  Sections = ['Action','Move','Observable']
  g = open(basename+'/TimeAnalysis.dat','w')
  for Section in Sections:
    try:
      names = FetchStrings(files,'Observables/TimeAnalysis/'+Section+'Names')
      data = FetchDataArray(files,'Observables/TimeAnalysis/'+Section+'Time',-1,len(names))
      for i in range(0,len(names)):
        stats = []
        for datum in data[i]:
          stats.append(Stats.stats(datum))
        [Mean,StdDev,Kappa,StdErr] = Stats.UnweightedAvg(stats)
        print names[i], Mean, StdDev, Kappa, StdErr
        g.write('%s %f %f %f %f\n' % (names[i], Mean, StdDev, Kappa, StdErr))
    except:
      print 'Error in Time Anaylsis,', Section
  g.close()

# Structure Factor
def ProcessStructureFactor(files,basename,nSpecies):
  print '*** Processing Structure Factor ***'
  nCombos = nSpecies*nSpecies - choose(nSpecies,2)
  for iSpecies in range(0,nCombos):
    if nSpecies == 1:
      speciesString = ''
    else:
      speciesString = '_'+str(iSpecies)
      print 'StructureFactor'+speciesString
    try:
      xs = FetchData(files,'Observables/StructureFactor'+speciesString+'/x',0,0)
      xs = xs[0].round(decimals=8)
      Species1 = FetchStrings(files,'Observables/StructureFactor'+speciesString+'/Species1')[0]
      Species2 = FetchStrings(files,'Observables/StructureFactor'+speciesString+'/Species2')[0]
      ys = FetchDataArray(files,'Observables/StructureFactor'+speciesString+'/y',-1,len(xs))
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
      g = open(basename+'/StructureFactor'+speciesString+'.dat','w')
      for i in range(0,len(ks)):
        print ks[i], Means[i], StdDevs[i], Kappas[i], StdErrs[i]
        g.write('%f %f %f\n' % (ks[i], Means[i], StdErrs[i]))
      g.close()
    except:
      print 'Error in Structure Factor', iSpecies

# Pair Correlation
def ProcessPairCorrelation(files,basename,nSpecies):
  print '*** Processing Pair Correlation ***'
  nCombos = nSpecies*nSpecies - choose(nSpecies,2)
  for iSpecies in range(0,nCombos):
    if nSpecies == 1:
      speciesString = ''
    else:
      speciesString = '_'+str(iSpecies)
      print 'PairCorrelation'+speciesString
    try:
      rs = FetchData(files,'Observables/PairCorrelation'+speciesString+'/x',0,0)
      rs = rs[0]
      Species1 = FetchStrings(files,'Observables/PairCorrelation'+speciesString+'/Species1')[0]
      Species2 = FetchStrings(files,'Observables/PairCorrelation'+speciesString+'/Species2')[0]
      grs = FetchDataArray(files,'Observables/PairCorrelation'+speciesString+'/y',-1,len(rs))
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
      g = open(basename+'/PairCorrelation'+speciesString+'.dat','w')
      for i in range(0,len(rs)):
        if rs[i] <= Box[0][0]/2:
          print rs[i], Means[i], StdDevs[i], Kappas[i], StdErrs[i]
          g.write('%f %f %f\n' % (rs[i], Means[i], StdErrs[i]))
      g.close()
    except:
      print 'Error in Pair Correlation'


# Get Base File Name
if (sys.argv[1]=='-f'):
  print "NOTE: First argument is -f"
  basename=sys.argv[2]
  listFile=file(basename)
  filesList=listFile.readlines()
  for counter in range(0,len(filesList)):
    filesList[counter]=filesList[counter][:-1]
  print "The file list is ",filesList
  basename=basename+"d"
else:
  try:
    fullpath = sys.argv[1].split('/')
    basename = fullpath[2]
    basedir = fullpath[0]+'/'+fullpath[1]+'/'
  except:
    fullpath
    basename = sys.argv[1]
    basedir = '../rawdata/'
  print basedir, basename
  splitBaseName=string.split(basename,'.')
  if splitBaseName[-1]=='h5':
    basename=string.join(splitBaseName[0:-2],'.')
print 'Basename:', basename

# Get Start Cut
try:
  StartCut = int(sys.argv[2])
  StartCutOverride = 1
except:
  StartCut = 0
  StartCutOverride = 1 #HACK
print 'Start Cut:', StartCut

# Get number of species
try:
  nSpecies = int(sys.argv[3])
except:
  nSpecies = 1

# Get max files
try:
  maxFiles = int(sys.argv[4])
except:
  maxFiles = 100000000

# Make Files List
files = []
clonestr = ''
i = 0
while 1:
  proc = 0
  while 1:
    filename = basedir+basename+clonestr+'.'+str(proc)+'.h5'
    if not os.path.exists(filename) or proc >= maxFiles:
      break
    files.append(filename)
    proc += 1
  clonestr = '.'+str(i)
  i += 1
  if proc == 0 or proc > maxFiles:
    break

# File List Info
nproc = len(files)
print 'Found', nproc, 'Files'
if nproc == 0:
  print 'No Files!'
  exit()

if not os.path.exists(basename):
  print 'Creating New Directory', basename
  os.makedirs(basename)


def CheckAcceptance(files,basename,nSpecies):
  if nSpecies == 1:
    speciesString = ''
  else:
    speciesString = '_0'
  Moves = ['BisectionBlock'+speciesString+'/AcceptRatio'] # Will need to be more general eventually
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
      if Mean == 0.0:
        return 0
    except:
      print 'Error in CheckAcceptance'
      return 0
    i += 1
  return 1

print '*** Checking for dead runs ***'
g = open(basename+'.dead','w')
count = 0
toremove = []
for f in files:
  if not CheckAcceptance([f],basename,nSpecies):
    print f, 'is bad'
    toremove.append(count)
    g.write(f+'\n')
  count += 1
toremove.reverse()
map(files.pop,toremove)

# Process Observables
ProcessMoves(files,basename)
ProcessEnergy(files,basename)
#ProcessSign(files,basename)
#ProcessEnergySign(files,basename)
#AdjustValuesBySign(files,basename)
ProcessTimeAnalysis(files,basename)
#ProcessStructureFactor(files,basename,nSpecies)
ProcessPairCorrelation(files,basename,nSpecies)
