import MPIClass
mpi = MPIClass.MPI(False)
import sys, os
import getopt, string
from math import sqrt
from numpy import unique, transpose
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
                if len(data.shape) == 2 and len(data[0]) == 1: # Kind of hacky for VTailShortRange
                    AllData[i].append(data[0])
                if len(data.shape) == 1:
                    if len(data) > StartCut:
                        AllData[i].append(data[StartCut:])
                    elif len(data) == 1:
                        AllData[i].append(data[:])
                    else:
                        print 'Warning: StartCut > len(data)!'
            except:
                AllData[i].append([0.])
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
                #    print 'Warning: StartCut > len(data)!', filename, Section
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
            data = f[Section][StartCut:]
            nCol = len(data)
            AllData.append(data)
        except:
            continue
        count += 1
    return AllData

# Grabs data from given Section
def FetchDataArrayChunk(files,Section,StartCut,ChunkNum,ChunkSize):
    AllData = []
    count = 0
    for f in files:
        try:
            data = f[Section][(ChunkSize*ChunkNum)+StartCut:ChunkSize*(ChunkNum+1)]
            nCol = len(data)
            AllData.append(data)
        except:
            pass
        count += 1
    return AllData

def FetchStrings(files,SectionBase):
    #f = h5.File(files[0],'r')
    f = files[0]
    strings = f[SectionBase]
    return strings[:]

def ProcessCentroid(files,basename,StartCut,clone):
    print '* Processing Centroid'
    ProcessCentroidDist(files,basename,StartCut,clone)
    ProcessCentroidSpread(files,basename,StartCut,clone)
    ProcessCentroidVals(files,basename,StartCut,clone)

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
            for d in range(0,D):
                g.write('%i%i %e %e \n' % (i, d, ValStats[i][d][0], ValStats[i][d][3]))
            for d1 in range(0,D):
                for d2 in range(0,D):
                    g.write('%i%i%i %e %e \n' % (i, d1, d2, VecStats[i][d1][d2][0], VecStats[i][d1][d2][3]))
        g.close()
    except:
        print 'Error in', Section

def ProcessCentroidVals(files,basename,StartCut,clone):
    Section = 'Centroid'
    try:
        xVals = FetchDataArray2(files,'Observables/'+Section+'/xVals',StartCut)
        yzVals = FetchDataArray2(files,'Observables/'+Section+'/yzVals',StartCut)
        N = len(xVals[0][0][0]) # Particles
        xValStats,yzValStats = [],[]
        for i in range(0,N): # particle
            xValStats.append([])
            yzValStats.append([])
            for d in range(0,2): # in between and outside
                xValStats[i].append([])
                for f in range(0,len(xVals)): # file
                    for m in range(0,len(xVals[f])): # measurement
                        xValStats[i][d].append(xVals[f][m][i][d])
                xValStats[i][d] = Stats.stats(xValStats[i][d])
            for f in range(0,len(yzVals)): # file
                for m in range(0,len(yzVals[f])): # measurement
                    yzValStats[i].append(yzVals[f][m][i])
            yzValStats[i] = Stats.stats(yzValStats[i])
        g = open(basename+'/'+Section+'xVals'+clone+'.dat','w')
        for i in range(0,N):
            for d in range(0,2):
                g.write('%i%i %e %e \n' % (i, d, xValStats[i][d][0], xValStats[i][d][3]))
        g.close()
        g = open(basename+'/'+Section+'yzVals'+clone+'.dat','w')
        for i in range(0,N):
            g.write('%i %e %e \n' % (i, yzValStats[i][0], yzValStats[i][3]))
        g.close()
    except:
        print 'Error in', Section

def GetSubSections(files,SectionName):
    #f = h5.File(files[0],'r')
    f = files[0]
    try:
        Section = f[SectionName]
        Sections = Section.keys()
        if 'Description' in Sections:
            Sections.remove('Description')
        if 'Type' in Sections:
            Sections.remove('Type')
        return Sections
    except:
        sys.stderr.write('Unable to open '+SectionName+'\n')
        return []

def GetObservable(basename,startCuts,section,cloneStrs):
    sectionInfo = [basename,startCuts,section,cloneStrs]
    if "TimeAnalysis" in section:
        return TimeAnalysis(*sectionInfo), 0
    if "SectorEnergy" in section:
        return SectorEnergy(*sectionInfo), 0
    elif "Energy" in section:
        return Energy(*sectionInfo), 1
    #elif "Centroid" in section:
    #    ProcessCentroid(h5files,basename,StartCuts[i],cloneStrs[i])
    elif "SectorCount" in section:
        return SectorCount(*sectionInfo), 0
    elif "CycleCount" in section:
        return CycleCount(*sectionInfo), 0
    elif "Sign" in section:
        return Sign(*sectionInfo), 0
    elif "ImportanceWeight" in section:
        return ImportanceWeight(*sectionInfo), 0
    elif "PairCorrelation" in section:
        return PairCorrelation(*sectionInfo), 1
    #elif "nofr" in section:
    #    return NofR(*sectionInfo), 1
    elif "StructureFactor" in section:
        return StructureFactor(*sectionInfo), 1
    else:
        return Observable(*sectionInfo), 0

class Observable:
    def __init__(self, basename, startCut, section, cloneStr):
        self.basename = basename
        self.startCut = startCut
        self.section = section
        self.cloneStr = cloneStr

    def Process(self, files):
        try:
            print 'Processing', self.section
            xs, ys = self.GetData(files)
            if (type(xs) is int and xs == 0) and (type(ys) is int and ys == 0):
                return 0
            stats = self.DoStatistics(xs,ys)
            self.WriteToFile(xs,stats)
            return 1
        except:
            sys.stderr.write('Error while processing '+self.section+'\n')
            return 0

    def GetData(self, files):
        print 'GetData not defined for', self.section
        return 0, 0

    def DoStatistics(self, xs, ys):
        stats = []
        for (x,y) in zip(xs,ys):
            statsx = []
            for yi in y:
                statsx.append(Stats.stats(yi))
            stats.append(Stats.UnweightedAvg(statsx))
        return stats

    def WriteToFile(self, xs, stats):
        g = open(self.basename+'/'+self.section+self.cloneStr+'.dat','w')
        for (x,stat) in zip(xs,stats):
            g.write(''.join(str(x).split())+' '+' '.join([str(x) for x in stat])+'\n')
        g.close()

def GetMove(basename,startCuts,section,cloneStrs):
    sectionInfo = [basename,startCuts,section,cloneStrs]
    if "BisectionBlock" in section:
        return BisectionBlock(*sectionInfo), 0
    elif "RefSliceShift" in section:
        return RefSliceShift(*sectionInfo), 0
    elif "RefSlice" in section:
        return RefSlice(*sectionInfo), 0
    else:
        return Observable(*sectionInfo), 0

class BisectionBlock(Observable):
    def GetData(self, files):
        Sections = ['AcceptRatio']
        data,Sections = FetchDataList(files,'Moves/BisectionBlock',Sections,self.startCut,0)
        return Sections, data

class RefSlice(Observable):
    def GetData(self, files):
        Sections = ['AcceptRatio']
        data,Sections = FetchDataList(files,'Moves/RefSlice',Sections,self.startCut,0)
        return Sections, data

class RefSliceShift(Observable):
    def GetData(self, files):
        Sections = ['AcceptRatio']
        data,Sections = FetchDataList(files,'Moves/RefSliceShift',Sections,self.startCut,0)
        return Sections, data

class Energy(Observable):
    def GetData(self, files):
        Sections = GetSubSections(files,'Observables/Energy')
        if 'SectorEnergy' in Sections:
            Sections.pop(Sections.index('SectorEnergy'))
        if 'PermEnergy' in Sections:
            Sections.pop(Sections.index('PermEnergy'))
        data,Sections = FetchDataList(files,'Observables/Energy',Sections,self.startCut,0)
        return Sections, data

class Sign(Observable):
    def GetData(self, files):
        Sections = ['Total']
        data,Sections = FetchDataList(files,'Observables/Sign',Sections,self.startCut,0)
        return Sections, data

class ImportanceWeight(Observable):
    def GetData(self, files):
        Sections = ['Total']
        data,Sections = FetchDataList(files,'Observables/ImportanceWeight',Sections,self.startCut,0)
        return Sections, data

class TimeAnalysis(Observable):
    def GetData(self, files):
        Sections = ['Action','Move','Observable']
        labels, vals = [],[]
        for Section in Sections:
            names = FetchStrings(files,'Observables/TimeAnalysis/'+Section+'Names')
            data = FetchDataArray(files,'Observables/TimeAnalysis/'+Section+'Time',self.startCut)
            for (name,datum) in zip(names,data):
                labels.append(name)
                vals.append(datum)
        return labels, vals

class StructureFactor(Observable):
    def GetData(self, files):
        xs = FetchData(files,'Observables/'+self.section+'/x',0,0)
        xs = xs[0].round(decimals=8)
        Species1 = FetchStrings(files,'Observables/'+self.section+'/Species1')[0]
        Species2 = FetchStrings(files,'Observables/'+self.section+'/Species2')[0]
        ys = FetchDataArray(files,'Observables/'+self.section+'/y',self.startCut)
        ks = unique(xs) # Pick out unique k values
        Sk = []
        for i in range(0,len(ks)):
            Sk.append([])
            for f in range(0,len(ys[0])):
                Sk[i].append([])
        for i in range(0,len(ks)):
            for j in range(0,len(xs)):
                if ks[i] == xs[j]:
                    for f in range(0,len(ys[j])):
                        try:
                            Sk[i][f] = Sk[i][f] + list(ys[j][f]) # Combine data for unique k values
                        except:
                            print i, j, f
                            sys.exit()
        return ks, Sk

class PairCorrelation(Observable):
    def GetData(self, files):
        rs = FetchData(files,'Observables/'+self.section+'/x',0,0)
        rs = rs[0]
        Species1 = FetchStrings(files,'Observables/'+self.section+'/Species1')[0]
        Species2 = FetchStrings(files,'Observables/'+self.section+'/Species2')[0]
        grs = FetchDataArray(files,'Observables/'+self.section+'/y',self.startCut)
        Box = FetchData(files,'System/Box',0,0)
        return rs, grs

class NofR(Observable):
    def GetData(self, files):
        rs = FetchData(files,'Observables/'+Section+'/x',0,0)
        rs = rs[0]
        nrs = FetchDataArray(files,'Observables/'+Section+'/y',StartCut)
        Box = FetchData(files,'System/Box',0,0)
        return rs, nrs

class SectorCount(Observable):
    def GetData(self, files):
        data = FetchDataArray2(files,'Observables/CycleCount/SectorCount',self.startCut)
        Data = {}
        nFiles = len(data)
        nSamples = 0
        for i in range(0,nFiles):
            for j in range(self.startCut,len(data[i])):
                label = int(data[i][j][0])
                val = data[i][j][1]
                nSamples += val
                if label in Data:
                    Data[label].append(val)
                else:
                    Data[label] = [val]
        return sorted(Data.iterkeys()), Data

    def DoStatistics(self, xs, ys):
        stats = []
        tot = 0
        for label in xs:
            ys[label] = sum(ys[label])
            tot += ys[label]
        for label in xs:
            ys[label] = float(ys[label])/tot
            stats.append([ys[label],0.,tot])
        return stats

class SectorEnergy(Observable):
    def GetData(self, files):
        Data = {}
        for f in files:
            try:
                AllData = f['Observables/Energy/SectorEnergy'][self.startCut:]
            except KeyError:
                AllData = f['Observables/Energy/PermEnergy'][self.startCut:]
            for data in AllData:
                [label, val, N, err] = list(data)
                label = int(label)
                if label in Data:
                    Data[label].append([val,err,N])
                else:
                    Data[label] = [[val,err,N]]
        return sorted(Data.iterkeys()), Data

    def DoStatistics(self, xs, ys):
        stats = []
        for label in xs:
            [mean,err,N] = Stats.UnweightedReAvg(ys[label])
            stats.append([mean,err,N])
        return stats

class CycleCount(Observable):
    def GetData(self, files):
        data = FetchDataArray(files,'Observables/CycleCount/y',self.startCut)
        return range(len(data)), data

def AdjustValuesByImportanceWeight(basename,Sections,clone=''):
    print '* Adjusting Values By Importance Weight'
    try:
        f = open(basename+'/ImportanceWeight'+clone+'.dat','r')
    except:
        f = open(basename+'/Sign'+clone+'.dat','r')
    for line in f:
        line = line.split()
        sgn = float(line[1])
        if len(line) == 5:
            dsgn = float(line[4])
        else:
            dsgn = float(line[2])
    f.close()
    for Section in Sections:
        try:
            g = open(basename+'/'+Section+clone+'.dat','r')
            h = open(basename+'/Adj'+Section+clone+'.dat','w')
            for line in g:
                line = line.split()
                label = line[0]
                x = float(line[1])
                if len(line) == 5:
                    dx = float(line[4])
                else:
                    dx = float(line[2])
                z = x
                dz = dx
                if (sgn != 0) and (label != 'VTailLongRange') and (label != 'VTailShortRange'):
                    z = x/sgn
                    if (x != 0) and (dx != float('inf')) and (dsgn != float('inf')):
                        dz = abs(z)*sqrt((dx/x)**2 + (dsgn/sgn)**2)
                #print label, x, dx
                h.write('%s %e %e\n' % (label, z, dz))
            g.close()
            h.close()
        except:
            print 'Error in', Section

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
                y = float(line[1])
                if len(line) == 5: # label y sig(y) kappa(y) err(y)
                    dy = float(line[4])
                else:
                    dy = float(line[2])
                if len(line) == 4:
                    label = getnumber(line[0])
                    N = float(line[3])
                    if label in data:
                        data[label].append([y,dy,N,i])
                    else:
                        data[label] = [[y,dy,N,i]]
                else:
                    label = line[0]
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
                if Section == 'SectorCount':
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
    parser.add_option("-t", "--startclone", action="store", type="int", dest="StartClone", default=0, help="Starting clone [default: %default]")
    parser.add_option("-y", "--endclone", action="store", type="int", dest="EndClone", default=0, help="Ending clone [default: %default]")
    parser.add_option("-r", "--excluderun", action="callback", callback=cb, dest="ExcludedRuns", default=[], help="Exclude specific runs (e.g. the first one, 0) [default: %default]")
    parser.add_option("-q", "--startingrun", action="store", type="int", dest="StartingRun", default=0, help="Start at specific run (e.g. the first one, 0) [default: %default]")
    parser.add_option("-b", "--blocksize", action="store", type="int", dest="BlockSize", default=1, help="Block size [default: %default]")
    (opts, args) = parser.parse_args()
    if mpi.rank == 0:
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
            basename = string.join(splitBaseName[0:-2],'.')
            SingleFile = True

    i = 0
    while i < opts.StartingRun:
        opts.ExcludedRuns.append(str(i))
        i += 1
    if mpi.rank == 0:
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

    # Use clone range
    if opts.EndClone != opts.StartClone:
        clone = opts.StartClone
        run = 0
        while str(run) in opts.ExcludedRuns:
            run += 1
        while clone != opts.EndClone + 1:
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
    if mpi.rank == 0:
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
        files = [[basedir+basename+'.'+str(run)+'.'+clones[0]+'.h5']]*cloneCount

    # Create directory
    try:
        print 'Creating New Directory', basename
        os.makedirs(basename)
    except OSError:
        pass

    # Run Through Clones
    cloneStrs = ['_'+clone for clone in clones]
    cloneCountPerProc = len(clones) / mpi.size
    for j in range(0,cloneCountPerProc):
        i = mpi.rank * cloneCountPerProc + j
        print '*** Working on clone', clones[i], '***'

        # Open h5 files
        h5files = []
        for f in files[i]:
            try:
                h5files.append(h5.File(f,'r'))
            except:
                sys.stderr.write('Error with: '+f+'\n')
                sys.exit()

        # File List Info
        Sections = []
        nrun = len(h5files)
        print 'Found', nrun, 'file(s).'
        if nrun == 0:
            sys.stderr.write('ERROR: No Files!\n')
            sys.exit()
        else:
            Sections = GetSubSections(h5files,'Observables')

        # Special Observables
        if 'SectorEnergy' in GetSubSections(h5files,'Observables/Energy'):
            Sections.append('SectorEnergy')
        if 'PermEnergy' in GetSubSections(h5files,'Observables/Energy'):
            Sections.append('SectorEnergy')
        if 'SectorCount' in GetSubSections(h5files,'Observables/CycleCount'):
            Sections.append('SectorCount')

        # Process Observables
        AllSections, AdjSections = [], []
        for Section in Sections:
            O, adj = GetObservable(basename,StartCuts[i],Section,cloneStrs[i])
            if opts.DoAnalysis:
                if O.Process(h5files):
                    AllSections.append(Section)
            else:
                AllSections.append(Section)
            if adj:
                AdjSections.append(Section)

        # Process Moves
        Sections = GetSubSections(h5files,'Moves')
        for Section in Sections:
            M, adj = GetMove(basename,StartCuts[i],Section,cloneStrs[i])
            if opts.DoAnalysis:
                if M.Process(h5files):
                    AllSections.append(Section)
            else:
                AllSections.append(Section)

        # Adjust by sign
        if opts.TrackSign:
            AdjustValuesByImportanceWeight(basename,AdjSections,cloneStrs[i])

        # Close files
        for file in h5files:
            file.close()
    mpi.barrier()

    # Combine Clones
    if mpi.rank == 0:
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
            AdjustValuesByImportanceWeight(basename,AdjSections)

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
