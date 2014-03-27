import sys
import os
from numpy import *
import h5py as h5

def cds(group,name,data):
  a = zeros(1,float)
  try:
    if not data.dtype == a.dtype:
      a = array(data)
    else:
      a = data
  except:
    a = array([data])
  ds = group.create_dataset(name,a.shape,a.dtype)
  ds[:] = data

def check(s, sCheck):
  if(s != sCheck):
    print "MISMATCH: expected",sCheck,"got",s
    return(False)
  return(True);


class LongRangeParser:
  def __init__(self,baseName,outfilename,kcut):
     self.fileName = baseName
     self.a = h5.File(outfilename)
     self.kcut = kcut
  def ReadYK(self):
     vals = loadtxt(self.fileName+"yk")
     self.b = self.a.create_group("LongRange")
     cds(self.b,'Type','yk')
     cds(self.b,'kcut',self.kcut)
     cds(self.b,'kPoints',vals[:,0])
     cds(self.b,'u_k',vals[:,1])
     potgen_inFile=open(self.fileName+"dm")
     a=potgen_inFile.readlines()
     mass1=float(a[1].split()[2])
     mass2=float(a[2].split()[2])
     cds(self.b,'mass1',mass1)
     cds(self.b,'mass2',mass2)
     ndim=int(a[5].split()[5])
     box=zeros(ndim,float)
     cds(self.b,'ndim',ndim)
     box[0]=float(a[5].split()[-1])
     for dim in range(ndim):
       #box[dim]=float(a[5].split()[8+dim])
       box[dim]=box[0] # HACK ONLY CUBIC BOXES
     cds(self.b,'Box',box)
     print box
     print mass1,mass2
  def Done(self):
    self.a.close()


class Squarer2HDFParser:
  def __init__(self,filename,UseVimage):
    self.f = open(filename,'r')
    self.basename = filename[:-2]
    outfilename = self.basename + 'h5'
    self.a = h5.File(outfilename,'w')
    self.numFits = 0
    self.SpitPotential = True
    self.Ucounter = 0
    self.dUcounter = 0
    self.samplingcounter = 0
    self.spec1 = ''
    self.spec2 = ''
    self.UseVimage = UseVimage

  def Done(self):
    print "Done."
    self.f.close()
    self.a.close()

  def next(self):
    w = ''
    c = self.f.read(1)
    # check for EOF
    if(c == ''):
      print "ENCOUNTERED EOF"
      sys.exit()
      return w
    empty = True
    while(empty):
      while(c!=' ' and c!='\t' and c!='\n' and c!=''):
        w += c
        c = self.f.read(1)
        empty = False
      if(empty):
        c = self.f.read(1)
    #print "parsed",w
    return w

  def find(self, target):
    s = self.next()
    while(s != target):
      s = self.next()

  def ProcessSquarerInfo(self):
    ### collect squarer run info
    self.b = self.a.create_group("Squarer")
    g = self.next()
    check(g,'UNITS')
    self.c = self.b.create_group("Units")
    cds(self.c,'Temp',self.next())
    cds(self.c,'Length',self.next())
    g = self.next()
    check(g,'TYPE')
    self.c = self.b.create_group("Type_0")
    self.spec1 = self.next()
    cds(self.c,'Species',self.spec1)
    cds(self.c,'Lambda',float(self.next()))
    g = self.next()
    check(g,'TYPE')
    self.c = self.b.create_group("Type_1")
    self.spec2 = self.next()
    cds(self.c,'Species',self.spec2)
    cds(self.c,'Lambda',float(self.next()))

    #### get important stats for remainder of read
    self.find("SQUARER")
    self.next()
    self.next()
    self.next()
    self.numFits = int(self.next())
    cds(self.b,'NumFits',self.numFits)
    self.find("POT")
    if self.UseVimage:
      self.next()
      self.next()
      self.kcut=float(self.next())
      print "KCUT IS",self.kcut
      cds(self.b,'kcut',self.kcut)
      self.find("VIMAGE")
      self.vimage=float(self.next())
      print "VIMAGE IS",self.vimage
      cds(self.b,'vimage',self.vimage)
  ### end SquarerInfo ###

  ### get potential
  def ProcessPotential(self):
    self.find("RANK")
    print "RANK"
    self.next()
    size = int(self.next())
    self.find("BEGIN")
    self.next()
    self.next()
    u = zeros(size) + 0.
    for i in range(0,size):
      u[i] = float(self.next())
    self.b = self.a.create_group("Potential")
    cds(self.b,'Data',u)
    # dump potential to ASCII
    if(self.SpitPotential):
      pout = open(self.basename + 'potential.dat', 'w')
      for i in range(0,len(u)):
        pout.write(str(u[i]) + '\n')
      pout.close()
  ### end ProcessPotential ###

  def ProcessU(self):
    self.find("RANK")
    rnk = int(self.next())
    check(rnk,3)
    numPts = int(self.next())
    numUkj = int(self.next())
    numTau = int(self.next())

    self.find("GRID")
    self.next()
    gridType = self.next()
    start = float(self.next())
    end = float(self.next())
    self.find("GRID")
    self.next()
    checkgrid = self.next()
    check(checkgrid, "LOG")
    loTau = float(self.next())
    hiTau = float(self.next())

    self.find("BEGIN")
    self.next()
    self.next()
    self.next()
    self.next()
    NMax = int(self.next())
    self.next()
    self.next()
    derv = int(self.next())
    #check(derv,order)

    # load array from .dm file
    Ukj = zeros([numPts, numUkj, numTau]) + 0.
    for cT in range(0,numTau):
      for cU in range(0,numUkj):
        for cG in range(0,numPts):
          Ukj[cG, cU, cT] = float(self.next())

    Taus = zeros(numTau) + 0.
    tau0 = loTau/2
    for t in range(0,numTau):
      tau0 *= 2.
      Taus[t] = tau0

    print 'NMax:', NMax, 'derv:', derv
    SectionTitle = 'Ukj' + str(self.Ucounter)
    self.Ucounter += 1
    self.b = self.a.create_group(SectionTitle)
    self.c = self.b.create_group("Grid")
    cds(self.c,"NumGridPoints",numPts)
    cds(self.c,"Type",gridType)
    cds(self.c,"Start",start)
    cds(self.c,"End",end)
    cds(self.b,"NumUkj",numUkj)
    cds(self.b,"NumTau",numTau)
    cds(self.b,"NMax",NMax)
    cds(self.b,"Derv",derv)
    cds(self.b,"Rank",rnk)
    cds(self.b,"Taus",Taus)
    cds(self.b,"Data",Ukj)
  ### end ProcessU ###

  def ProcessdU_dBeta(self):
    self.find("RANK")
    rnk = int(self.next())
    check(rnk,3)
    numPts = int(self.next())
    numUkj = int(self.next())
    numTau = int(self.next())

    self.find("GRID")
    self.next()
    gridType = self.next()
    start = float(self.next())
    end = float(self.next())
    self.find("GRID")
    self.next()
    checkgrid = self.next()
    check(checkgrid, "LOG")
    loTau = float(self.next())
    hiTau = float(self.next())

    self.find("BEGIN")
    self.next()
    self.next()
    self.next()
    self.next()
    NMax = int(self.next())
    self.next()
    self.next()
    derv = int(self.next())
    #check(derv,order)

    # load array from .dm file
    Ukj = zeros([numPts, numUkj, numTau]) + 0.
    for cT in range(0,numTau):
      for cU in range(0,numUkj):
        for cG in range(0,numPts):
          Ukj[cG, cU, cT] = float(self.next())

    Taus = zeros(numTau) + 0.
    tau0 = loTau/2
    for t in range(0,numTau):
      tau0 *= 2.
      Taus[t] = tau0

    print 'NMax:', NMax, 'derv:', derv
    SectionTitle = 'dUkjdBeta' + str(self.dUcounter)
    self.dUcounter += 1
    self.b = self.a.create_group(SectionTitle)
    self.c = self.b.create_group("Grid")
    cds(self.c,"NumGridPoints",numPts)
    cds(self.c,"Type",gridType)
    cds(self.c,"Start",start)
    cds(self.c,"End",end)
    cds(self.b,"NumUkj",numUkj)
    cds(self.b,"NumTau",numTau)
    cds(self.b,"NMax",NMax)
    cds(self.b,"Derv",derv)
    cds(self.b,"Rank",rnk)
    cds(self.b,"Taus",Taus)
    cds(self.b,"Data",Ukj)
  ### end ProcessdU_dBeta ###

  def ProcessSampling(self):
    self.find("RANK")
    rnk = int(self.next())
    check(rnk,3)
    numPts = int(self.next())
    numUkj = int(self.next())
    numTau = int(self.next())

    self.find("GRID")
    self.next()
    checkgrid = self.next()
    check(checkgrid, "LOG")
    loTau = float(self.next())
    hiTau = float(self.next())

    self.find("BEGIN")
    self.next()
    self.next()
    derv = int(self.next())
    print 'derv:', derv

    # load array from .dm file
    Ukj = zeros([numPts, numUkj, numTau]) + 0.
    for cT in range(0,numTau):
      for cU in range(0,numUkj):
        for cG in range(0,numPts):
          Ukj[cG, cU, cT] = float(self.next())

    Taus = zeros(numTau) + 0.
    tau0 = loTau/2
    for t in range(0,numTau):
      tau0 *= 2.
      Taus[t] = tau0

    SectionTitle = 'Sampling_' + str(self.samplingcounter)
    self.samplingcounter += 1
    self.b = self.a.create_group(SectionTitle)
    cds(self.b,"NumUkj",numUkj)
    cds(self.b,"NumTau",numTau)
    cds(self.b,"Derv",derv)
    cds(self.b,"Taus",Taus)
    cds(self.b,"Data",Ukj)
  ### end ProcessSampling ###

### end class Squarer2HDFParser ###

def Parse(dmfile):
  basename = dmfile[:-2]
  infilename=basename+"yk"
  UseVimage = False
  if os.path.exists(infilename):
    UseVimage = True

  print 'Squarer Parsing'
  sq = Squarer2HDFParser(dmfile,UseVimage)
  print "Process Squarer Info"
  sq.ProcessSquarerInfo()
  print "Process Potential"
  sq.ProcessPotential()
  print "Process U"
  for sec in range(0,sq.numFits + 1):
    sq.ProcessU()
  print "Process dU/dBeta"
  for sec in range(0,sq.numFits + 1):
    sq.ProcessdU_dBeta()
  print "Species are",sq.spec1,sq.spec2
  print "Process Sampling",
  if(sq.spec1 == sq.spec2):
    sq.ProcessSampling()
    sq.ProcessSampling()
  else:
    print "Skipping sampling table for different species"
  sq.Done()

  outfilename =basename + 'h5'
  infilename=basename+"yk"
  if os.path.exists(infilename):
    print "About to process the long range file",infilename
    #should check to see if it's consistent somehow?
    longRangeParse=LongRangeParser(basename,outfilename,1.0)#sq.kcut)
    longRangeParse.ReadYK()
    longRangeParse.Done()

def usage():
  print "Usage:  %s dmfile.dm" % os.path.basename(sys.argv[0])

def main(argv=None):
  if argv is None:
    argv = sys.argv
  if "-h" in argv or "--help" in argv:
    usage()
    sys.exit(2)

  Parse(argv[1])

if __name__ == "__main__":
  sys.exit(main())
