print 'pylab'
import pylab
print 'IO'
from IO import *
print 'sys'
import sys
print 'numpy'
from numpy import *
print 'done'



def dToe(myString):
  newString=""
  for i in myString:
    if i!='d' and i!="D":
        newString=newString+i
    else:
      newString=newString+"e"
  return newString
      
class LongRangeParser:
  def __init__(self,baseName,outfilename,kcut):
     self.fileName=baseName
     self.a=IOSectionClass()
     self.a.OpenFile(outfilename)
     self.kcut=kcut
  def ReadYK(self):
     vals=pylab.mlab.load(self.fileName+"yk")
     self.a.NewSection("LongRange")
     self.a.WriteVar("Type","yk")
     self.a.WriteVar("kcut",self.kcut)
     kPoints=vals[:,0].copy()
     uk=vals[:,1].copy()
     self.a.WriteVar("kPoints",kPoints)
     self.a.WriteVar("u_k",uk)
     potgen_inFile=open(self.fileName+"in")
     a=potgen_inFile.readlines()
     mass1=float(dToe(a[1].split()[2]))
     mass2=float(dToe(a[2].split()[2]))
     self.a.WriteVar("mass1",mass1)
     self.a.WriteVar("mass2",mass2)
     ndim=int(a[5].split()[5])
     box=zeros(ndim,float)
     self.a.WriteVar("ndim",ndim)

     for dim in range(ndim):
       box[dim]=float(dToe(a[5].split()[8+dim]))
     self.a.WriteVar("Box",box)
     print box
     print mass1,mass2
  def Done(self):
    self.a.FlushFile()
    self.a.CloseFile()


class Squarer2HDFParser:
  def __init__(self,filename):
    self.f = open(filename,'r')
    self.basename = filename[:-2]
    outfilename = self.basename + 'h5'
    self.a=IOSectionClass()
    self.a.NewFile(outfilename)
    self.numFits = 0
    self.SpitPotential = True
    self.Ucounter = 0 
    self.dUcounter = 0
    self.spec1 = ''
    self.spec2 = ''

  def Done(self):
    print "Done."
    self.f.close()
    self.a.FlushFile()
    self.a.CloseFile()

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
    print "ProcessSquarerInfo"
    ### collect squarer run info
    print self.a
    print self.a.NewSection
    self.a.NewSection("Squarer")
    self.a.NewSection("Units")
    g = self.next()
    check(g,'UNITS')
    self.a.WriteVar("Temp",self.next())
    self.a.WriteVar("Length",self.next())
    self.a.CloseSection()
    g = self.next() 
    check(g,'TYPE')
    self.a.NewSection("Type")
    self.spec1 = self.next()
    self.a.WriteVar("Species",self.spec1)
    self.a.WriteVar("Lambda",float(dToe(self.next())))
    self.a.CloseSection()
    g = self.next() 
    check(g,'TYPE')
    self.a.NewSection("Type")
    self.spec2 = self.next()
    self.a.WriteVar("Species",self.spec2)
    self.a.WriteVar("Lambda",float(dToe(self.next())))
    self.a.CloseSection()
    
    ### get important stats for remainder of read
    self.find("SQUARER")
    self.next()
    self.next()
    self.next()
    self.numFits = int(self.next())
    self.a.WriteVar("NumFits",self.numFits)
    self.find("POT")
    self.next()
    self.next()
    self.kcut=float(self.next())
    print "KCUT IS",self.kcut
    print "VIMAGE"
    self.find("VIMAGE")
    print "VIMAGE"
    self.vimage=float(self.next())
    self.a.WriteVar("vimage",self.vimage)
    self.a.CloseSection()
    #self.a.WriteVar("NumFits",self.numFits)
    self.a.FlushFile()
  ### end SquarerInfo ###

  ### get potential
  def ProcessPotential(self):
    print "ProcessPotential"
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
    self.a.NewSection("Potential")
    self.a.WriteVar("Data",u)
    self.a.CloseSection()
    self.a.FlushFile()
    # dump potential to ASCII
    if(self.SpitPotential):
      pout = open(self.basename + 'potential.dat', 'w')
      for i in range(0,len(u)):
        pout.write(str(u[i]) + '\n')
      pout.close()
  ### end ProcessPotential ###

  def ProcessU(self):
    print "ProcessU"
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
   
    print NMax,derv
    SectionTitle = 'Ukj' + str(self.Ucounter)
    self.Ucounter += 1
    #SectionTitle = ''
    #print "Sectiontitle reset"
    #if(order == 1):
    #  self.Ucounter += 1
    #  print "setting sectiontitle",self.SectionTitle,"to"
    #  self.SectionTitle = 'Ukj' + str(self.Ucounter)
    #  print order,self.SectionTitle
    #elif(order == 2):
    #  self.dUcounter += 1
    #  print "setting sectiontitle",self.SectionTitle,"to"
    #  self.SectionTitle = 'dUkj_dBeta' + str(self.dUcounter)
    #  print order,self.SectionTitle
    #else:
    #  print "What's going on? order is",order
    #  sys.exit()
    #print order,self.SectionTitle
    print "Going to write section",SectionTitle
    self.a.NewSection(SectionTitle)
    self.a.NewSection("Grid")
    self.a.WriteVar("NumGridPoints",numPts)
    self.a.WriteVar("Type",gridType)
    self.a.WriteVar("Start",start)
    self.a.WriteVar("End",end)
    self.a.CloseSection()
    self.a.WriteVar("NumUkj",numUkj)
    self.a.WriteVar("NumTau",numTau)
    self.a.WriteVar("NMax",NMax)
    self.a.WriteVar("Derv",derv)
    self.a.WriteVar("Rank",rnk)
    self.a.WriteVar("Taus",Taus)
    self.a.WriteVar("Data",Ukj)
    self.a.CloseSection()
    self.a.FlushFile()
  ### end ProcessU ###

  def ProcessdU_dBeta(self):
    print "ProcessdU_dBeta"
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
   
    print NMax,derv
    SectionTitle = 'dUkjdBeta' + str(self.dUcounter)
    self.dUcounter += 1
    #SectionTitle = ''
    #print "Sectiontitle reset"
    #if(order == 1):
    #  self.Ucounter += 1
    #  print "setting sectiontitle",self.SectionTitle,"to"
    #  self.SectionTitle = 'Ukj' + str(self.Ucounter)
    #  print order,self.SectionTitle
    #elif(order == 2):
    #  self.dUcounter += 1
    #  print "setting sectiontitle",self.SectionTitle,"to"
    #  self.SectionTitle = 'dUkj_dBeta' + str(self.dUcounter)
    #  print order,self.SectionTitle
    #else:
    #  print "What's going on? order is",order
    #  sys.exit()
    #print order,self.SectionTitle
    print "Going to write section",SectionTitle
    self.a.NewSection(SectionTitle)
    self.a.NewSection("Grid")
    self.a.WriteVar("NumGridPoints",numPts)
    self.a.WriteVar("Type",gridType)
    self.a.WriteVar("Start",start)
    self.a.WriteVar("End",end)
    self.a.CloseSection()
    self.a.WriteVar("NumUkj",numUkj)
    self.a.WriteVar("NumTau",numTau)
    self.a.WriteVar("NMax",NMax)
    self.a.WriteVar("Derv",derv)
    self.a.WriteVar("Rank",rnk)
    self.a.WriteVar("Taus",Taus)
    self.a.WriteVar("Data",Ukj)
    self.a.CloseSection()
    self.a.FlushFile()
  ### end ProcessdU_dBeta ###
  
  def ProcessSampling(self):
    print "ProcessSampling",
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
    print derv
  
    # load array from .dm file
    Ukj = zeros([numPts, numUkj, numTau]) + 0.
    for cT in range(0,numTau):
      for cU in range(0,numUkj):
        for cG in range(0,numPts):
          #print derv,cT,"/",numTau,cU,"/",numUkj,cG,"/",numPts,
          Ukj[cG, cU, cT] = float(self.next())
  
    Taus = zeros(numTau) + 0.
    tau0 = loTau/2
    for t in range(0,numTau):
      tau0 *= 2.
      Taus[t] = tau0
    
    self.a.NewSection("Sampling")
    self.a.WriteVar("Taus",Taus)
    self.a.WriteVar("NumUkj",numUkj)
    self.a.WriteVar("NumTau",numTau)
    self.a.WriteVar("Derv",derv)
    self.a.WriteVar("Data",Ukj)
    self.a.CloseSection()
    self.a.FlushFile()
  ### end ProcessSampling ###

### end class Squarer2HDFParser ###

def check(s, sCheck):
  if(s != sCheck):
    print "MISMATCH: expected",sCheck,"got",s
    return(False)
  return(True);

### main ###
def main():
  if(len(sys.argv)!=2):
    print "Usage: python dmparse.py squarer_output_file.dm"
    sys.exit()
  print 'Squarer Parser'
  sq = Squarer2HDFParser(sys.argv[1])
  print 'Process Squarer Info'
  sq.ProcessSquarerInfo()
  print 'Process Potential'
  sq.ProcessPotential()
  for sec in range(0,sq.numFits + 1):
    sq.ProcessU()
  for sec in range(0,sq.numFits + 1):
    sq.ProcessdU_dBeta()
  print "Species are",sq.spec1,sq.spec2
  if(sq.spec1 == sq.spec2):
    sq.ProcessSampling()
    sq.ProcessSampling()
  else:
    print "Skipping sampling table for different species"
  sq.Done()
  
  basename = sys.argv[1][:-2]
  outfilename =basename + 'h5'
  infilename=basename+"yk"
  if os.path.exists(infilename):
      print "About to process the long range file",infilename
      #should check to see if it's consistent somehow?
      longRangeParse=LongRangeParser(basename,outfilename,sq.kcut)
      longRangeParse.ReadYK()
      longRangeParse.Done()

main()
