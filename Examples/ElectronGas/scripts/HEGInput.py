import sys
import os
import subprocess
from math import log,sqrt,pi
import HF
import dmparse

def GenSquarerInput(fileNameSquarer,N,D,T,L,lam,eps,M,nMax,nTemp,nSquare,UseMinTau=True):

  if nTemp > 8:
    print 'WARNING: Max nTemp is 8!'
    nTemp = 8
  MinTau = 1./(T*M)
  MaxTau = MinTau*(2**(nTemp-1))

  #rMax = 0.75*sqrt(3.)*L
  rMax = L

  print 'Creating '+fileNameSquarer+'.in'
  f = open(fileNameSquarer+'.in','w')
  f.write('UNITS r a\n')
  f.write('TYPE e %f\n' % (lam))
  f.write('TYPE e %f\n' % (lam))
  f.write('GRID 450 LOG 0.0001 %f\n' % (rMax))
  f.write('SQUARER %f %i %i 3 30 %i\n' % (1./MaxTau,nTemp,D,nSquare))
  Box = ' '.join([str(L) for i in range(0,D)])
  UsedTau = MinTau if UseMinTau else MaxTau
  f.write('POT COULOPT %f %f %f %i 0.D0 1.0 %s\n' %((L/2.0) - sqrt(lam*UsedTau),nMax*2.0*pi/L,eps,D,Box))
  f.close()

def BuildPairAction(prefix,lam,D):
  f = open(prefix+'.PairAction','w')

  f.write('Section (Fits)')
  f.write('\n{')
  f.write('\n  string Type="DavidFit";')
  f.write('\n  int NumOffDiagonalTerms = 3;')
  f.write('\n  Section (Particle1)')
  f.write('\n  {')
  f.write('\n    string Name = "e";')
  f.write('\n    double lambda ='+str(lam)+';')
  f.write('\n    int Ndim = '+str(D)+';')
  f.write('\n  }')
  f.write('\n  Section (Particle2)')
  f.write('\n  {')
  f.write('\n    string Name = "e";')
  f.write('\n    double lambda ='+str(lam)+';')
  f.write('\n    int Ndim = '+str(D)+';')
  f.write('\n  }')
  f.write('\n  string Daviddmfile = "inputs/' + prefix +'.h5";')
  f.write('\n}')
  f.close()

  f = open(prefix+'.PairAction','r')
  for line in f:
    print line
  f.close()

def GenPIMCInput(fileNamePIMC,fileNameSquarer,N,D,T,L,lam,eps,M,nMax,pol,ParticleType,TrackSign,restart,PPC):

  # Get Statistics of particles
  Statistics = ''
  if ParticleType == 0:
    Statistics = 'BOSON'
    NodeType = 'NONE'
  elif ParticleType == 1:
    Statistics = 'FERMION'
    NodeType = 'FREE'
  elif ParticleType == 2:
    Statistics = 'BOLZTMANNON'
    NodeType = 'NONE'
  else:
    print 'Unknown Particle Type:', ParticleType
    sys.exit()

  # Whether or not to count permutation sectors
  CountPerms = 0
  if TrackSign:
    if Statistics == 'FERMION':
      CountPerms = 0
    else:
      CountPerms = 1
  else:
    CountPerms = 0

  # How to initiate simulation
  if restart:
    InitPaths = 'RESTART'
  else:
    InitPaths = 'BCC'

  # Particle lists
  if pol:
    Particles = [['eUp','e',N,lam,Statistics,NodeType,InitPaths]]
  else:
    assert(not N%2)
    Particles = [['eUp','e',N/2,lam,Statistics,NodeType,InitPaths],['eDown','e',N/2,lam,Statistics,NodeType,InitPaths]]

  # This is squarer's kmax
  kMax = nMax*2.0*pi/L

  # For the structure factor, we actually want it much larger
  #kMax *= 2

  # Build PIMC++ Input Files
  beta = 1/T
  tau = beta/M
  if M > 2*PPC:
    MaxLevels = int(log(M/(2*PPC))/log(2))
  else:
    MaxLevels = 0
    print "Warning: M < 2*PPC !"
    print "Setting MaxLevels = 0"
  print tau, MaxLevels, M, PPC

  # Generate Input String
  inputString = """
Section (Output)
{
  string OutFileBase = \"rawdata/"""+fileNamePIMC+"""\";
  bool Restart = true;
}

Section (Parallel)
{
  int ProcsPerClone = """+str(PPC)+""";
}

Section (System)
{
  double tau = """+str(tau)+""";
  int NumTimeSlices = """+str(M)+""";
  Array<double,1> Box("""+str(D)+""") = ["""+','.join([str(L) for i in range(0,D)])+"""];
  Array<bool,1> IsPeriodic("""+str(D)+""") = ["""+','.join(['true' for i in range(0,D)])+"""];
  Section (Particles)
  {"""
  for [Name,Type,N,Lambda,Statistics,NodeType,InitPaths] in Particles:
    inputString += """
    Section (Species)
    {
      string Name = \""""+Name+"""\";
      string Type = \""""+Type+"""\";
      double lambda = """+str(Lambda)+""";
      string Statistics = \""""+Statistics+"""\";"""
    if Statistics == 'FERMION':
      inputString += """
      string NodeType = \""""+NodeType+"""\";"""
    inputString += """
      int NumParticles = """+str(N)+""";
      int NumDim = """+str(D)+""";
      string InitPaths = \""""+InitPaths+"""\";
      string File = "rawdata/"""+fileNamePIMC+"""";
    }"""
  inputString += """
  }
}

Section (Actions)
{
  int NumImages = 1;
  int MaxLevels = """+str(MaxLevels)+""";
  Array<string,1> PairActionFiles(1) = ["inputs/"""+fileNameSquarer+""".PairAction"];

  Section (Action)
  {
    string Name = "ShortRange";
    string Type = "ShortRange";
  }

  Section (Action)
  {
    string Name = "DavidLongRange";
    string Type = "DavidLongRange";
    double kCutoff = """+str(kMax)+""";
  }"""
  for [Name,Type,N,Lambda,Statistics,NodeType,InitPaths] in Particles:
    if Statistics == 'FERMION':
      inputString += """
  Section (NodalAction)
  {
    string Type = \""""+NodeType+ """\";
    string Species = \""""+Name+"""\";
    bool StoreNodeDist = true;
    bool StoreNodeDet = true;
    bool UseNoDist = false;
    bool UseHybridDist=true;
  }"""
  inputString += """
}

Section (Observables)
{
  Section (Observable)
  {
    string Type = "Energy";
    string Name = "Energy";
    string Description = "Total Energy";
    int Frequency = 1;
    double HistStart = 0.0;
    double HistEnd = 1.0;
    int HistPoints = 2;
    bool TrackSign = """
  if TrackSign:
    inputString += """true;"""
  else:
    inputString += """false;"""
    inputString += """
    bool CountPerms = """
  if CountPerms:
    inputString += """true;"""
  else:
    inputString += """false;"""
  inputString += """
  }"""
  if CountPerms:
    inputString += """
  Section (Observable)
  {
    string Type = "CycleCount";
    string Name = "CycleCount";
    int Frequency = 1;
  }"""
  if TrackSign:
    inputString += """
  Section (Observable)
  {
    string Type = "Sign";
    string Name = "Sign";
    string Description = "Total Energy";
    int Frequency = 1;
  }"""
  count = 0
  for [Name1,Type1,N2,Lambda1,Statistics1,NodeType1,InitPaths1] in Particles:
    for [Name2,Type2,N2,Lambda2,Statistics2,NodeType2,InitPaths2] in Particles[count:]:
      inputString += """
  Section (Observable)
  {
    string Type = "PairCorrelation";
    string Name = "PairCorrelation"""+Name1+Name2+"""\";
    string Description = "Pair Correlation";
    string Species1 = \""""+Name1+"""\";
    string Species2 = \""""+Name2+"""\";
    Section (Grid)
    {
      string Type = "Linear";
      double start = 0.0;
      int NumPoints = 100;
    }
    int Frequency = 1;
  }

  Section (Observable)
  {
    string Type = "StructureFactor";
    string Name = "StructureFactor"""+Name1+Name2+"""\";
    string Description = "Structure Factor";
    double kCutoff = """+str(kMax)+""";
    string Species1 = \""""+Name1+"""\";
    string Species2 = \""""+Name2+"""\";
    int Frequency = 1;
  }"""
      if (Statistics1 == 'FERMION' and Statistics2 == 'FERMION'):
        inputString += """
  Section (Observable)
  {
    string Type = "RefPairCorrelation";
    string Name = "RefPairCorrelation"""+Name1+Name2+"""\";
    string Description = "Ref Pair Correlation";
    string Species1 = \""""+Name1+"""\";
    string Species2 = \""""+Name2+"""\";
    Section (Grid)
    {
      string Type = "Linear";
      double start = 0.0;
      int NumPoints = 100;
    }
    int Frequency = 1;
  }

  Section (Observable)
  {
    string Type = "RefStructureFactor";
    string Name = "RefStructureFactor"""+Name1+Name2+"""\";
    string Description = "Structure Factor";
    double kCutoff = """+str(kMax)+""";
    string Species1 = \""""+Name1+"""\";
    string Species2 = \""""+Name2+"""\";
    int Frequency = 1;
  }"""
    count += 1
  inputString += """
  Section (Observable)
  {
    bool AllClones = true;
    string Type = "PathDump";
    string Name = "PathDump";
    int TemporalFrequency = """+str(1700)+""";
  }

  Section (Observable)
  {
    string Type = "TimeAnalysis";
    string Name = "TimeAnalysis";
    int Frequency = 1;
  }

}

Section (Moves)
{"""
  for [Name,Type,N,Lambda,Statistics,NodeType,InitPaths] in Particles:
    inputString += """

  Section (Move)
  {
    string Type = "BisectionBlock";
    string Name = "BisectionBlock"""+Name+"""\";
    //Array<string,1> HigherLevelActions(2) = ["ShortRange","DavidLongRange"];
    Array<string,1> SamplingActions(2) = ["ShortRange","DavidLongRange"];"""
    if (Statistics == 'BOSON' or Statistics == 'FERMION'):
      inputString += """
    string PermuteType = "TABLE";"""
      if (Statistics == 'FERMION'):
        inputString += """
    Array<double,1> Gamma(4) = [1.0,0.0,1.0,0.0];"""
      else:
        inputString += """
    Array<double,1> Gamma(4) = [1.0,1.0,1.0,1.0];"""
    else:
      inputString += """
    string PermuteType = "NONE";"""
    inputString += """
    double epsilon = 1e-10;
    string Species = \""""+Name+"""\";
    int NumLevels = """+str(MaxLevels)+""";
    int StepsPerBlock = """+str(N)+""";
  }"""

    if (NodeType != 'NONE'):
      inputString += """
  Section (Move)
  {
    string Type = "RefSlice";
    string Name = "RefSlice"""+Name+"""\";
    string PermuteType = "NONE";
    int NumLevels = """+str(MaxLevels)+""";
    string Species = \""""+Name+"""\";
    Array<string,1> SamplingActions(2) = ["ShortRange","DavidLongRange"];
  }"""
  inputString += """
  Section (Move)
  {
    string Type = "ShiftMove";
    string Name = "Shift";
  }
}

Section (Algorithm)
{

  Section (Loop)
  {
    int Steps = 10000000;

    Section (Loop)
    {
      int Steps = 10;
      bool Equilibrate = false;"""
  count = 0
  for [Name1,Type1,N1,Lambda1,Statistics1,NodeType1,InitPaths1] in Particles:
    for [Name2,Type2,N2,Lambda2,Statistics2,NodeType2,InitPaths2] in Particles[count:]:
      inputString += """
      Section (Observe) {string Name = "PairCorrelation"""+Name1+Name2+"""\";}
      Section (Observe) {string Name = "StructureFactor"""+Name1+Name2+"""\";}"""
      if (Statistics1 == 'FERMION' and Statistics2 == 'FERMION'):
        inputString += """
      Section (Observe) {string Name = "RefPairCorrelation"""+Name1+Name2+"""\";}
      Section (Observe) {string Name = "RefStructureFactor"""+Name1+Name2+"""\";}"""
    inputString += """
      Section (Move) {string Name = "BisectionBlock"""+Name1+"""\";}"""
    if (Statistics1 == 'FERMION'):
      inputString += """
      Section (Loop)
      {
        int Steps = """+str(N1)+""";
        bool Equilibrate = false;
        Section (Move) {string Name = "RefSlice"""+Name1+"""\";}
      }"""
    count += 1
  if CountPerms:
    inputString += """
      Section (Observe) {string Name = "CycleCount";}"""
  if TrackSign:
      inputString += """
      Section (Observe) {string Name = "Sign";}"""
  inputString += """
      Section (Observe) {string Name = "Energy";}
      Section (Observe) {string Name = "PathDump";}
      Section (Observe) {string Name = "TimeAnalysis";}
      Section (Move) {string Name = "Shift";}
    }
    Section (WriteData){}
  }

}"""

  # Write Input String to File
  f = open(fileNamePIMC+'.in','w')
  f.write(inputString)
  f.close()

def usage():
  print "Usage: python %s prefix N D rs ToTF pol nMax nTemp RefToTF RefM MinM ParticleType TrackSign restart ProcPerClone ProcPerNode" % os.path.basename(sys.argv[0])
  sys.exit(2)

def main(argv=None):
  if argv is None:
    argv = sys.argv
  if "-h" in argv or "--help" in argv:
    usage()

  print argv

  try:
    prefix = str(argv[1])
    N = int(argv[2])
    D = int(argv[3])
    rs = float(argv[4])
    ToTF = float(argv[5])
    pol = int(argv[6])
    nMax = int(argv[7])
    nTemp = int(argv[8])
    RefToTF = float(argv[9])
    RefM = int(argv[10])
    MinM = int(argv[11])
    ParticleType = int(argv[12]) # 2 - boltzmannon, 1 - fermion, 0 - boson
    TrackSign = int(argv[13])
    restart = int(argv[14])
    PPC = int(argv[15])
    PPNODE = int(argv[16])
  except:
    usage()

  print ':::Calculate Parameters:::'
  M = int((RefToTF/float(ToTF))*RefM) # Number of Time Slices
  nSquare = 14
  if M < MinM:
    nSquare -= int(log(MinM/M,2))
    M = MinM
  lam = HF.CalcLam(rs)
  eps = HF.CalcEps(rs)
  L = HF.CalcL(N,D)
  T = HF.CalcT(rs,ToTF,pol,D)
  print 'M:',M,', lam:',lam,', eps:',eps,', L:',L,', T:',T

  fileNameSquarer = prefix+'-'+str(rs)+'-'+str(RefToTF)+'-'+str(pol)+'-'+str(N)+'-'+str(RefM)+'-'+str(M)+'.sq'

  print ':::Generating POTGEN/SQUARER inputs:::'
  GenSquarerInput(fileNameSquarer,N,D,T,L,lam,eps,M,nMax,nTemp,nSquare)

  print ':::Running POTGEN:::'
  subprocess.call([os.getenv("HOME")+'/bin/potgen_lr',fileNameSquarer])

  print ':::Running SQUARER:::'
  subprocess.call([os.getenv("HOME")+'/bin/squarer',fileNameSquarer])

  print ':::Parsing Density Matrix:::'
  dmparse.Parse(fileNameSquarer+'.dm')

  print ':::Building PairAction File:::'
  BuildPairAction(fileNameSquarer,lam,D)

  print ':::Generating PIMC++ input:::'
  if not PPC:
    if M > 256:
      PPC = M/256 # Max Slice per Processor
      if PPNODE % PPC != 0 or PPC > PPNODE:
        PPC = PPNODE
    else:
      PPC = 1
  fileNamePIMC = prefix+'-'+str(rs)+'-'+str(ToTF)+'-'+str(pol)+'-'+str(N)+'-'+str(M)
  GenPIMCInput(fileNamePIMC,fileNameSquarer,N,D,T,L,lam,eps,M,nMax,pol,ParticleType,TrackSign,restart,PPC)
  g = open(fileNamePIMC+'.ins','w')
  g.write(fileNamePIMC+'.in '+str(PPC)+'\n')
  g.close()

if __name__ == "__main__":
  sys.exit(main())
