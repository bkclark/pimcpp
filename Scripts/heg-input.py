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

def BuildPairAction(prefix,lam):
  f = open(prefix+'.PairAction','w')

  f.write('Section (Fits)')
  f.write('\n{')
  f.write('\n  string Type="DavidFit";')
  f.write('\n  int NumOffDiagonalTerms = 3;')
  f.write('\n  Section (Particle1)')
  f.write('\n  {')
  f.write('\n    string Name = "e";')
  f.write('\n    double lambda ='+str(lam)+';')
  f.write('\n    int Ndim = 3;')
  f.write('\n  }')
  f.write('\n  Section (Particle2)')
  f.write('\n  {')
  f.write('\n    string Name = "e";')
  f.write('\n    double lambda ='+str(lam)+';')
  f.write('\n    int Ndim = 3;')
  f.write('\n  }')
  f.write('\n  string Daviddmfile = "inputs/' + prefix +'.h5";')
  f.write('\n}')
  f.close()

  f = open(prefix+'.PairAction','r')
  for line in f:
    print line
  f.close()

def GenPIMCInput(fileNamePIMC,fileNameSquarer,N,D,T,L,lam,eps,M,nMax,pol,ParticleType,TrackSign,NEquilibrate,restart,PPC):

  if ParticleType == 1:
    fermion = 1
    boson = 0
  else:
    fermion = 0
    if ParticleType == 0:
      boson = 1
    else:
      boson = 0 # meaning boltzmannon

  if TrackSign and fermion:
    CountPerms = 0
  elif TrackSign and not fermion:
    CountPerms = 1
  else:
    CountPerms = 0

  MaxClones = 0 # HACKED VALUE !!!

  NUp = N
  if not pol:
    assert(not N%2)
    NUp = N/2
    NDown = N/2

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

  f = open(fileNamePIMC+'.in','w')
  f.write('bool OpenLoops = false;\n')
  f.write('int NEquilibrate = '+str(NEquilibrate)+';\n')
  f.write('\n')
  f.write('Section (Output)\n')
  f.write('{\n')
  f.write('  string OutFileBase = "rawdata/'+fileNamePIMC+'";\n')
  f.write('  bool Restart = true;\n')
  f.write('}\n')
  f.write('\n')
  f.write('Section (Parallel)\n')
  f.write('{\n')
  f.write('  int ProcsPerClone = '+str(PPC)+';\n')
  f.write('}\n')
  f.write('\n')
  f.write('Section (System)\n')
  f.write('{\n')
  f.write('  double tau = '+str(tau)+';\n')
  f.write('  int NumTimeSlices = '+str(M)+';\n')
  f.write('  Array<double,1> Box('+str(D)+') = [')
  for i in range(0,D-1):
    f.write(str(L)+',')
  f.write(str(L)+'];\n')
  f.write('  Array<bool,1> IsPeriodic('+str(D)+') = [')
  for i in range(0,D-1):
    f.write('true,')
  f.write('true];\n')
  f.write('  Section (Particles)\n')
  f.write('  {\n')
  f.write('    Section (Species)\n')
  f.write('    {\n')
  f.write('      string Name = "eUp";\n')
  f.write('      string Type = "e";\n')
  f.write('      double lambda = '+str(lam)+';\n')
  if fermion:
    f.write('      string Statistics = "FERMION";\n')
    f.write('      string NodeType = "FREE";\n')
  elif boson:
    f.write('      string Statistics = "BOSON";\n')
  else:
    f.write('      string Statistics = "BOLTZMANNON";\n')
  f.write('      int NumParticles = '+str(NUp)+';\n')
  f.write('      int NumDim = '+str(D)+';\n')
  if restart:
    f.write('      string InitPaths = "RESTART";\n')
    f.write('      int MaxClones = '+str(MaxClones)+';\n')
    #f.write('      bool ParallelFileRead = false;\n')
  else:
    f.write('      string InitPaths = "BCC";\n')
  f.write('      string File = "rawdata/'+fileNamePIMC+'";\n')
  f.write('    }\n')
  if not pol:
    f.write('    Section (Species)\n')
    f.write('    {\n')
    f.write('      string Name = "eDown";\n')
    f.write('      string Type = "e";\n')
    f.write('      double lambda = '+str(lam)+';\n')
    if fermion:
      f.write('      string Statistics = "FERMION";\n')
      f.write('      string NodeType = "FREE";\n')
    else:
      f.write('      string Statistics = "BOSON";\n')
    f.write('      int NumParticles = '+str(NDown)+';\n')
    f.write('      int NumDim = '+str(D)+';\n')
    if restart:
      f.write('      string InitPaths = "RESTART";\n')
      f.write('      int MaxClones = '+str(MaxClones)+';\n')
      #f.write('      bool ParallelFileRead = false;\n')
    else:
      f.write('      string InitPaths = "BCC";\n')
    f.write('      string File ="rawdata/'+fileNamePIMC+'";\n')
    f.write('    }\n')
  f.write('  }\n')
  f.write('}\n')
  f.write('\n')
  f.write('Section (Actions)\n')
  f.write('{\n')
  f.write('  int NumImages = 1;\n')
  f.write('  int MaxLevels = '+str(MaxLevels)+';\n')
  f.write('  Array<string,1> PairActionFiles(1) = ["inputs/'+fileNameSquarer+'.PairAction"];\n')
  f.write('\n')
  f.write('  Section (Action)\n')
  f.write('  {\n')
  f.write('    string Name = "ShortRange";\n')
  f.write('    string Type = "ShortRange";\n')
  f.write('  }\n')
  f.write('\n')
  f.write('  Section (Action)\n')
  f.write('  {\n')
  f.write('    string Name = "DavidLongRange";\n')
  f.write('    string Type = "DavidLongRange";\n')
  f.write('    double kCutoff = '+str(kMax)+';\n')
  f.write('  }\n')
  if fermion:
    f.write('  Section (NodalAction)\n')
    f.write('  {\n')
    f.write('    string Type = "FREE";\n')
    f.write('    string Species = "eUp";\n')
    f.write('    bool StoreNodeDist = true;\n')
    f.write('    bool StoreNodeDet = true;\n')
    f.write('    bool UseNoDist = true;\n')
    f.write('  }\n')
    f.write('\n')
    if not pol:
      f.write('  Section (NodalAction)\n')
      f.write('  {\n')
      f.write('    string Type = "FREE";\n')
      f.write('    string Species = "eDown";\n')
      f.write('    bool StoreNodeDist = true;\n')
      f.write('    bool StoreNodeDet = true;\n')
      f.write('    bool UseNoDist = true;\n')
      f.write('  }\n')
      f.write('\n')
  f.write('}\n')
  f.write('\n')
  f.write('Section (Observables)\n')
  f.write('{\n')
  f.write('  Section (Observable)\n')
  f.write('  {\n')
  f.write('    string Type = "Energy";\n')
  f.write('    string Name = "Energy";\n')
  f.write('    string Description = "Total Energy";\n')
  f.write('    int Frequency = 1;\n')
  f.write('    double HistStart = 0.0;\n')
  f.write('    double HistEnd = 1.0;\n')
  f.write('    int HistPoints = 2;\n')
  if TrackSign:
    f.write('    bool TrackSign = true;\n')
  if CountPerms:
    f.write('    bool CountPerms = true;\n')
  f.write('  }\n')
  f.write('\n')
  if CountPerms:
    f.write('  Section (Observable)\n')
    f.write('  {\n')
    f.write('    string Type = "CycleCount";\n')
    f.write('    string Name = "CycleCount";\n')
    f.write('    int Frequency = 1;\n')
    f.write('  }\n')
    f.write('\n')
  if TrackSign:
    f.write('  Section (Observable)\n')
    f.write('  {\n')
    f.write('    string Type = "Sign";\n')
    f.write('    string Name = "Sign";\n')
    f.write('    string Description = "Total Energy";\n')
    f.write('    int Frequency = 1;\n')
    f.write('  }\n')
    f.write('\n')
  f.write('  Section (Observable)\n')
  f.write('  {\n')
  f.write('    string Type = "PairCorrelation";\n')
  f.write('    string Name = "PairCorrelationUpUp";\n')
  f.write('    string Description = "Pair Correlation";\n')
  f.write('    string Species1 = "eUp";\n')
  f.write('    string Species2 = "eUp";\n')
  f.write('    Section (Grid)\n')
  f.write('    {\n')
  f.write('      string Type = "Linear";\n')
  f.write('      double start = 0.0;\n')
  f.write('      int NumPoints = 100;\n')
  f.write('    }\n')
  f.write('    int Frequency = 1;\n')
  f.write('  }\n')
  f.write('\n')
  if fermion:
    f.write('  Section (Observable)\n')
    f.write('  {\n')
    f.write('    string Type = "RefPairCorrelation";\n')
    f.write('    string Name = "RefPairCorrelationUpUp";\n')
    f.write('    string Description = "Ref Pair Correlation";\n')
    f.write('    string Species1 = "eUp";\n')
    f.write('    string Species2 = "eUp";\n')
    f.write('    Section (Grid)\n')
    f.write('    {\n')
    f.write('      string Type = "Linear";\n')
    f.write('      double start = 0.0;\n')
    f.write('      int NumPoints = 100;\n')
    f.write('    }\n')
    f.write('    int Frequency = 1;\n')
    f.write('  }\n')
    f.write('\n')
  if not pol:
    f.write('  Section (Observable)\n')
    f.write('  {\n')
    f.write('    string Type = "PairCorrelation";\n')
    f.write('    string Name = "PairCorrelationUpDown";\n')
    f.write('    string Description = "Pair Correlation";\n')
    f.write('    string Species1 = "eUp";\n')
    f.write('    string Species2 = "eDown";\n')
    f.write('    Section (Grid)\n')
    f.write('    {\n')
    f.write('      string Type = "Linear";\n')
    f.write('      double start = 0.0;\n')
    f.write('      int NumPoints = 100;\n')
    f.write('    }\n')
    f.write('    int Frequency = 1;\n')
    f.write('  }\n')
    f.write('\n')
    f.write('  Section (Observable)\n')
    f.write('  {\n')
    f.write('    string Type = "PairCorrelation";\n')
    f.write('    string Name = "PairCorrelationDownDown";\n')
    f.write('    string Description = "Pair Correlation";\n')
    f.write('    string Species1 = "eDown";\n')
    f.write('    string Species2 = "eDown";\n')
    f.write('    Section (Grid)\n')
    f.write('    {\n')
    f.write('      string Type = "Linear";\n')
    f.write('      double start = 0.0;\n')
    f.write('      int NumPoints = 100;\n')
    f.write('    }\n')
    f.write('    int Frequency = 1;\n')
    f.write('  }\n')
    f.write('\n')
    if fermion:
      f.write('  Section (Observable)\n')
      f.write('  {\n')
      f.write('    string Type = "RefPairCorrelation";\n')
      f.write('    string Name = "RefPairCorrelationUpDown";\n')
      f.write('    string Description = "Ref Pair Correlation";\n')
      f.write('    string Species1 = "eUp";\n')
      f.write('    string Species2 = "eDown";\n')
      f.write('    Section (Grid)\n')
      f.write('    {\n')
      f.write('      string Type = "Linear";\n')
      f.write('      double start = 0.0;\n')
      f.write('      int NumPoints = 100;\n')
      f.write('    }\n')
      f.write('    int Frequency = 1;\n')
      f.write('  }\n')
      f.write('\n')
      f.write('  Section (Observable)\n')
      f.write('  {\n')
      f.write('    string Type = "RefPairCorrelation";\n')
      f.write('    string Name = "RefPairCorrelationDownDown";\n')
      f.write('    string Description = "Ref Pair Correlation";\n')
      f.write('    string Species1 = "eDown";\n')
      f.write('    string Species2 = "eDown";\n')
      f.write('    Section (Grid)\n')
      f.write('    {\n')
      f.write('      string Type = "Linear";\n')
      f.write('      double start = 0.0;\n')
      f.write('      int NumPoints = 100;\n')
      f.write('    }\n')
      f.write('    int Frequency = 1;\n')
      f.write('  }\n')
      f.write('\n')
  f.write('  Section (Observable)\n')
  f.write('  {\n')
  f.write('    bool AllClones = true;\n')
  f.write('    string Type = "PathDump";\n')
  f.write('    string Name = "PathDump";\n')
  f.write('    int TemporalFrequency = '+str(1700)+';\n')
  f.write('  }\n')
  f.write('\n')
  f.write('  Section (Observable)\n')
  f.write('  {\n')
  f.write('    string Type = "StructureFactor";\n')
  f.write('    string Name = "StructureFactorUpUp";\n')
  f.write('    string Description = "Structure Factor";\n')
  f.write('    double kCutoff = '+str(kMax)+';\n')
  f.write('    string Species1 = "eUp";\n')
  f.write('    string Species2 = "eUp";\n')
  f.write('    int Frequency = 1;\n')
  f.write('  }\n')
  f.write('\n')
  if fermion:
    f.write('  Section (Observable)\n')
    f.write('  {\n')
    f.write('    string Type = "RefStructureFactor";\n')
    f.write('    string Name = "RefStructureFactorUpUp";\n')
    f.write('    string Description = "Ref Structure Factor";\n')
    f.write('    double kCutoff = '+str(kMax)+';\n')
    f.write('    string Species1 = "eUp";\n')
    f.write('    string Species2 = "eUp";\n')
    f.write('    int Frequency = 1;\n')
    f.write('  }\n')
    f.write('\n')
  if not pol:
    f.write('  Section (Observable)\n')
    f.write('  {\n')
    f.write('    string Type = "StructureFactor";\n')
    f.write('    string Name = "StructureFactorUpDown";\n')
    f.write('    string Description = "Structure Factor";\n')
    f.write('    double kCutoff = '+str(kMax)+';\n')
    f.write('    string Species1 = "eUp";\n')
    f.write('    string Species2 = "eDown";\n')
    f.write('    int Frequency = 1;\n')
    f.write('  }\n')
    f.write('\n')
    f.write('  Section (Observable)\n')
    f.write('  {\n')
    f.write('    string Type = "StructureFactor";\n')
    f.write('    string Name = "StructureFactorDownDown";\n')
    f.write('    string Description = "Structure Factor";\n')
    f.write('    double kCutoff = '+str(kMax)+';\n')
    f.write('    string Species1 = "eDown";\n')
    f.write('    string Species2 = "eDown";\n')
    f.write('    int Frequency = 1;\n')
    f.write('  }\n')
    f.write('\n')
    if fermion:
      f.write('  Section (Observable)\n')
      f.write('  {\n')
      f.write('    string Type = "RefStructureFactor";\n')
      f.write('    string Name = "RefStructureFactorUpDown";\n')
      f.write('    string Description = "Ref Structure Factor";\n')
      f.write('    double kCutoff = '+str(kMax)+';\n')
      f.write('    string Species1 = "eUp";\n')
      f.write('    string Species2 = "eDown";\n')
      f.write('    int Frequency = 1;\n')
      f.write('  }\n')
      f.write('\n')
      f.write('  Section (Observable)\n')
      f.write('  {\n')
      f.write('    string Type = "RefStructureFactor";\n')
      f.write('    string Name = "RefStructureFactorDownDown";\n')
      f.write('    string Description = "Ref Structure Factor";\n')
      f.write('    double kCutoff = '+str(kMax)+';\n')
      f.write('    string Species1 = "eDown";\n')
      f.write('    string Species2 = "eDown";\n')
      f.write('    int Frequency = 1;\n')
      f.write('  }\n')
      f.write('\n')
  f.write('  Section (Observable)\n')
  f.write('  {\n')
  f.write('    string Type = "TimeAnalysis";\n')
  f.write('    string Name = "TimeAnalysis";\n')
  f.write('    int Frequency = 1;\n')
  f.write('  }\n')
  f.write('\n')
  f.write('}\n')
  f.write('\n')
  f.write('Section (Moves)\n')
  f.write('{\n')
  f.write('\n')
  f.write('  Section (Move)\n')
  f.write('  {\n')
  f.write('    string Type = "Displace";\n')
  f.write('    string Name = "DisplaceUp";\n')
  f.write('    double Sigma = '+str(L/8.)+';\n')
  f.write('    Array<string,1> ActiveSpecies(1) = ["eUp"];\n')
  f.write('    int NumToMove = 1;\n')
  f.write('    bool MoveAll = true;\n')
  f.write('    Array<string,1> SamplingActions(2) = ["ShortRange","DavidLongRange"];\n')
  f.write('  }\n')
  f.write('\n')
  if not pol:
    f.write('  Section (Move)\n')
    f.write('  {\n')
    f.write('    string Type = "Displace";\n')
    f.write('    string Name = "DisplaceDown";\n')
    f.write('    double Sigma = '+str(sigma)+';\n')
    f.write('    Array<string,1> ActiveSpecies(1) = ["eDown"];\n')
    f.write('    int NumToMove = 1;\n')
    f.write('    bool MoveAll = true;\n')
    f.write('    Array<string,1> SamplingActions(2) = ["ShortRange","DavidLongRange"];\n')
    f.write('  }\n')
    f.write('\n')
  f.write('  Section (Move)\n')
  f.write('  {\n')
  f.write('    string Type = "BisectionBlock";\n')
  f.write('    string Name = "BisectionBlockUp";\n')
  f.write('    //Array<string,1> HigherLevelActions(2) = ["ShortRange","DavidLongRange"];\n')
  f.write('    Array<string,1> SamplingActions(2) = ["ShortRange","DavidLongRange"];\n')
  if boson or fermion:
    f.write('    string PermuteType = "TABLE";\n')
    if fermion:
      f.write('    Array<double,1> Gamma(4) = [1.0,0.0,1.0,0.0];\n')
    else:
      f.write('    Array<double,1> Gamma(4) = [1.0,1.0,1.0,1.0];\n')
  else:
    f.write('    string PermuteType = "NONE";\n')
  f.write('    double epsilon = 1e-5;\n')
  f.write('    string Species = "eUp";\n')
  f.write('    int NumLevels = '+str(MaxLevels)+';\n')
  f.write('    int StepsPerBlock = '+str(NUp)+';\n')
  f.write('  }\n')
  f.write('\n')
  if not pol:
    f.write('  Section (Move)\n')
    f.write('  {\n')
    f.write('    string Type = "BisectionBlock";\n')
    f.write('    string Name = "BisectionBlockDown";\n')
    f.write('    //Array<string,1> HigherLevelActions(2) = ["ShortRange","DavidLongRange"];\n')
    f.write('    Array<string,1> SamplingActions(2) = ["ShortRange","DavidLongRange"];\n')
    if boson or fermion:
      f.write('    string PermuteType = "TABLE";\n')
      if fermion:
        f.write('    Array<double,1> Gamma(4) = [1.0,0.0,1.0,0.0];\n')
      else:
        f.write('    Array<double,1> Gamma(4) = [1.0,1.0,1.0,1.0];\n')
    else:
      f.write('    string PermuteType = "NONE";\n')
    f.write('    double epsilon = 1e-5;\n')
    f.write('    string Species = "eDown";\n')
    f.write('    int NumLevels = '+str(MaxLevels)+';\n')
    f.write('    int StepsPerBlock = '+str(NDown)+';\n')
    f.write('  }\n')
    f.write('\n')
  f.write('  Section (Move)\n')
  f.write('  {\n')
  f.write('    string Type = "ShiftMove";\n')
  f.write('    string Name = "Shift";\n')
  f.write('  }\n')
  f.write('\n')
  if fermion:
    f.write('  Section (Move)\n')
    f.write('  {\n')
    f.write('    string Type = "RefSlice";\n')
    f.write('    string Name = "RefSliceUp";\n')
    f.write('    string PermuteType = "NONE";\n')
    f.write('    int NumLevels = '+str(MaxLevels)+';\n')
    f.write('    string Species = "eUp";\n')
    f.write('    Array<string,1> SamplingActions(2) = ["ShortRange","DavidLongRange"];\n')
    f.write('  }\n')
    f.write('\n')
    if not pol:
      f.write('  Section (Move)\n')
      f.write('  {\n')
      f.write('    string Type = "RefSlice";\n')
      f.write('    string Name = "RefSliceDown";\n')
      f.write('    string PermuteType = "NONE";\n')
      f.write('    int NumLevels = '+str(MaxLevels)+';\n')
      f.write('    string Species = "eDown";\n')
      f.write('    Array<string,1> SamplingActions(2) = ["ShortRange","DavidLongRange"];\n')
      f.write('  }\n')
      f.write('\n')
  f.write('}\n')
  f.write('\n')
  f.write('Section (Algorithm)\n')
  f.write('{\n')
  f.write('\n')
  f.write('  Section (Loop){\n')
  f.write('    int Steps = '+str(NEquilibrate)+';\n')
  f.write('    bool Equilibrate = true;\n')
  f.write('\n')
  f.write('    Section (Move) {string Name = "DisplaceUp";}\n')
  if not pol:
    f.write('    Section (Move) {string Name = "DisplaceDown";}\n')
  f.write('    Section (Move) {string Name = "BisectionBlockUp";}\n')
  if not pol:
    f.write('    Section (Move) {string Name = "BisectionBlockDown";}\n')
  if fermion:
    f.write('    Section (Loop) {\n')
    f.write('      int Steps = '+str(NUp)+';\n')
    f.write('      bool Equilibrate = true;\n')
    f.write('      Section (Move) {string Name = "RefSliceUp";}\n')
    if not pol:
      f.write('      Section (Move) {string Name = "RefSliceDown";}\n')
    f.write('    }\n')
  f.write('    Section (Observe) {string Name = "TimeAnalysis";}\n')
  f.write('    Section (Move) {string Name = "Shift";}\n')
  f.write('  }\n')
  f.write('\n')
  f.write('  Section (Loop){\n')
  f.write('    int Steps = 10000000;\n')
  f.write('\n')
  f.write('    Section (Loop){\n')
  if fermion:
    f.write('      int Steps = 10;\n')
  else:
    f.write('      int Steps = 1000;\n')
  f.write('      bool Equilibrate = false;\n')
  f.write('\n')
  f.write('      Section (Observe) {string Name = "Energy";}\n')
  f.write('      Section (Observe) {string Name = "PairCorrelationUpUp";}\n')
  if fermion:
    f.write('      Section (Observe) {string Name = "RefPairCorrelationUpUp";}\n')
  if not pol:
    f.write('      Section (Observe) {string Name = "PairCorrelationUpDown";}\n')
    f.write('      Section (Observe) {string Name = "PairCorrelationDownDown";}\n')
    if fermion:
      f.write('      Section (Observe) {string Name = "RefPairCorrelationUpDown";}\n')
      f.write('      Section (Observe) {string Name = "RefPairCorrelationDownDown";}\n')
  f.write('      Section (Observe) {string Name = "StructureFactorUpUp";}\n')
  if fermion:
    f.write('      Section (Observe) {string Name = "RefStructureFactorUpUp";}\n')
  if not pol:
    f.write('      Section (Observe) {string Name = "StructureFactorUpDown";}\n')
    f.write('      Section (Observe) {string Name = "StructureFactorDownDown";}\n')
    if fermion:
      f.write('      Section (Observe) {string Name = "RefStructureFactorUpDown";}\n')
      f.write('      Section (Observe) {string Name = "RefStructureFactorDownDown";}\n')
  f.write('      Section (Move) {string Name = "DisplaceUp";}\n')
  if not pol:
    f.write('      Section (Move) {string Name = "DisplaceDown";}\n')
  f.write('      Section (Move) {string Name = "BisectionBlockUp";}\n')
  if not pol:
    f.write('      Section (Move) {string Name = "BisectionBlockDown";}\n')
  if fermion:
    f.write('      Section (Loop) {\n')
    f.write('        int Steps = '+str(NUp)+';\n')
    f.write('        bool Equilibrate = false;\n')
    f.write('        Section (Move) {string Name = "RefSliceUp";}\n')
    if not pol:
      f.write('        Section (Move) {string Name = "RefSliceDown";}\n')
    f.write('      }\n')
  if CountPerms:
    f.write('      Section (Observe) {string Name = "CycleCount";}\n')
  if TrackSign:
    f.write('      Section (Observe) {string Name = "Sign";}\n')
  f.write('      Section (Observe) {string Name = "PathDump";}\n')
  f.write('      Section (Observe) {string Name = "TimeAnalysis";}\n')
  f.write('      Section (Move) {string Name = "Shift";}\n')
  f.write('    }\n')
  f.write('    Section (WriteData){}\n')
  f.write('  }\n')
  f.write('\n')
  f.write('}\n')
  f.close()

def usage():
  print "Usage: python %s prefix N D rs ToTF pol nMax nTemp RefToTF RefM MinM ParticleType TrackSign NEquilibrate restart ProcPerClone ProcPerNode" % os.path.basename(sys.argv[0])
  sys.exit(2)

def main(argv=None):
  if argv is None:
    argv = sys.argv
  if "-h" in argv or "--help" in argv:
    usage()

  try:
    prefix = str(sys.argv[1])
    N = int(sys.argv[2])
    D = int(sys.argv[3])
    rs = float(sys.argv[4])
    ToTF = float(sys.argv[5])
    pol = int(sys.argv[6])
    nMax = int(sys.argv[7])
    nTemp = int(sys.argv[8])
    RefToTF = float(sys.argv[9])
    RefM = int(sys.argv[10])
    MinM = int(sys.argv[11])
    ParticleType = int(sys.argv[12]) # 2 - boltzmannon, 1 - fermion, 0 - boson
    TrackSign = int(sys.argv[13])
    NEquilibrate = int(sys.argv[14])
    restart = int(sys.argv[15])
    PPC = int(sys.argv[16])
    PPNODE = int(sys.argv[17])
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
  BuildPairAction(fileNameSquarer,lam)

  print ':::Generating PIMC++ input:::'
  if not PPC:
    if M > 256:
      PPC = M/256 # Max Slice per Processor
      if PPNODE % PPC != 0 or PPC > PPNODE:
        PPC = PPNODE
    else:
      PPC = 1
  fileNamePIMC = prefix+'-'+str(rs)+'-'+str(ToTF)+'-'+str(pol)+'-'+str(N)+'-'+str(M)
  GenPIMCInput(fileNamePIMC,fileNameSquarer,N,D,T,L,lam,eps,M,nMax,pol,ParticleType,TrackSign,NEquilibrate,restart,PPC)
  g = open(fileNamePIMC+'.ins','w')
  g.write(fileNamePIMC+'.in '+str(PPC)+'\n')
  g.close()

if __name__ == "__main__":
  sys.exit(main())
