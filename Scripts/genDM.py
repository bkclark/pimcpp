import sys
import os
import subprocess
from math import log,sqrt,pi
from numpy import logspace, log10, loadtxt
import dmparse

def GetUnique(a):
  seen = set()
  return [x for x in a if str(x) not in seen and not seen.add(str(x))]

def GenPotgenInput(prefix,unit1,unit2,type1,type2,lam1,lam2,Z1Z2,L,D,tau,gridType,nGrid,rMin,rMax,rCut,kCut,nTemp,breakup):
    # Determine temperatures
    if nTemp > 8:
        print 'WARNING: Max nTemp is 8!'
        nTemp = 8
    minTau = tau
    maxTau = minTau*(2**(nTemp-1))
    nSquare = 14 + nTemp-1 # total number of squarings to reach lowest temperature

    print 'Creating '+prefix+'.in'
    f = open(prefix+'.in','w')
    f.write('UNITS '+unit1+' '+unit2+'\n')
    f.write('TYPE '+type1+' %f\n' % (lam1))
    f.write('TYPE '+type2+' %f\n' % (lam2))
    f.write('GRID %i %s %f %f\n' % (nGrid,gridType,rMin,rMax))
    f.write('SQUARER %f %i %i 3 30 %i\n' % (1./maxTau,nTemp,D,nSquare))
    boxString = ' '.join([str(L) for i in range(D)])
    if breakup == 2:
        f.write('POT COUL %f %f 0.D0\n' % (rMax,Z1Z2))
    if breakup == 1:
        f.write('POT COULOPT %f %f %f %i 0.D0 1.0 %s\n' % (rCut,kCut,Z1Z2,D,boxString))
    elif breakup == 0:
        f.write('POT COULLR %f %f %f %i 0.D0 1.0 %s\n' % (rCut,kCut,Z1Z2,D,boxString))
    f.close()

def RunEwald(prefix,type1,type2,lam1,lam2,Z1Z2,L,D,tau,gridType,nGrid,rMin,rMax,rCut,kCut,breakup):
    if breakup!=2:
        subprocess.call(['ewald',str(L),str(kCut),str(rMin),str(rCut),str(Z1Z2),str(nGrid),str(breakup),str(20)])
        # Write .yk file
        Vk0 = 0.0
        kData = loadtxt('kData.txt')
        ks = kData[:,0]
        Vks = kData[:,1]
        kData = GetUnique(kData.tolist())
        kData.sort()
        Vk0 = kData[0][1]
        print 'Vk0', Vk0
        f = open(prefix+'.yk','w')
        for kDatum in kData:
            f.write('  %.10E'%kDatum[0]+'       %.10E'%kDatum[1]+'\n')
        f.close()
    else:
        f = open('rData.txt','w')
        if gridType=="LINEAR":
            rs = linspace(rMin,rMax,num=nGrid,endpoint=True)
        elif gridType=="LOG":
            rs = logspace(log10(rMin),log10(rMax),num=nGrid,endpoint=True)
        else:
            print 'Unrecognized grid'
        for r in rs:
            f.write('%.10E %.10E %.10E\n' % (r,Z1Z2/r,0.0))
        f.close()

    # Write .dm and .in file for squarer
    g = open(prefix+'.dm','w')
    f = open(prefix+'.in','r')
    for line in f:
        g.write(line)
    f.close()
    if breakup!=2:
        g.write('\n VIMAGE ')
        if (D == 2):
            g.write(str(-3.90026492*Z1Z2/L))
        elif (D == 3):
            g.write(str(-2.837297479*Z1Z2/L))
    g.write('\n POTTAIL 0.0')
    g.write('\n RANK 2 '+str(nGrid)+' 1')
    g.write('\n GRID 1 '+gridType+' '+str(rMin)+' '+str(rMax))
    g.write('\n LABEL 1 r')
    g.write('\n BEGIN potential 0')
    rData = loadtxt('rData.txt')
    count = 0
    g.write('\n  ')
    for [r,VShort,VLong] in rData:
        g.write('%.10E'%(VShort)+'  ')
        count += 1
        if count % 5 == 0:
            g.write('\n  ')
    g.write('\n')
    g.close()

def GenPairActionInput(prefix,type1,lam1,type2,lam2,D):
    f = open(prefix+'.PairAction','w')
    f.write('Section (Fits)')
    f.write('\n{')
    f.write('\n  string Type="DavidFit";')
    f.write('\n  int NumOffDiagonalTerms = 3;')
    f.write('\n  Section (Particle1)')
    f.write('\n  {')
    f.write('\n    string Name = "'+type1+'";')
    f.write('\n    double lambda ='+str(lam1)+';')
    f.write('\n    int Ndim = '+str(D)+';')
    f.write('\n  }')
    f.write('\n  Section (Particle2)')
    f.write('\n  {')
    f.write('\n    string Name = "'+type2+'";')
    f.write('\n    double lambda ='+str(lam2)+';')
    f.write('\n    int Ndim = '+str(D)+';')
    f.write('\n  }')
    f.write('\n  string Daviddmfile = "inputs/'+prefix+'.h5";')
    f.write('\n}')
    f.close()

def GenDM(suffix,unit1,unit2,type1,type2,lam1,lam2,Z1Z2,L,D,tau,gridType,nGrid,rMin,rCut,nCut,nTemp,breakup,OurEwald):
    print type1, lam1, type2, lam2, 'Z1Z2', Z1Z2, 'tau', tau

    prefix = type1+'-'+type2+suffix+'.sq'

    # Box and maximum pair distance, rMax
    box = [L]*D
    rMax = 0.
    for l in box:
        rMax += l*l/4.
    rMax = sqrt(rMax)

    # Cutoff in r-space
    lam = max(lam1,lam2)
    rCut = L/2. - sqrt(lam*tau)

    # Cutoff in k-space
    kCut = nCut*2.*pi/L

    # Generate Potential
    print '**** Creating Pair Potential ****'
    GenPotgenInput(prefix,unit1,unit2,type1,type2,lam1,lam2,Z1Z2,L,D,tau,gridType,nGrid,rMin,rMax,rCut,kCut,nTemp,breakup)

    if OurEwald:
        print ':::Running Ewald:::'
        RunEwald(prefix,type1,type2,lam1,lam2,Z1Z2,L,D,tau,gridType,nGrid,rMin,rMax,rCut,kCut,breakup)
    else:
        print ':::Running POTGEN:::'
        if breakup!=2:
            subprocess.call([os.getenv("HOME")+'/bin/potgen_lr',prefix])
        else:
            subprocess.call([os.getenv("HOME")+'/bin/potgen_sr',prefix])

    # Squarer
    print '**** Performing Squaring Procedure ****'
    subprocess.call(['squarer',prefix])

    # Density Matrix Parser
    print '**** Parsing Density Matrix ****'
    dmparse.main(['',prefix+'.dm'])

    # Build PairAction File
    print '**** Creating PairAction File ****'
    GenPairActionInput(prefix,type1,lam1,type2,lam2,D)

def usage():
  print "Usage: python %s prefix unit1 unit2 type1 type2 lam1 lam2 Z1Z2 L D tau gridType nGrid rMin rCut nCut nTemp breakup" % os.path.basename(sys.argv[0])
  sys.exit(2)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    if "-h" in argv or "--help" in argv:
        usage()

    try:
        suffix = str(sys.argv[1])
        unit1 = sys.argv[2]
        unit2 = sys.argv[3]
        type1 = sys.argv[4]
        type2 = sys.argv[5]
        lam1 = float(sys.argv[6])
        lam2 = float(sys.argv[7])
        Z1Z2 = float(sys.argv[8])
        L = float(sys.argv[9])
        D = int(sys.argv[10])
        tau = float(sys.argv[11])
        gridType = sys.argv[12]
        nGrid = int(sys.argv[13])
        rMin = float(sys.argv[14])
        rCut = float(sys.argv[15])
        nCut = int(sys.argv[16])
        nTemp = int(sys.argv[17])
        breakup = int(sys.argv[18])
        OurEwald = int(sys.argv[19])
    except:
        usage()

    GenDM(suffix,unit1,unit2,type1,type2,lam1,lam2,Z1Z2,L,D,tau,gridType,nGrid,rMin,rCut,nCut,nTemp,breakup,OurEwald)

if __name__ == "__main__":
  sys.exit(main())
