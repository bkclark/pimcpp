import os
import sys
from numpy import Inf, exp, pi, log, sqrt
#from scipy.integrate import quad
#from scipy.optimize import fsolve

def getC(rs,ToTF,pol):
  return (2./3.)*ToTF**(-3./2.)

def integrand(x,v,n):
  return (x**v)/(1+exp(x-n))

def I(v,n):
  return quad(integrand, 0, Inf, args=(v,n))[0]

def LeqR(n,v,C):
  return I(v,n)-C

def getn(C):
  return fsolve(LeqR, 1.0, args=(0.5,C), xtol=1.49012e-12)

def EFb(rs,beta,n,pol):
  if pol:
    return ((rs**3)/(3*pi)) * (1./(beta**(5./2.))) * I(3./2.,n)
  else:
    return ((2.*rs**3)/(3*pi)) * (1./(beta**(5./2.))) * I(3./2.,n)

def EF(rs,t,pol):
  C = getC(rs,t,pol)
  n = getn(C)
  beta = 1./CalcT(rs,t,pol)
  return EFb(rs,beta,n,pol)

def EF0(rs,pol):
  # (3/5)Efermi
  if pol:
    return 3.508/rs**2
  else:
    return 2.2099/rs**2

def EXintegrand(x):
  return I(-0.5,x)**2

def EXb(rs,beta,n,pol):
  if pol:
    return -(rs**3)/(6.*pi**2) * (1./beta**2) * quad(EXintegrand, -Inf, n)[0]
  else:
    return -(rs**3)/(3.*pi**2) * (1./beta**2) * quad(EXintegrand, -Inf, n)[0]

def EX(rs,t,pol):
  C = getC(rs,t,pol)
  n = getn(C)
  beta = 1./CalcT(rs,t,pol)
  return EXb(rs,beta,n,pol)

def EX0(rs,pol):
  # (3/4pi)sqrt(Efermi)
  if pol:
    return -1.1545/rs
  else:
    return -0.916331/rs

def CalcL(N,D):
  if D == 3:
    L = (4.*pi*N/3.)**(1./3.)
  elif D == 2:
    L = (pi*N)**(1./2.)
  else:
    print 'Bad dimensions', D
    sys.exit()
  return L

def CalcTF(rs,pol,D):
  if D == 3:
    # Fermi Temp (see Martin pg 103 (only he uses atomic units))
    if pol:
      TF = (9.*pi/2.)**(2./3.) / (rs**2)
    else:
      TF = (9.*pi/4.)**(2./3.) / (rs**2)
  elif D == 2:
    # Fermi Temp (see Martin pg 103 (only he uses atomic units))
    if pol:
      TF = 4./(rs**2)
    else:
      TF = 1./(rs**2)

  return TF

def CalcLam(rs):
  return 1./(rs**2)

def CalcEps(rs):
  return 2./rs

def CalcT(rs,ToTF,pol,D):
  TF = CalcTF(rs,pol,D)
  TRy = ToTF * TF # Rydberg
  return TRy

def CalcG(rs,ToTF,pol):
  T = CalcT(rs,ToTF,pol)
  return 2./(rs*T)

def CalcToTF(rs,TRy,pol):
  TF = CalcTF(rs,pol)
  ToTF = TRy/TF # Rydberg
  return ToTF

def CalcE(rs,ToTF,pol):

  TRy = CalcT(rs,ToTF,pol)
  beta = 1/TRy

  C = getC(rs,ToTF,pol)
  n = getn(C)

  ef = EFb(rs,beta,n,pol)
  ef0 = EF0(rs,pol)
  ex = EXb(rs,beta,n,pol)
  ex0 = EX0(rs,pol)
  return [ef,ef0,ex,ex0]

def usage():
  print "Usage:  %s rs ToTF pol" % os.path.basename(sys.argv[0])

def main(argv=None):
  if argv is None:
    argv = sys.argv
  if "-h" in argv or "--help" in argv:
    usage()
    sys.exit(2)

  rs = float(sys.argv[1])
  ToTF = float(sys.argv[2])
  pol = int(sys.argv[3])
  print 'rs:',rs
  print 'pol:',pol

  TRy = CalcT(rs,ToTF,pol)
  print 'T:',TRy
  print 'Beta:',1/TRy

  [ef,ef0,ex,ex0] = CalcE(rs,ToTF,pol)
  print 'EF:',ef,'EF0',ef0
  print 'EX:',ex,'EX0',ex0
  print 'E0:', ef0+ex0

if __name__ == "__main__":
  sys.exit(main())
