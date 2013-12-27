import os
import sys
from math import pi

def z0(n,lam,B,L,d):
  return (L**d)/((4.*pi*n*B*lam)**(d/2.))

def dz0(n,lam,B,L,d):
  return -2.*d*pi*n*lam*(L**d)/((4.*pi*n*B*lam)**((d/2.)+1.))

def z(N,lam,B,L,d,Fermi=0):
  if N == 0:
    return 1.
  else:
    tot = 0.
    for n in range(1,N+1):
      tot += ((-1)**((n-1)*Fermi)) * z0(n,lam,B,L,d) * z(N-n,lam,B,L,d,Fermi)
    return tot/N

def dz(N,lam,B,L,d,Fermi=0):
  if N == 0:
    return 0
  else:
    tot = 0.
    for n in range(1,N+1):
      tot += ((-1)**((n-1)*Fermi)) * (dz0(n,lam,B,L,d)*z(N-n,lam,B,L,d,Fermi) + z0(n,lam,B,L,d)*dz(N-n,lam,B,L,d,Fermi))
    return tot/N

def e(N,lam,B,L,d,Fermi=0):
  return -(1./z(N,lam,B,L,d,Fermi))*dz(N,lam,B,L,d,Fermi)

def usage():
  print "Usage:  %s N lam B L d (Fermi)" % os.path.basename(sys.argv[0])

def main(argv=None):
  if argv is None:
    argv = sys.argv
  if "-h" in argv or "--help" in argv:
    usage()
    sys.exit(2)

  try:
    N = int(sys.argv[1])
    lam = float(sys.argv[2])
    B = float(sys.argv[3])
    L = float(sys.argv[4])
    d = int(sys.argv[5])
  except:
    usage()

  try:
    Fermi = int(sys.argv[6])
  except:
    Fermi = 0

  print 'E:', e(N,lam,B,L,d,Fermi)

if __name__ == "__main__":
  sys.exit(main())

