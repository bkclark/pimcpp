import os
import sys
from math import asinh,sinh,tanh

def wt(t,w):
  return (2./t)*asinh(w*t/2.)

def z0(B,t,w,d):
  return (1./(2.*sinh(B*wt(t,w)/2.)))**(d)

def dz0(B,t,w,d):
  return -(wt(t,w)*d/2.)*(1./tanh(B*wt(t,w)/2.))*z0(B,t,w,d)

def z(N,B,t,w,d,Fermi=0):
  if N == 0:
    return 1.
  else:
    tot = 0.
    for n in range(1,N+1):
      tot += ((-1)**((n-1)*Fermi)) * z0(n*B,t,w,d) * z(N-n,B,t,w,d,Fermi)
    return tot/N

def dz(N,B,t,w,d,Fermi=0):
  if N == 0:
    return 0
  else:
    tot = 0.
    for n in range(1,N+1):
      tot += ((-1)**((n-1)*Fermi)) * (n*dz0(n*B,t,w,d)*z(N-n,B,t,w,d,Fermi) + z0(n*B,t,w,d)*dz(N-n,B,t,w,d,Fermi))
    return tot/N

def e(N,B,t,w,d,Fermi=0):
  return -(1./z(N,B,t,w,d,Fermi))*dz(N,B,t,w,d,Fermi)

def usage():
  print "Usage:  %s N B t w d (Fermi)" % os.path.basename(sys.argv[0])

def main(argv=None):
  if argv is None:
    argv = sys.argv
  if "-h" in argv or "--help" in argv:
    usage()
    sys.exit(2)

  try:
    N = int(sys.argv[1])
    B = float(sys.argv[2])
    t = float(sys.argv[3])
    w = float(sys.argv[4])
    d = int(sys.argv[5])
  except:
    usage()

  try:
    Fermi = int(sys.argv[6])
  except:
    Fermi = 0

  print 'E:', e(N,B,t,w,d,Fermi)

if __name__ == "__main__":
  sys.exit(main())

