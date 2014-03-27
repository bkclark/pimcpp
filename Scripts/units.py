import sys
import os
import math

def GetUnits(rs,ToTF,pol,N,M):
  # INPUTS
  print 'N:', N
  print 'rs:', rs
  print 'ToTF:', ToTF
  print 'M:', M
  print 'polarized:', pol

  # CONSTANTS IN UNITS (a,Ry)
  lam = 1./(rs**2)
  eps = 2./rs
  # Fermi Temp (see Martin pg 103 (only he uses atomic units))
  if pol:
    TF = (9. * math.pi / 2.)**(2./3.) / (rs**2)
  else:
    TF = (9. * math.pi / 4.)**(2./3.) / (rs**2)
  print 'lam:', lam, 'eps:', eps, 'TF:', TF

  # Length Scale Conversions
  La = (4. * math.pi * N / 3.)**(1./3.) # Wigner-Seitz radii
  a0 = 0.529 # Angstroms
  a0cm = 5.29*(10**-9)
  LAng = La*rs*a0
  Lcm = La*rs*a0cm
  print 'L:', La, 'a, ', LAng, 'Ang, ', Lcm, 'cm'

  # Density Conversions
  NoV = (3.*math.pi/4.)*((rs*a0cm)**(-3))
  mekg = 9.10938188*(10**-31)
  meg = mekg * 1000
  print 'N/V:', NoV, 'cm^-3'
  print 'm/V:', meg*NoV, 'g/cm^-3'

  # Energy Scale Conversions
  TRy = ToTF * TF # Rydberg
  beta = 1/TRy
  tau = beta/M
  print 'beta:', beta, 'Ry^-1, tau:', tau
  TH = TRy * 2. #Hartree
  RyeV = 13.6056925330
  TeV = TRy * RyeV # electron volts
  eVK = 11604.45
  TK = TeV * eVK # Kelvin
  print 'T:', TRy, 'Ry, ', TH, 'H, ', TeV, 'eV, ', TK, 'K'

  return La,TRy,tau,lam,eps

def usage():
  print "Usage:  %s rs ToTF pol N M" % os.path.basename(sys.argv[0])
  sys.exit(2)

def main(argv=None):
  if argv is None:
    argv = sys.argv
  if "-h" in argv or "--help" in argv:
    usage()

  # INPUTS
  try:
    rs = float(sys.argv[1])
    ToTF = float(sys.argv[2])
    pol = int(sys.argv[3])
    N = int(sys.argv[4])
    M = int(sys.argv[5])
  except:
    usage()

  GetUnits(rs,ToTF,pol,N,M)

if __name__ == "__main__":
  sys.exit(main())
