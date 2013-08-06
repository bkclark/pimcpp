/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2013  B. Clark, K. Esler, E. Brown   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
// http://code.google.com/p/pimcplusplus/                  //
/////////////////////////////////////////////////////////////

#include "EAMClass.h"
#include "../PathDataClass.h"


EAMPotentialClass::EAMPotentialClass(PathDataClass &pathData) : 
  ActionBaseClass (pathData)
{
  // All parameters given in paper; see header file
  lMax = 5;
  mMax = 3;
  c.resize(6);
  c = 0.64085/12, -6.83764/12, 26.75616/12, -47.16495/12, 36.18925/12, -8.60834/12;
  //  cerr << "EAM INITIALIZED coeffs " << c << endl;
  s.resize(3);
  s = 12.0, 6.0, 24.0;
  Ec = 3.39;
  phi0 = 0.1318;

  // in angstrom
  //r0 = 2.8638;

  // converted to bohr
  r0 = 5.41361;

  alpha = 4.60;
  beta = 7.10;
  gamma = 7.34759;
  delta = 7.35;

  rn = 1.75 * r0;
  rc = 1.95 * r0;

  aob = alpha/beta;
}

void EAMPotentialClass::Read(IOSectionClass& in)
{
  cerr << "EAM READ" << endl;
  rho.resize(PathData.Path.NumTimeSlices(), PathData.Path.NumParticles());
  rho = 0.0;
  conversion = 1.0;
  cerr << "EAM default units are eV and bohr" << endl;
  if(in.ReadVar("Prefactor",conversion)) {
    cerr << "EAM using user-defined conversion factor " << conversion << endl;
  }
}


double EAMPotentialClass::SingleAction (int slice1, int slice2,
			      const Array<int,1> &changedParticles,
			      int level)
{
  double U = ComputeEnergy(slice1, slice2, changedParticles, level);
  int skip = (1<<level);
  double levelTau = Path.tau * (double)skip;
  //cerr << "EAM returning U " << U << " * tau " << levelTau << " = " << levelTau*U << " over slices " << slice1 << " " << slice2 << " for ptcls " << changedParticles << endl;
  return(levelTau*U);
}


double EAMPotentialClass::d_dBeta (int slice1, int slice2,  int level)
{
  Array<int,1> changedParticles(PathData.Path.NumParticles());
  for(int i=0; i<PathData.Path.NumParticles(); i++)
    changedParticles(i) = i;
  double U = ComputeEnergy(slice1, slice2-1, changedParticles, level);
  return(U);
}

double EAMPotentialClass::ComputeEnergy (int slice1, int slice2,
			      const Array<int,1> &changedParticles, int level)
{
  double U = 0.0;
  UpdateRho(slice1, slice2, changedParticles);

  for(int slice = slice1; slice<=slice2; slice++) {
    for (int counter=0; counter<Path.DoPtcl.size(); counter++)
      Path.DoPtcl(counter)=true;
    // SLOW!! including all terms from RHO; should be done more efficiently
    for(int ptcl1=0; ptcl1<PathData.Path.NumParticles(); ptcl1++) {
      U += conversion * F(rho(slice, ptcl1));
    }

    for(int index1=0; index1<changedParticles.size(); index1++) {
      int ptcl1 = changedParticles(index1);
      Path.DoPtcl(ptcl1) = false;

      //U += conversion * F(rho(slice, ptcl1));

      for(int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++) {
        if(ptcl2 != ptcl1 && Path.DoPtcl(ptcl2)) {
          //cerr << "EAM pairs " << ptcl1 << " " << ptcl2 << endl;
          dVec r;
          double rmag;
	      	PathData.Path.DistDisp(slice, ptcl1, ptcl2, rmag, r);
          U += conversion * 0.5 * phi(rmag);
        }
      }
    }
  }
  return (U);
}

string EAMPotentialClass::GetName()
{
  return "Al_EAM";
}

void EAMPotentialClass::UpdateRho(int slice1, int slice2, const Array<int,1>& activeP)
{
  for(int slice=slice1; slice <= slice2; slice++) {
    //for(int i=0; i<activeP.size(); i++) {
    for(int ptcl1=0; ptcl1<PathData.Path.NumParticles(); ptcl1++) {
      double newRho = 0.0;
      //int ptcl1 = activeP(i);
      for(int ptcl2=0; ptcl2<PathData.Path.NumParticles(); ptcl2++) {
        if(ptcl2 != ptcl1) {
          dVec r;
          double rmag;
	      	PathData.Path.DistDisp(slice, ptcl1, ptcl2, rmag, r);
          newRho += f(rmag);
        }
      }
      //cerr << "Updating rho for ptcl " << ptcl1 << " to " << newRho << endl;
      //cerr << "  size of rho is " << rho.size() << endl;
      //cerr << "  prev value was " << rho(ptcl1) << endl;
      rho(slice, ptcl1) = newRho;
    }
  }
}

// Smooth truncation
double EAMPotentialClass::q(double r)
{
  if(r <= rn)
    return 1.0;
  else if (r >= rc)
    return 0.0;
  else {
    double x = (r - rn)/(rc - rn);
    double qofr = (1-x)*(1-x)*(1-x)*(1 + 3*x + 6*x*x);
    return qofr;
  }
}
 
double EAMPotentialClass::f(double r)
{
  double newF = 0.0;
  double r0r = r0/r;
  for(int l=0; l<=lMax; l++) {
    newF += c(l) * pow(r0r, l);
  }
  return q(r)*newF;
}

double EAMPotentialClass::phi(double r)
{
  double rr0 = r/r0;
  double newphi = -phi0 * (1 + delta*(rr0-1)) * exp(-gamma * (rr0-1));
  return q(r)*newphi;
}

double EAMPotentialClass::F(double r)
{
  double newF = -Ec * (1 - aob * log(r)) * pow(r, aob);

  for(int m=1; m<=mMax; m++) {
    double sqM = sqrt(double(m));
    newF += 0.5*phi0 * s(m-1) * exp(-(sqM-1)*gamma)
      * (1 + (sqM-1)*delta - sqM*delta/beta*log(r))
      * pow(r,(sqM*gamma/beta));
  }

  return newF;
}
