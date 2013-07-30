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

#include "../PathDataClass.h"
#include "HarmonicPotential.h"
#include "ctime"
#include "sys/time.h"


void HarmonicPotentialClass::Read(IOSectionClass& in)
{
  TotalTime=0;
  TimeSpent=0;
  if (!in.ReadVar ("omega",omega))
    omega = 1.0;
}


HarmonicPotentialClass::HarmonicPotentialClass(PathDataClass &pathData)
  : ActionBaseClass (pathData)
{
  omega = 1.0;
}


double HarmonicPotentialClass::dUdR(int slice, int ptcl, int level)
{
  return 0.0;
}


double HarmonicPotentialClass::d2UdR2(int slice, int ptcl1, int ptcl2, int level)
{
  assert(1==2);
  return 0.0;
}


double HarmonicPotentialClass::SingleAction (int slice1, int slice2, const Array<int,1> &changedParticles, int level)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  double TotalU = 0.0;
  PathClass &Path = PathData.Path;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau * (1<<level);

  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl = changedParticles(ptclIndex);
    double TauOmegaOverFourLambda = levelTau*omega/(4.*PathData.Path.ParticleSpecies(ptcl).lambda);
    for (int slice=slice1; slice<slice2; slice+=skip) {
      double rmag2;
      dVec r = PathData.Path(slice, ptcl);
      PathData.Path.PutInBox(r);
      PathData.Path.MagSquared(r,rmag2);
      TotalU += TauOmegaOverFourLambda * rmag2;
    }
  }
  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);
  return (TotalU);
}


double HarmonicPotentialClass::d_dBeta (int slice1, int slice2, int level)
{
  double TotalU = 0.0;
  PathClass &Path = PathData.Path;
  int skip = 1<<level;
  // double levelTau = Path.tau* (1<<level);

  for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++){
    double OmegaOverFourLambda = omega/(4.*PathData.Path.ParticleSpecies(ptcl).lambda);
    for (int slice=slice1; slice<slice2; slice+=skip){
      double rmag2;
      dVec r = PathData.Path(slice, ptcl);
      PathData.Path.PutInBox(r);
      PathData.Path.MagSquared(r,rmag2);
      TotalU += OmegaOverFourLambda * rmag2;
    }
  }
  return (TotalU);
}


string HarmonicPotentialClass::GetName()
{
  return "HarmonicPotentialClass";
}
