/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
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
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "../PathDataClass.h"
#include "ExternalPotential.h"
#include "ctime"
#include "sys/time.h"

///This has to be called after pathdata knows how many
///particles it has
void ExternalPotential::Read(IOSectionClass& in)
{
  TotalTime=0;
  TimeSpent=0;
}

ExternalPotential::ExternalPotential(PathDataClass &pathData)
  : ActionBaseClass (pathData)
{
  extConst = 1.0;
}

double ExternalPotential::dUdR(int slice, int ptcl, int level)
{
  return 0.0;
}

double ExternalPotential::d2UdR2(int slice, int ptcl1, int ptcl2, int level)
{
  assert(1==2);
  return 0.0;
}

double ExternalPotential::SingleAction (int slice1, int slice2, const Array<int,1> &changedParticles, int level)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  double TotalU=0.0;
  PathClass &Path = PathData.Path;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau * (1<<level);

  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = changedParticles(ptcl1Index);
    int species1=Path.ParticleSpeciesNum(ptcl1);

    for (int slice=slice1; slice<slice2; slice+=skip){
      double rmag2;
      double U;
      dVec r = PathData.Path(slice, ptcl1);
      PathData.Path.PutInBox(r);
      PathData.Path.MagSquared(r,rmag2);
      U = 0.5 * levelTau * extConst * rmag2;
      TotalU+=U;
    }
  }
  gettimeofday(&end,   &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);
  return (TotalU);
}

double ExternalPotential::d_dBeta (int slice1, int slice2, int level)
{
  double TotalU=0.0;
  PathClass &Path = PathData.Path;
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);

  for (int ptcl1Index=0; ptcl1Index<Path.NumParticles(); ptcl1Index++){
    int ptcl1 = ptcl1Index;
    int species1=Path.ParticleSpeciesNum(ptcl1);

    for (int slice=slice1; slice<slice2; slice+=skip){
      double rmag2;
      double U;
      dVec r = PathData.Path(slice, ptcl1);
      PathData.Path.PutInBox(r);
      PathData.Path.MagSquared(r,rmag2);
      U = 0.5 * extConst * rmag2;
      TotalU+=U;
    }
  }
  return (TotalU);
}

////Gradient of the action only works in 3d!!!!!!!
//void ExternalPotential::GradAction(int slice1, int slice2, 
//			    const Array<int,1> &ptcls, int level,
//			    Array<dVec,1> &gradVec)
//{
//
//  PathClass &Path = PathData.Path;
//  int skip = (1<<level);
//  assert (gradVec.size() == ptcls.size());
//  for (int pi=0; pi<ptcls.size(); pi++) {
//    int ptcl1 = ptcls(pi);
//    int species1 = Path.ParticleSpeciesNum(ptcl1);
//    for (int ptcl2=0; ptcl2<Path.NumParticles(); ptcl2++) {
//      int species2 = Path.ParticleSpeciesNum(ptcl2);
//      PairActionFitClass &PA=*(PathData.Actions.PairMatrix(species1,species2));
//      if (ptcl1 != ptcl2) {
//	for (int slice=slice1; slice<slice2; slice += skip) {
//	  //	  double check=0.0;
//	  dVec r, rp;
//	  double rmag, rpmag, du_dq, du_dz;
//	  Path.DistDisp(slice, slice+skip, ptcl1, ptcl2, rmag, rpmag, r, rp);
//	  //	  check += dUdR(slice,ptcl1,ptcl2,level)+dUdR(slice+skip,ptcl1,ptcl2,level);
//	  double q = 0.5*(rmag+rpmag);
//	  double z = (rmag-rpmag);
//	  double s2 = dot (r-rp, r-rp);
//	  PA.Derivs(q,z,s2,level,du_dq, du_dz);
//	  dVec rhat  = (1.0/rmag)*r;
//	  dVec rphat = (1.0/rpmag)*rp;
//	  
//	  double g1 = 1.0;
//	  double g2 = 1.0;
//	  if (UseLowVariance) {
//	    g1 = g(rmag);
//	    g2 = g(rpmag);
//	  }
//	  gradVec(pi) -= (g1*(0.5*du_dq + du_dz)*rhat + 
//			  g2*(0.5*du_dq - du_dz)*rphat);
//	  //	  cerr<<"VALS: "<<(g1*(0.5*du_dq + du_dz)*rhat + 
//	  //			   g2*(0.5*du_dq - du_dz)*rphat)<<" "<<check<<endl;
//	  // gradVec(pi) -= (0.5*du_dq*(rhat+rphat) + du_dz*(rhat-rphat));
//	  /// Now, subtract off long-range part that shouldn't be in
//	  /// here 
//	  if (PA.IsLongRange() && PathData.Actions.UseLongRange)
//	    gradVec(pi) += 0.5*(PA.Ulong(level).Deriv(rmag)*g1*rhat+
//				PA.Ulong(level).Deriv(rpmag)*g2*rphat);
//	}
//      }
//    }
//  }
//
//  //  cerr<<"Gradient not working in 2d"<<endl;
//  //  assert(1==2);
//
//}


string
ExternalPotential::GetName()
{
  return "ExternalPotential";
}
