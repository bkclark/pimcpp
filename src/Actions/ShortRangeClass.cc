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
#include "ShortRangeClass.h"
#include "ShortRangeOnClass.h"
#include "ctime"
#include "../MatrixOps/MatrixOps.h"
#include "sys/time.h"


///This has to be called after pathdata knows how many particles it has
void ShortRangeClass::Read(IOSectionClass& in)
{
  TotalTime = 0;
  TimeSpent = 0;
}


ShortRangeClass::ShortRangeClass(PathDataClass &pathData, Array<PairActionFitClass* ,2> &pairMatrix) :
  ToCheck(pathData,pairMatrix), ActionBaseClass (pathData), PairMatrix(pairMatrix),
  m(2), NumBasisFuncs(4), Router(0.6), UseLowVariance(/*true*/false)
{
  Setup_ck();
}


void ShortRangeClass::Setup_ck()
{
  ck.resize(NumBasisFuncs);
  Array<double,1> h(NumBasisFuncs);
  Array<double,2> S(NumBasisFuncs);
  double Rtojplus1 = Router*Router;
  for (int j=1; j<=NumBasisFuncs; j++) {
    h(j-1) = Rtojplus1/(1.0+j);
    Rtojplus1 *= Router;
    for (int k=1; k<=NumBasisFuncs; k++)
      S(k-1,j-1) = pow(Router, m+k+j+1)/(m+k+j+1.0);
  }
  GJInverse (S);
  ck = 0.0;
  for (int i=0; i < NumBasisFuncs; i++)
    for (int j=0; j < NumBasisFuncs; j++)
      ck(i) += S(i,j) * h(j);
}


inline double ShortRangeClass::g(double r)
{
  if (r > Router)
    return 1.0;
  double rtokplusm = r;
  for (int i=0; i<m; i++)
    rtokplusm *= r;

  double gval = 0.0;
  for (int k=1; k<=NumBasisFuncs; k++) {
    gval += ck(k-1) * rtokplusm;
    rtokplusm *= r;
  }
  return gval;
}


double ShortRangeClass::dUdR(int slice,int ptcl1, int ptcl2, int level)
{
  int species1=Path.ParticleSpeciesNum(ptcl1);
  int species2=Path.ParticleSpeciesNum(ptcl2);
  PairActionFitClass *PA = PairMatrix(species1, species2);
  assert(((DavidPAClass*)PA)->SamplingTableRead);
  double dist;
  dVec disp;
  PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
  return ((DavidPAClass*)PA)->dUdRTimesSigma(dist,level);
  
}


double ShortRangeClass::dUdR_movers(int slice,int ptcl1, int ptcl2, int level)
{
  int species1=Path.ParticleSpeciesNum(ptcl1);
  int species2=Path.ParticleSpeciesNum(ptcl2);
  PairActionFitClass *PA = PairMatrix(species1, species2);
  assert(((DavidPAClass*)PA)->SamplingTableRead);
  double dist;
  dVec disp;
  PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
  return ((DavidPAClass*)PA)->dUdRTimesSigma_movers(dist,level);
}


double ShortRangeClass::d2UdR2(int slice,int ptcl1, int ptcl2, int level)
{
  int species1=Path.ParticleSpeciesNum(ptcl1);
  int species2=Path.ParticleSpeciesNum(ptcl2);
  PairActionFitClass *PA = PairMatrix(species1, species2);
  assert(((DavidPAClass*)PA)->SamplingTableRead);
  double dist;
  dVec disp;
  PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
  return ((DavidPAClass*)PA)->dUdRTimesSigma(dist,level);
}


double ShortRangeClass::d2UdR2_movers(int slice,int ptcl1, int ptcl2, int level)
{
  int species1=Path.ParticleSpeciesNum(ptcl1);
  int species2=Path.ParticleSpeciesNum(ptcl2);
  PairActionFitClass *PA = PairMatrix(species1, species2);
  assert(((DavidPAClass*)PA)->SamplingTableRead);
  double dist;
  dVec disp;
  PathData.Path.DistDisp(slice,ptcl1,ptcl2,dist,disp);
  return ((DavidPAClass*)PA)->dUdRTimesSigma_movers(dist,level);
}


double ShortRangeClass::SingleAction (int slice1, int slice2, const Array<int,1> &changedParticles, int level)
{
  struct timeval start, end;
  struct timezone tz;
  double TotalU = 0.0;
  PathClass &Path = PathData.Path;
  // First, sum the pair actions
  Path.DoPtcl = true;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau * (1<<level);
  gettimeofday(&start, &tz);
  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++) {
    int ptcl1 = changedParticles(ptcl1Index);
    Path.DoPtcl(ptcl1) = false;
    int species1 = Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0; ptcl2<Path.NumParticles(); ptcl2++) {
      if (Path.DoPtcl(ptcl2)) {
        int species2 = Path.ParticleSpeciesNum(ptcl2);
        PairActionFitClass &PA = *(PairMatrix(species1, species2));
        for (int slice=slice1; slice<slice2; slice+=skip) {
          dVec r, rp;
          double rmag, rpmag;
          PathData.Path.DistDisp(slice, slice+skip, ptcl1, ptcl2, rmag, rpmag, r, rp);
          double s2 = dot(r-rp, r-rp);
          double q = 0.5*(rmag + rpmag);
          double z = (rmag - rpmag);
          double U;
          U = PA.U(q,z,s2,level);
          if(isnan(U)) {
            cerr << "U "<< U << " " << q << " " << z << " " << s2 << " " << level << endl;
            abort();
          }
          // double UP = 0.5*(PA.U(rmag,0,0, level)+PA.U(rpmag,0,0,level))
          /// Subtract off long-range part from short-range action
          if (PA.IsLongRange() && PathData.Actions.UseLongRange)
            U -= 0.5* (PA.Ulong(level)(rmag) + PA.Ulong(level)(rpmag));
          TotalU+=U;
          //  cerr<<" "<<U<<" "<<-0.5*levelTau*(1.0/rmag+1.0/rpmag)<<endl;
          //  TotalU+=-0.5*levelTau*(1.0/rmag+1.0/rpmag);
        }
      }
    }
  }
  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);
  return (TotalU);
}


//Doesn't deal correctly with the long range term possibly
double ShortRangeClass::SingleActionForcedPairAction (int slice1, int slice2, PairActionFitClass &PA)
{
  double TotalU=0.0;
  PathClass &Path = PathData.Path;
  // First, sum the pair actions
  for (int ptcl=0;ptcl<Path.DoPtcl.size();ptcl++)
    Path.DoPtcl(ptcl)=true;
  int level=0;
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptcl1=0;ptcl1<Path.NumParticles();ptcl1++){
    int species1=Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0;ptcl2<ptcl1;ptcl2++) {
      int species2=Path.ParticleSpeciesNum(ptcl2);
      for (int slice=slice1;slice<slice2;slice+=skip){
	dVec r, rp;
	double rmag, rpmag;
	PathData.Path.DistDisp(slice, slice+skip, ptcl1, ptcl2,
			       rmag, rpmag, r, rp);
	double s2 = dot (r-rp, r-rp);
	double q = 0.5 * (rmag + rpmag);
	double z = (rmag - rpmag);
	double U;
	U = PA.U(q,z,s2, level);
	// Subtract off long-range part from short-range action
	if (PA.IsLongRange() && PathData.Actions.UseLongRange)
	  U -= 0.5* (PA.Ulong(level)(rmag) + PA.Ulong(level)(rpmag));
	TotalU += U;
      }
    }
  }
  return (TotalU);
}


double ShortRangeClass::d_dBetaForcedPairAction (int slice1, int slice2, PairActionFitClass &pa)
{
  cerr<<"Calling me "<<endl;
  double levelTau=Path.tau;
  int level=0;
  int skip = 1<<level;
  // Add constant part.  Note: we should really check the number of
  // dimensions. 
  double dU = 0.0;
  for (int ptcl1=0; ptcl1<PathData.NumParticles(); ptcl1++) {
    int species1=Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
      for (int slice=slice1;slice<slice2;slice+=skip){
	dVec r, rp;
	double rmag, rpmag;
	PathData.Path.DistDisp(slice,slice+skip,ptcl1,ptcl2,rmag,rpmag,r,rp);
	
	double s2 = dot(r-rp, r-rp);
	double q = 0.5*(rmag+rpmag);
	double z = (rmag-rpmag);
	

	dU += pa.dU(q, z, s2, level)*
	  PathData.Path.ParticleExist(slice,ptcl1)*
	  PathData.Path.ParticleExist(slice,ptcl2)*
	  PathData.Path.ParticleExist(slice+skip,ptcl1)*
	  PathData.Path.ParticleExist(slice+skip,ptcl2);

	// Subtract off long-range part from short-range action
	if (pa.IsLongRange() && PathData.Actions.UseLongRange)
	  dU -= 0.5*(pa.dUlong(level)(rmag)+pa.dUlong(level)(rpmag));
      }
    }
  }
  return dU;
}


double ShortRangeClass::d_dBeta (int slice1, int slice2, int level)
{
  PathClass &Path = PathData.Path;

  double levelTau=Path.tau;
  int skip = 1<<level;
  // Add constant part.  Note: we should really check the number of dimensions.
  double dU = 0.0;
  for (int ptcl1=0; ptcl1<PathData.NumParticles(); ptcl1++) {
    int species1=Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
      int species2 = Path.ParticleSpeciesNum(ptcl2);
      PairActionFitClass& pa = *(PairMatrix(species1,species2));
      for (int slice=slice1; slice<slice2; slice+=skip){
        dVec r, rp;
        double rmag, rpmag;
        Path.DistDisp(slice,slice+skip,ptcl1,ptcl2,rmag,rpmag,r,rp);

        double s2 = dot(r-rp, r-rp);
        double q = 0.5*(rmag+rpmag);
        double z = (rmag-rpmag);
        if (Path.WormOn){
          dU += pa.dU(q,z,s2,level)*
            Path.ParticleExist(slice,ptcl1)*
            Path.ParticleExist(slice,ptcl2)*
            Path.ParticleExist(slice+skip,ptcl1)*
            Path.ParticleExist(slice+skip,ptcl2);
        } else
          dU += pa.dU(q,z,s2,level);
        // Subtract off long-range part from short-range action
        if (pa.IsLongRange() && PathData.Actions.UseLongRange)
          dU -= 0.5*(pa.dUlong(level)(rmag)+pa.dUlong(level)(rpmag));
      }
    }
  }

  return dU;
}


//Gradient of the action only works in 3d!!!!!!!
void ShortRangeClass::GradAction(int slice1, int slice2,  const Array<int,1> &ptcls, int level, Array<dVec,1> &gradVec)
{

  PathClass &Path = PathData.Path;
  int skip = (1<<level);
  assert (gradVec.size() == ptcls.size());
  for (int pi=0; pi<ptcls.size(); pi++) {
    int ptcl1 = ptcls(pi);
    int species1 = Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0; ptcl2<Path.NumParticles(); ptcl2++) {
      int species2 = Path.ParticleSpeciesNum(ptcl2);
      PairActionFitClass &PA=*(PathData.Actions.PairMatrix(species1,species2));
      if (ptcl1 != ptcl2) {
	for (int slice=slice1; slice<slice2; slice += skip) {
	  //	  double check=0.0;
	  dVec r, rp;
	  double rmag, rpmag, du_dq, du_dz;
	  Path.DistDisp(slice, slice+skip, ptcl1, ptcl2, rmag, rpmag, r, rp);
	  //	  check += dUdR(slice,ptcl1,ptcl2,level)+dUdR(slice+skip,ptcl1,ptcl2,level);
	  double q = 0.5*(rmag+rpmag);
	  double z = (rmag-rpmag);
	  double s2 = dot (r-rp, r-rp);
	  PA.Derivs(q,z,s2,level,du_dq, du_dz);
	  dVec rhat  = (1.0/rmag)*r;
	  dVec rphat = (1.0/rpmag)*rp;
	  
	  double g1 = 1.0;
	  double g2 = 1.0;
	  if (UseLowVariance) {
	    g1 = g(rmag);
	    g2 = g(rpmag);
	  }
	  gradVec(pi) -= (g1*(0.5*du_dq + du_dz)*rhat + 
			  g2*(0.5*du_dq - du_dz)*rphat);
	  //	  cerr<<"VALS: "<<(g1*(0.5*du_dq + du_dz)*rhat + 
	  //			   g2*(0.5*du_dq - du_dz)*rphat)<<" "<<check<<endl;
	  // gradVec(pi) -= (0.5*du_dq*(rhat+rphat) + du_dz*(rhat-rphat));
	  /// Now, subtract off long-range part that shouldn't be in
	  /// here 
	  if (PA.IsLongRange() && PathData.Actions.UseLongRange)
	    gradVec(pi) += 0.5*(PA.Ulong(level).Deriv(rmag)*g1*rhat+
				PA.Ulong(level).Deriv(rpmag)*g2*rphat);
	}
      }
    }
  }

  //  cerr<<"Gradient not working in 2d"<<endl;
  //  assert(1==2);

}


string ShortRangeClass::GetName()
{
  return "ShortRange";
}
