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
#include "KineticClass.h"

///This has to be called after pathdata knows how many
///particles it has
void KineticClass::Read(IOSectionClass& in)
{
  NumImages = 0;
  in.ReadVar("NumImages",NumImages);
}

KineticClass::KineticClass(PathDataClass &pathData ) : 
  ActionBaseClass (pathData)
{
  TimeSpent = 0.0;
}

double 
KineticClass::SingleAction (int slice1, int slice2,
			    const Array<int,1> &changedParticles, int level)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  double TotalK = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl = changedParticles(ptclIndex);
    int species=Path.ParticleSpeciesNum(ptcl);
    double lambda = Path.Species(species).lambda;
    if (lambda != 0){
      double FourLambdaTauInv=1.0/(4.0*Path.Species(species).lambda*levelTau);
      for (int slice=slice1; slice < slice2;slice+=skip) {
        dVec vel;
	vel = PathData.Path.Velocity(slice, slice+skip, ptcl);
	//vel= PathData.Path(slice+skip,ptcl)-PathData.Path(slice,ptcl);
	
        double GaussProd = 1.0;
        for (int dim=0; dim<NDIM; dim++) {
	  double GaussSum=0.0;
	  for (int image=-NumImages; image<=NumImages; image++) {
	    double dist = vel[dim]+(double)image*Path.GetBox()[dim];
	    GaussSum += exp(-dist*dist*FourLambdaTauInv);
	  }
	  GaussProd *= GaussSum;
        }
	//	if (changedParticles.size()==1)
	//	if (ptcl<3)
	//	  TotalK -= PathData.Path.ExistsCoupling*log(GaussProd);    
	//	else 
	//	if (changedParticles.size()>1)
	//	  TotalK -= 0;
	//	else
	////////////////////////////////	  TotalK -= log(GaussProd);   
	if (PathData.Path.WormOn){
	  if (PathData.Path.ParticleExist(slice,ptcl)*PathData.Path.ParticleExist(slice+skip,ptcl)!=0.0)
	    TotalK -= log(GaussProd);
	}
	else {
	  TotalK -= log(GaussProd);
	}
	//TotalK += dot(vel,vel)*FourLambdaTauInv; 
	//	KineticVal(slice)+=dot(vel,vel)*FourLambdaTauInv; 
      }
    }
  }
  //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
  //  for (int counter=0;counter<KineticVal.size();counter++){
  //    cerr<<"My kinetic link "<<counter<<" is "<<KineticVal(counter)<<endl;
  //  }
  //  if (PathData.Path.ExistsCoupling>0)
  //    TotalK-=100*PathData.Path.ExistsCoupling*log(0.9);
    //    TotalK+=pow(0.8,100*PathData.Path.ExistsCoupling);
  //  cerr<<"My total K is "<<TotalK<<endl;
  //
  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);

  return (TotalK);
}

double 
KineticClass::SingleActionForcedTau (int slice1, int slice2,
				     const Array<int,1> &changedParticles, 
				     int level,
				     double forcedTau)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  double TotalK = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = forcedTau* (1<<level);
  for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
    int species=Path.ParticleSpeciesNum(ptcl);
    double lambda = Path.Species(species).lambda;
    if (lambda != 0){
      double FourLambdaTauInv=1.0/(4.0*Path.Species(species).lambda*levelTau);
      for (int slice=slice1; slice < slice2;slice+=skip) {
        dVec vel;
	vel = PathData.Path.Velocity(slice, slice+skip, ptcl);
	
        double GaussProd = 1.0;
        for (int dim=0; dim<NDIM; dim++) {
	  double GaussSum=0.0;
	  for (int image=-NumImages; image<=NumImages; image++) {
	    double dist = vel[dim]+(double)image*Path.GetBox()[dim];
	    GaussSum += exp(-dist*dist*FourLambdaTauInv);
	  }
	  GaussProd *= GaussSum;
        }
	TotalK -= log(GaussProd);
      }
    }
  }

  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);

  return (TotalK);
}



double KineticClass::d_dBetaForcedTau (int slice1, int slice2,
				       int level,
				       double forcedTau)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  double spring=0.0;
  double levelTau=ldexp(forcedTau, level);
  spring  = 0.0;  
  int skip = 1<<level;
  for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++) {
    // Do free-particle part
    int speciesNum  = Path.ParticleSpeciesNum(ptcl);
    SpeciesClass &species = Path.Species(speciesNum);
    double lambda = species.lambda;
    if (lambda != 0.0) {
      double FourLambdaTauInv = 1.0/(4.0*lambda*levelTau);
      for (int slice=slice1; slice<slice2; slice+=skip) {
	spring += (0.5*NDIM)/levelTau; //*Path.ParticleExist(slice,ptcl);
	dVec vel;
	vel = PathData.Path.Velocity(slice, slice+skip, ptcl);
	double Z = 1.0;
	dVec GaussSum=0.0;
	dVec numSum=0.0;
	for (int dim=0; dim<NDIM; dim++) {
	  for (int image=-NumImages; image<=NumImages; image++) {
	    double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	    double d2overFLT = dist*dist*FourLambdaTauInv;
	    double expPart = exp(-d2overFLT);
	    GaussSum[dim] += expPart;
	    numSum[dim] += -d2overFLT/levelTau* expPart;
	  }
	  Z *= GaussSum[dim];
	}
	double scalarnumSum=0.0;
	for (int dim=0;dim<NDIM;dim++) {
	  dVec numProd=1.0;
	  for (int dim2=0;dim2<NDIM;dim2++) {
	    if (dim2!=dim)
	      numProd[dim] *= GaussSum[dim2];
	    else 
	      numProd[dim] *=  numSum[dim2];
	  }
	  scalarnumSum += numProd[dim];
	} 
	spring += scalarnumSum/Z;
      }
    }
  }
  if (isnan(spring)){
    cerr<<"NAN!"<<endl;
    PathData.Path.PrintRealSlices();
  }

  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);

  return spring;
}





double KineticClass::d_dBeta (int slice1, int slice2,
			      int level)
{
  struct timeval start, end;
  struct timezone tz;
  gettimeofday(&start, &tz);

  double spring=0.0;
  // ldexp(double x, int n) = x*2^n
  double levelTau=ldexp(Path.tau, level);
//   for (int i=0; i<level; i++) 
//     levelTau *= 2.0;
  spring  = 0.0;  
  int skip = 1<<level;
  for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++) {
    // Do free-particle part
    int speciesNum  = Path.ParticleSpeciesNum(ptcl);
    SpeciesClass &species = Path.Species(speciesNum);
    double lambda = species.lambda;
    if (lambda != 0.0) {
      double FourLambdaTauInv = 1.0/(4.0*lambda*levelTau);
      for (int slice=slice1; slice<slice2; slice+=skip) {
	spring += (0.5*NDIM)/levelTau; //*Path.ParticleExist(slice,ptcl);
	dVec vel;
	vel = PathData.Path.Velocity(slice, slice+skip, ptcl);
	double Z = 1.0;
	dVec GaussSum=0.0;
	dVec numSum=0.0;
	for (int dim=0; dim<NDIM; dim++) {
	  for (int image=-NumImages; image<=NumImages; image++) {
	    double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	    double d2overFLT = dist*dist*FourLambdaTauInv;
	    double expPart = exp(-d2overFLT);
	    GaussSum[dim] += expPart;
	    numSum[dim] += -d2overFLT/levelTau* expPart;
	    //	    cerr<<"dist is "<<dist<<" "<<endl;
	  }
	  ///	  cerr<<"GaussSum[dim] is "<<GaussSum[dim]<<" "<<endl;
	  Z *= GaussSum[dim];
	}
	double scalarnumSum=0.0;
	for (int dim=0;dim<NDIM;dim++) {
	  dVec numProd=1.0;
	  for (int dim2=0;dim2<NDIM;dim2++) {
	    if (dim2!=dim)
	      numProd[dim] *= GaussSum[dim2];
	    else 
	      numProd[dim] *=  numSum[dim2];
	  }
	  scalarnumSum += numProd[dim];
	} //cerr << "Z = " << Z << " scalarnumSum = " << scalarnumSum << endl;
	spring += scalarnumSum/Z;
	//	cerr << "spring = " << spring << endl;
	//	cerr << "Z      = " << Z      << endl;
	//	cerr << "scalarnumsum      = " << scalarnumSum      << endl;
	//	cerr << "ptcl is  " <<ptcl<<endl;
	//	cerr<<  "slice is " <<slice<<endl;
      }
    }

  }
  //cerr << "I'm returning kinetic energy " << TotalK << endl;
//   cerr<<"Returning kinetic energy"<<endl;
//   for (int counter=0;counter<KineticVal.size();counter++){
//     cerr<<"My kinetic link "<<counter<<" is "<<KineticVal(counter)<<endl;
//   }
  if (isnan(spring)){
    cerr<<PathData.Path.CloneStr << " NAN in Kinetic "<<endl;
    //PathData.Path.PrintRealSlices();
  }

  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);

  return spring;
}


string
KineticClass::GetName()
{
  return "Kinetic";
}
