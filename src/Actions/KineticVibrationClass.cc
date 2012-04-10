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
#include "KineticVibrationClass.h"
#include "../Moves/MoveUtils.h"
#include <Common/SpecialFunctions/HermitePoly.h>

///This has to be called after pathdata knows how many
///particles it has
void KineticVibrationClass::Read(IOSectionClass& in)
{
  NumImages = 0;
  in.ReadVar("NumImages",NumImages);
  assert(in.ReadVar("Lambda",lambda));
  cerr << "KineticVibration Action Class: right now this is specialized for water molecules; generalization may not be reliable without modifications" << endl;
  //Array<string,1> speciesList;
  //Array<string,1> LambdaList;
  //assert(in.ReadVar("SpeciesList",speciesList));
  //assert(in.ReadVar("LambdaList",LambdaList));
  //assert(LambdaList.size() == SpeciesList.size());
  //cerr << "KineticVibration read; active species are " << speciesList << endl;
  //DoSpecies.resize(PathData.Path.NumSpecies());
  //Lambdas.resize(PathData.Path.NumSpecies());
  //Lambdas = 0.0;
  //DoSpecies = false;
  //for(int s=0; s<speciesList.size(); s++) {
  //  int sNum = PathData.Path.SpeciesNum(speciesList(s));
  //  DoSpecies(sNum) = true;
  //  Lambdas(sNum) = atof(LambdaList(s).c_str());

  //}
  //cerr << "DoSpecies is " << DoSpecies << endl;
  //cerr << "Lambdas is " << Lambdas << endl;
}

KineticVibrationClass::KineticVibrationClass(PathDataClass &pathData ) : 
  ActionBaseClass (pathData)
{
}

double 
KineticVibrationClass::SingleAction (int slice1, int slice2,
			    const Array<int,1> &changedParticles, int level)
{
  double TotalK = 0.0;
  // need to map activePtcls to active molecules - probably could be done more cleanly
  Array<bool,1> activeMol;
  activeMol.resize(PathData.Mol.NumMol());
  activeMol = false;
  for(int pindex=0; pindex<changedParticles.size(); pindex++) {
    int myMol = PathData.Mol(changedParticles(pindex));
    activeMol(myMol) = true;
  }
  vector<int> molRoster(0);
  for(int mindex=0; mindex<activeMol.size(); mindex++) {
    if(activeMol(mindex))
      molRoster.push_back(mindex);
  }

  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  double FourLambdaTauInv=1.0/(4.0*lambda*levelTau);
  int numChangedMol = molRoster.size();
  for (int molIndex=0; molIndex<numChangedMol; molIndex++){
    int mol = molRoster[molIndex];
    for (int slice=slice1; slice < slice2;slice+=skip) {
      Array<int,1> molMembers;
      PathData.Mol.MembersOf(molMembers, mol);
      dVec O1 = PathData.Path(slice, molMembers(0));
      dVec O2 = PathData.Path(slice+skip, molMembers(0));
      dVec p1a = PathData.Path(slice, molMembers(1)) - O1;
      dVec p1b = PathData.Path(slice, molMembers(2)) - O1;
      dVec p2a = PathData.Path(slice+skip, molMembers(1)) - O2;
      dVec p2b = PathData.Path(slice+skip, molMembers(2)) - O2;
      // must impose minimum image for molecules at boundaries
      PathData.Path.PutInBox(p1a);
      PathData.Path.PutInBox(p1b);
      PathData.Path.PutInBox(p2a);
      PathData.Path.PutInBox(p2b);
      double p1aMag = Mag(p1a);
      double p1bMag = Mag(p1b);
      double p2aMag = Mag(p2a);
      double p2bMag = Mag(p2b);
      double dra = p1aMag - p2aMag;
      double drb = p1bMag - p2bMag;
      double rAvg = 0.25*(p1aMag + p1bMag + p2aMag + p2bMag);
      p1a = Normalize(p1a);
      p1b = Normalize(p1b);
      p2a = Normalize(p2a);
      p2b = Normalize(p2b);
      double phi1 = GetAngle(p1a,p1b);
      double phi2 = GetAngle(p2a,p2b);
      double dArc = rAvg * (phi2 - phi1);
      double distSq = dArc*dArc + dra*dra + drb*drb;

      // now get the effective vel right!
      //dVec vel;
      //vel = PathData.Path.Velocity(slice, slice+skip, ptcl);
      double GaussProd = exp(-distSq * FourLambdaTauInv);
      //double GaussProd = 1.0;
      //for (int dim=0; dim<NDIM; dim++) {
      //  double GaussSum= 0.;
      //  //for (int image=-NumImages; image<=NumImages; image++) {
      //  //  double dist = vel[dim]+(double)image*Path.GetBox()[dim];
      //  //  GaussSum += exp(-dist*dist*FourLambdaTauInv);
      //  //}
      //  GaussProd *= GaussSum;
      //}
      TotalK -= log(GaussProd);
    }
  }
  return (TotalK);
}

double 
KineticVibrationClass::SingleActionForcedTau (int slice1, int slice2,
				     const Array<int,1> &changedParticles, 
				     int level,
				     double forcedTau)
{
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
  return (TotalK);
}



double KineticVibrationClass::d_dBetaForcedTau (int slice1, int slice2,
				       int level,
				       double forcedTau)
{
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
  return spring;
}





double KineticVibrationClass::d_dBeta (int slice1, int slice2,
			      int level)
{
  //// need to map activePtcls to active molecules - probably could be done more cleanly
  //Array<bool,1> activeMol;
  //activeMol.resize(PathData.Mol.NumMol());
  //activeMol = false;
  //for(int pindex=0; pindex<changedParticles.size(); pindex++) {
  //  int myMol = PathData.Mol(changedParticles(pindex));
  //  activeMol(myMol) = true;
  //}
  //vector<int> molRoster(0);
  //for(int mindex=0; mindex<activeMol.size(); mindex++) {
  //  if(activeMol(mindex))
  //    molRoster.push_back(mindex);
  //}

  double spring = 0.0;
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  double FourLambdaTauInv=1.0/(4.0*lambda*levelTau);
  //int numChangedMol = molRoster.size();
  for (int molIndex=0; molIndex<PathData.Mol.NumMol(); molIndex++){
    int mol = molIndex;
    for (int slice=slice1; slice < slice2;slice+=skip) {
      Array<int,1> molMembers;
      PathData.Mol.MembersOf(molMembers, mol);
      dVec O1 = PathData.Path(slice, molMembers(0));
      dVec O2 = PathData.Path(slice+skip, molMembers(0));
      dVec p1a = PathData.Path(slice, molMembers(1)) - O1;
      dVec p1b = PathData.Path(slice, molMembers(2)) - O1;
      dVec p2a = PathData.Path(slice+skip, molMembers(1)) - O2;
      dVec p2b = PathData.Path(slice+skip, molMembers(2)) - O2;
      //cerr << O1 << " " << p1a << " " << p1b << endl;
      //cerr << O2 << " " << p2a << " " << p2b << endl;
      // must impose minimum image for molecules at boundaries
      PathData.Path.PutInBox(p1a);
      PathData.Path.PutInBox(p1b);
      PathData.Path.PutInBox(p2a);
      PathData.Path.PutInBox(p2b);
      //cerr << p1a << " " << p1b << endl;
      //cerr << p2a << " " << p2b << endl;
      double p1aMag = Mag(p1a);
      double p1bMag = Mag(p1b);
      double p2aMag = Mag(p2a);
      double p2bMag = Mag(p2b);
      double dra = p1aMag - p2aMag;
      double drb = p1bMag - p2bMag;
      //double rAvg = 0.25*(Mag(p1a) + Mag(p1b) + Mag(p2a) + Mag(p2b));
      double rAvg = 0.25*(p1aMag + p1bMag + p2aMag + p2bMag);
      p1a = Normalize(p1a);
      p1b = Normalize(p1b);
      p2a = Normalize(p2a);
      p2b = Normalize(p2b);
      double phi1 = GetAngle(p1a,p1b);
      double phi2 = GetAngle(p2a,p2b);
      double dArc = rAvg * (phi2 - phi1);
      double distSq = dArc*dArc + dra*dra + drb*drb;
      //cerr << slice << " " << slice+skip << " " << dArc << " " << dra << " " << drb << endl;
      double GaussProd = exp(-distSq * FourLambdaTauInv);
      spring += distSq * FourLambdaTauInv/levelTau;
      // I don't think this is correct yet.
    }
  }
  //cerr << "returning " << spring << endl;
  return (spring);
}




string
KineticVibrationClass::GetName()
{
  return "Kinetic";
}

void KineticVibrationEigenFunctionClass::Read(IOSectionClass& in)
{
  assert(in.ReadVar("Lambda",lambda));
  assert(in.ReadVar("ROH",R_OH_0));
  assert(in.ReadVar("RHH",R_HH_0));
  assert(in.ReadVar("OmegaOH",omega_OH));
  assert(in.ReadVar("OmegaHH",omega_HH));
  assert(in.ReadVar("MaxH",maxH));
  cerr << "KineticVibration EigenFunction Expansion Action Class: right now this is specialized for water molecules; generalization may not be reliable without modifications" << endl;
}

KineticVibrationEigenFunctionClass::KineticVibrationEigenFunctionClass(PathDataClass &pathData ) : 
  ActionBaseClass (pathData)
{
}

double ComputeRho(int maxH, double x, double xprime, double omegaOverLambda, double tauE0) {
  double psi = 0.0;
  Array<double,1> HermiteX(maxH), HermiteXprime(maxH);
  HermitePoly(sqrt(omegaOverLambda)*x, HermiteX);
  HermitePoly(sqrt(omegaOverLambda)*xprime, HermiteXprime);
  double n = 1.0;
  for (int h=0; h<=maxH; h++) {
    if(h>0)
      n *= h;
    psi += sqrt(omegaOverLambda/M_PI) / (pow(2.,h)*n)
      * HermiteX(h) * HermiteXprime(h)
      * exp(-0.5 * omegaOverLambda * x*x)
      * exp(-0.5 * omegaOverLambda * xprime*xprime)
      * exp(-1 * (0.5 + h) * tauE0);
  }
  return psi;
}

double ComputeEigenEnergy(int maxH, double x, double xprime, double omegaOverLambda, double tauE0, double E0, double& rho, double& E) {
  double psi = 0.0;
  rho = 0.0;
  E = 0.0;
  Array<double,1> HermiteX(maxH), HermiteXprime(maxH);
  HermitePoly(sqrt(omegaOverLambda)*x, HermiteX);
  HermitePoly(sqrt(omegaOverLambda)*xprime, HermiteXprime);
  double n = 1.0;
  for (int h=0; h<=maxH; h++) {
    if(h>0)
      n *= h;
    psi = sqrt(omegaOverLambda/M_PI) / (pow(2.,h)*n)
      * HermiteX(h) * HermiteXprime(h)
      * exp(-0.5 * omegaOverLambda * x*x)
      * exp(-0.5 * omegaOverLambda * xprime*xprime)
      * exp(-1 * (0.5 + h) * tauE0);
    rho += psi;
    E += (0.5 + h) * E0 * psi;
  }
}

double 
KineticVibrationEigenFunctionClass::SingleAction (int slice1, int slice2,
			    const Array<int,1> &changedParticles, int level)
{
  double TotalK = 0.0;
  // need to map activePtcls to active molecules - probably could be done more cleanly
  Array<bool,1> activeMol;
  activeMol.resize(PathData.Mol.NumMol());
  activeMol = false;
  for(int pindex=0; pindex<changedParticles.size(); pindex++) {
    int myMol = PathData.Mol(changedParticles(pindex));
    activeMol(myMol) = true;
  }
  vector<int> molRoster(0);
  for(int mindex=0; mindex<activeMol.size(); mindex++) {
    if(activeMol(mindex))
      molRoster.push_back(mindex);
  }

  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  double FourLambdaTauInv=1.0/(4.0*lambda*levelTau);
  int numChangedMol = molRoster.size();
  // compute constants
  double prefactorOH = 0.5 * omega_OH / lambda;
  double prefactorHH = 0.5 * omega_HH / lambda;
  for (int molIndex=0; molIndex<numChangedMol; molIndex++){
    int mol = molRoster[molIndex];
    for (int slice=slice1; slice < slice2;slice+=skip) {
      Array<int,1> molMembers;
      PathData.Mol.MembersOf(molMembers, mol);
      dVec O1 = PathData.Path(slice, molMembers(0));
      dVec O2 = PathData.Path(slice+skip, molMembers(0));
      dVec p1a = PathData.Path(slice, molMembers(1)) - O1;
      dVec p1b = PathData.Path(slice, molMembers(2)) - O1;
      dVec p2a = PathData.Path(slice+skip, molMembers(1)) - O2;
      dVec p2b = PathData.Path(slice+skip, molMembers(2)) - O2;
      // must impose minimum image for molecules at boundaries
      PathData.Path.PutInBox(p1a);
      PathData.Path.PutInBox(p1b);
      PathData.Path.PutInBox(p2a);
      PathData.Path.PutInBox(p2b);
      double p1aMag = Mag(p1a);
      double p1bMag = Mag(p1b);
      double p2aMag = Mag(p2a);
      double p2bMag = Mag(p2b);
      double Rho = 1.0;
      Rho *= ComputeRho(maxH, (p1aMag - R_OH_0), (p2aMag - R_OH_0), prefactorOH, (levelTau*omega_OH)); 
      Rho *= ComputeRho(maxH, (p1bMag - R_OH_0), (p2bMag - R_OH_0), prefactorOH, (levelTau*omega_OH)); 

      double dr1 = Mag(p1a - p1b) - R_HH_0;
      double dr2 = Mag(p2a - p2b) - R_HH_0;
      Rho *= ComputeRho(maxH, dr1, dr2, prefactorHH, (levelTau*omega_HH)); 
      TotalK -= log(Rho);
    }
  }
  return (TotalK);
}

double 
KineticVibrationEigenFunctionClass::SingleActionForcedTau (int slice1, int slice2,
				     const Array<int,1> &changedParticles, 
				     int level,
				     double forcedTau)
{
  double TotalK = 0.0;
  //int numChangedPtcls = changedParticles.size();
  //int skip = 1<<level;
  //double levelTau = forcedTau* (1<<level);
  //for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++){
  //  int species=Path.ParticleSpeciesNum(ptcl);
  //  double lambda = Path.Species(species).lambda;
  //  if (lambda != 0){
  //    double FourLambdaTauInv=1.0/(4.0*Path.Species(species).lambda*levelTau);
  //    for (int slice=slice1; slice < slice2;slice+=skip) {
  //      dVec vel;
	//vel = PathData.Path.Velocity(slice, slice+skip, ptcl);
	//
  //      double GaussProd = 1.0;
  //      for (int dim=0; dim<NDIM; dim++) {
	//  double GaussSum=0.0;
	//  for (int image=-NumImages; image<=NumImages; image++) {
	//    double dist = vel[dim]+(double)image*Path.GetBox()[dim];
	//    GaussSum += exp(-dist*dist*FourLambdaTauInv);
	//  }
	//  GaussProd *= GaussSum;
  //      }
	//TotalK -= log(GaussProd);
  //    }
  //  }
  //}
  return (TotalK);
}

double KineticVibrationEigenFunctionClass::d_dBetaForcedTau (int slice1, int slice2,
				       int level,
				       double forcedTau)
{
  double spring=0.0;
  //double levelTau=ldexp(forcedTau, level);
  //spring  = 0.0;  
  //int skip = 1<<level;
  //for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++) {
  //  // Do free-particle part
  //  int speciesNum  = Path.ParticleSpeciesNum(ptcl);
  //  SpeciesClass &species = Path.Species(speciesNum);
  //  double lambda = species.lambda;
  //  if (lambda != 0.0) {
  //    double FourLambdaTauInv = 1.0/(4.0*lambda*levelTau);
  //    for (int slice=slice1; slice<slice2; slice+=skip) {
	//spring += (0.5*NDIM)/levelTau; //*Path.ParticleExist(slice,ptcl);
	//dVec vel;
	//vel = PathData.Path.Velocity(slice, slice+skip, ptcl);
	//double Z = 1.0;
	//dVec GaussSum=0.0;
	//dVec numSum=0.0;
	//for (int dim=0; dim<NDIM; dim++) {
	//  for (int image=-NumImages; image<=NumImages; image++) {
	//    double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	//    double d2overFLT = dist*dist*FourLambdaTauInv;
	//    double expPart = exp(-d2overFLT);
	//    GaussSum[dim] += expPart;
	//    numSum[dim] += -d2overFLT/levelTau* expPart;
	//  }
	//  Z *= GaussSum[dim];
	//}
	//double scalarnumSum=0.0;
	//for (int dim=0;dim<NDIM;dim++) {
	//  dVec numProd=1.0;
	//  for (int dim2=0;dim2<NDIM;dim2++) {
	//    if (dim2!=dim)
	//      numProd[dim] *= GaussSum[dim2];
	//    else 
	//      numProd[dim] *=  numSum[dim2];
	//  }
	//  scalarnumSum += numProd[dim];
	//} 
	//spring += scalarnumSum/Z;
  //    }
  //  }
  //}
  //if (isnan(spring)){
  //  cerr<<"NAN!"<<endl;
  //  PathData.Path.PrintRealSlices();
  //}
  return spring;
}

double KineticVibrationEigenFunctionClass::d_dBeta (int slice1, int slice2,
			      int level)
{
  double TotalK = 0.0;
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  double FourLambdaTauInv=1.0/(4.0*lambda*levelTau);
  // compute constants
  double prefactorOH = 0.5 * omega_OH / lambda;
  double prefactorHH = 0.5 * omega_HH / lambda;
  for (int molIndex=0; molIndex<PathData.Mol.NumMol(); molIndex++){
    int mol = molIndex;
    for (int slice=slice1; slice < slice2;slice+=skip) {
      Array<int,1> molMembers;
      PathData.Mol.MembersOf(molMembers, mol);
      dVec O1 = PathData.Path(slice, molMembers(0));
      dVec O2 = PathData.Path(slice+skip, molMembers(0));
      dVec p1a = PathData.Path(slice, molMembers(1)) - O1;
      dVec p1b = PathData.Path(slice, molMembers(2)) - O1;
      dVec p2a = PathData.Path(slice+skip, molMembers(1)) - O2;
      dVec p2b = PathData.Path(slice+skip, molMembers(2)) - O2;
      // must impose minimum image for molecules at boundaries
      PathData.Path.PutInBox(p1a);
      PathData.Path.PutInBox(p1b);
      PathData.Path.PutInBox(p2a);
      PathData.Path.PutInBox(p2b);
      double p1aMag = Mag(p1a);
      double p1bMag = Mag(p1b);
      double p2aMag = Mag(p2a);
      double p2bMag = Mag(p2b);
      double Rho = 1.0;
      double r, e;
      ComputeEigenEnergy(maxH, (p1aMag - R_OH_0), (p2aMag - R_OH_0), prefactorOH, (levelTau*omega_OH), omega_OH, r, e); 
      TotalK += e/r;
      Rho *= ComputeEigenEnergy(maxH, (p1bMag - R_OH_0), (p2bMag - R_OH_0), prefactorOH, (levelTau*omega_OH), omega_OH, r, e); 
      TotalK += e/r;
      double dr1 = Mag(p1a - p1b) - R_HH_0;
      double dr2 = Mag(p2a - p2b) - R_HH_0;
      Rho *= ComputeEigenEnergy(maxH, dr1, dr2, prefactorHH, (levelTau*omega_HH), omega_HH, r, e); 
      TotalK += e/r;
    }
  }
  return (TotalK);
}

string
KineticVibrationEigenFunctionClass::GetName()
{
  return "Kinetic";
}
