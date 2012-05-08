////////////////////////////////////////////////////////////
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


#include "LeviFlightStage.h"

void LeviFlightStageClass::Read(IOSectionClass &in)

{
//   string speciesName;
//   assert (in.ReadVar ("Species", speciesName));
//   SpeciesNum = PathData.Path.SpeciesNum (speciesName);
  
}
 
void LeviFlightStageClass::WriteRatio()
{ 
  AcceptRatioVar.Write((double)NumAccepted/(double)NumAttempted);
  AcceptRatioVar.Flush();
}

bool LeviFlightStageClass::Attempt(int &slice1, int &slice2, 
			      Array<int,1> &activeParticles,
			      double &prevActionChange)
{
  cerr<<"Attempting"<<endl;
  SetMode (NEWMODE);
  double sampleRatio=Sample(slice1,slice2,activeParticles);
  SetMode(OLDMODE);
  double oldAction=StageAction(slice1,slice2,activeParticles);
  SetMode(NEWMODE);
  double newAction =StageAction(slice1,slice2,activeParticles);
  double currActionChange=newAction-oldAction;
  double logAcceptProb=log(sampleRatio)-currActionChange+prevActionChange;
  cerr<<log(sampleRatio)<<" "<<currActionChange<<" "<<prevActionChange
      <<" "<<log(sampleRatio)-currActionChange+prevActionChange<<endl;
  bool toAccept = logAcceptProb>=log(PathData.Path.Random.Local()); /// Accept condition
  if (toAccept)
    NumAccepted++;
  NumAttempted++;
  prevActionChange=currActionChange;
  cerr<<"Done attempting"<<endl;
  return toAccept;
}


void LeviFlightStageClass::Accept()
{
  //do nothing for now
  
}

void LeviFlightStageClass::Reject()
{
  //do nothing for now

}

double 
LeviFlightStageClass::Sample(int &slice1,int &slice2,
			     Array<int,1> &activeParticles)
{
  PathClass &Path = PathData.Path;
  double logSampleProb=0.0;
  assert(activeParticles.size()==1);
  //you are about to close and this needs to be lcosed for the action
  assert(!PathData.Path.NowOpen);
  int ptcl=activeParticles(0);
  double lambda=Path.Species(Path.ParticleSpeciesNum(ptcl)).lambda;
  logSampleProb+=LeviFlight(slice1,slice2,ptcl,lambda);
  return logSampleProb;
}



  
/// Constructs a Levi flight beginning in the vec(0) and ending
/// in vec(N-1) if vec has length N.  This is a path which samples
/// exactly the free particle density matrix.
double 
LeviFlightStageClass::LeviFlight(int slice1,int slice2, int ptcl, double lambda)
{
  PathClass &Path = PathData.Path;
  double logSampleProb=0.0;
  int N = slice2-slice1+1;
  dVec dispToBuild=Path(slice2,ptcl)-Path(slice1,ptcl);
  PathData.Path.PutInBox(dispToBuild);
  Path(slice2,ptcl)=Path(slice1,ptcl)+dispToBuild;
  for (int slice=slice1+1; slice<slice2; slice++) {
    double di = (double)slice;
    double delta = (double)(N-(slice-slice1)-1);
    dVec center = (1.0/(delta+1.0))*(delta*Path(slice-1,ptcl) + Path(slice2,ptcl));

    double taueff = PathData.Path.tau*(1.0 - 1.0/(delta+1.0));
    double sigma2= (2.0*lambda*taueff);
    double sigma = sqrt (sigma2);
    dVec disp;
    Path.Random.LocalGaussianVec(sigma, disp);
    double dist2=dot(disp,disp);
    logSampleProb+=-0.5*dist2/sigma2;
    Path(slice,ptcl) = center+disp;
  }
  int  nSlices=slice2-slice1;
  double dist2=dot(dispToBuild,dispToBuild);
  //  return exp(-dist2/(4.0*lambda*nSlices*PathData.Path.tau));; //
  return exp(-logSampleProb);
}



