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


#include "EndStage.h"
void EndStageClass::Read(IOSectionClass &in)
{
  OnlyX=false;
  in.ReadVar("OnlyX",OnlyX);
  LocalStageClass::Read(in);
}

void EndStageClass::WriteRatio()
{
  //  cerr<<"About to write my ratio"<<endl;
  Array<double,1> acceptRatio(2);
  acceptRatio(0)=(double)AcceptRatio(0)/EndAttempts;
  acceptRatio(1)=(double)AcceptRatio(1)/EndAttempts;
  AcceptRatioVar.Write(acceptRatio);
  AcceptRatioVar.Flush();
  //  cerr<<"done writing my ratio"<<endl;
}


void EndStageClass::Accept()
{
  if (Open==HEAD)
    AcceptRatio(0)++;
  else 
    AcceptRatio(1)++;
  EndAttempts++;

}

void EndStageClass::Reject()
{
  EndAttempts++;

}


///Chooses the time slices and moves the join so that the join is in
///the correct place for that time slice.
void EndStageClass::ChooseTimeSlices(int &slice1,int &slice2)
{
  bool shiftWasNecessary=false;
  if (PathData.Path.Random.Local()>0.5)
    Open=HEAD;
  else
    Open=TAIL;
    
  //  Open=HEAD;
  //HACK!
  //  Open=TAIL;
  if (Open==HEAD){
    slice1=(int)PathData.Path.OpenLink;
    slice2=(1<<NumLevels)+slice1;
    
    if (PathData.Path.Communicator.NumProcs()>1)
      return;
    ///Shift the time slices so that from the head there are
    ///2^numlevels slices available  
    while (slice2>=PathData.Path.NumTimeSlices() || slice1==0){
      shiftWasNecessary=true;
      PathData.MoveJoin(0);
      PathData.ShiftData(2);
      PathData.Join=2;
      slice1=(int)PathData.Path.OpenLink;
      slice2=(1<<NumLevels)+slice1;
    }
  }
  else if (Open==TAIL){
    slice2=(int)PathData.Path.OpenLink;
    slice1=slice2-(1<<NumLevels);
    if (PathData.Path.Communicator.NumProcs()>1)
      return;

    ///Shift the time slices so that from the tail there are
    ///2^numlevels slices available  
    while (slice1<0 || slice2==PathData.Path.NumTimeSlices()-1){
      shiftWasNecessary=true;
      PathData.MoveJoin(0);
      PathData.ShiftData(2);
      PathData.Join=2;
      slice2=(int)PathData.Path.OpenLink;
      slice1=slice2-(1<<NumLevels);
    }
  }
  else {
    cerr<<"ERROR! We don't know whether to do head or tail?"<<endl;
    assert(1==2);
  }
  //  PathData.MoveJoin(slice2);
  if (shiftWasNecessary && PathData.Path.OrderN){
    for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++)
      PathData.Path.Cell.BinParticles(slice);
  }

    
}


///HACK! HACK! HACK!
double mysign(double num)
{
  if (num<0)
    return -1.0;
  else return 1.0;
}

double EndStageClass::Sample(int &slice1,int &slice2, 
	      Array<int,1> &activeParticles)
{
  //  cerr<<"I am sampling now"<<endl;
  int procWithRefSlice = PathData.Path.SliceOwner (PathData.Path.RefSlice);
  //  cerr<<"Entering end stage class" <<procWithRefSlice<<" "
  //      <<PathData.Path.RefSlice<<" "
  //      <<PathData.Path.Communicator.MyProc()<<endl;
  
  ///If you don't own the open link  
  ///then just choose some slices to work with
  if (procWithRefSlice != PathData.Path.Communicator.MyProc()) {
    int sliceSep = 1<<NumLevels;
    assert (sliceSep < PathData.Path.NumTimeSlices());
    int numLeft = PathData.Path.NumTimeSlices()-sliceSep;
    slice1 = PathData.Path.Random.LocalInt (numLeft);
    slice2 = slice1+sliceSep;
    ///BUG: Will only work for one species
    activeParticles.resize(1);
    activeParticles(0)=
      PathData.Path.Random.LocalInt(PathData.Path.NumParticles());
    PathData.MoveJoin(slice2);
    //    cerr<<"Exiting as an incorrect processor"<<endl;
    return 1.0;
  }
  ChooseTimeSlices(slice1,slice2);
    
  ///if you are running with time slice parallelization it doesn't
  ///shift the data since everyone doesn't have to know about it. It
  ///simply returns 1 allowing everyone to do their own bisection
  if ((slice1<=0 || slice2>=PathData.Path.NumTimeSlices()-1) &&
      PathData.Path.Communicator.NumProcs()>1){
        int sliceSep = 1<<NumLevels;
    assert (sliceSep < PathData.Path.NumTimeSlices());
    int numLeft = PathData.Path.NumTimeSlices()-sliceSep;
    slice1 = PathData.Path.Random.LocalInt (numLeft);
    slice2 = slice1+sliceSep;
    ///BUG: Will only work for one species
    activeParticles.resize(1);
    activeParticles(0)=
      PathData.Path.Random.LocalInt(PathData.Path.NumParticles());
    PathData.MoveJoin(slice2);
    //    cerr<<"Exiting as a correct processor with bad time slices"<<endl;
    return 1.0;
  }
  PathData.MoveJoin(slice2);
  ///The active particle should be a slice
  activeParticles.resize (1);
  activeParticles(0) = (int)PathData.Path.OpenPtcl;

  int changePtcl;
  int skipSign;
  dVec oldrdiff;
  dVec oldPos;
  if (Open==HEAD){
    // cerr<<"Setting to head"<<endl;
    changePtcl=(int)PathData.Path.OpenPtcl;
    oldrdiff = Path.Velocity(slice1,slice2,changePtcl);
    oldPos = PathData.Path(slice2,changePtcl);
  }
  else if (Open==TAIL){
    // cerr<<"Setting to Tail"<<endl;
    changePtcl=PathData.Path.NumParticles();
    oldrdiff = Path.Velocity(slice1, slice2, (int)PathData.Path.OpenPtcl);
    oldPos = PathData.Path(slice1,(int)PathData.Path.OpenPtcl);
  }
  else {
    cerr<<
      "I'm neither seeing a head or tail here! Warning! Warning!"<<endl;
    abort();
  }

  int skip = (1<<NumLevels);
  double levelTau = PathData.Path.tau*skip;
  double lambda=PathData.Path.ParticleSpecies(activeParticles(0)).lambda;
  double sigma2=(2.0*lambda*levelTau);
  double sigma=sqrt(sigma2);

  dVec Delta;
  Path.Random.LocalGaussianVec(sigma,Delta);
  PathData.Path.PutInBox(Delta);

  dVec newPos= oldPos + Delta;

  dVec newrdiff = Delta;

  double logOldSampleProb = -0.5*dot(oldrdiff,oldrdiff)/sigma2;
  double logSampleProb = -0.5*dot(newrdiff,newrdiff)/sigma2;

  PathData.Path.SetPos(PathData.Path.OpenLink,changePtcl,newPos);

  return exp(-logSampleProb+logOldSampleProb);

  // Throw in Box
  //
  //  cerr<<"Existing as a correct processor with correct time slices"<<endl;
  //return 1.0;
  //newPos(0)= oldPos(0)+
  //  //    PathData.Path.Random.Local()*(PathData.Path.GetBox()(0)/10.0)*
  //  //HACK
  //  PathData.Path.Random.Local()*(PathData.Path.GetBox()(0)/5.0)*
  //  mysign(PathData.Path.Random.Local()-0.5);
  //newPos(1)= oldPos(1)+
  //  PathData.Path.Random.Local()*(PathData.Path.GetBox()(1)/5.0)*
  //  mysign(PathData.Path.Random.Local()-0.5);
  ////newPos(2)=  oldPos(2)+
  ////  PathData.Path.Random.Local()*(PathData.Path.GetBox()(2)/5.0)*
  ////  mysign(PathData.Path.Random.Local()-0.5);
  //if (OnlyX){
  //  newPos(1)=oldPos(0);
  //  //newPos(2)=oldPos(2);
  //}
  //if (fabs(newPos(0))>500 || fabs(newPos(1))>500) { // || fabs(newPos(2))>500){
  //  cerr<<"ERROR! ERROR! ERROR!"<<newPos(0)
	//<<" "<<newPos(1)<<endl;//" "<<newPos(2)<<endl;
  //}
}
