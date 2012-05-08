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


#include "OpenStage.h"
#include "EndStage.h"

void OpenStageClass::Read(IOSectionClass &in)
{
  LocalStageClass::Read(in);
  assert(in.ReadVar("NumLevels",NumLevels));
}

void OpenStageClass::WriteRatio()
{
  //  cerr<<"About to write my ratio"<<endl;
  Array<double,1> acceptRatio(2);
  acceptRatio(0)=(double)AcceptRatio(0)/EndAttempts;
  acceptRatio(1)=(double)AcceptRatio(1)/EndAttempts;
  AcceptRatioVar.Write(acceptRatio);
  AcceptRatioVar.Flush();
  //  cerr<<"done writing my ratio"<<endl;
}


void OpenStageClass::Accept()
{
  cerr<<"I've accepted a worm switch"<<endl;
  EndAttempts++;
}

void OpenStageClass::Reject()
{
  cerr<<"Rejecting"<<endl;
  EndAttempts++;
}

dVec OpenStageClass::RandomBoxLocation()
{
  dVec myLoc;
  dVec boxSize=PathData.Path.GetBox();
  myLoc[0]=PathData.Path.Random.Local()*(2*boxSize[0])-boxSize[0];
  myLoc[1]=PathData.Path.Random.Local()*(2*boxSize[1])-boxSize[1];
#ifdef NDIM==3
  myLoc[2]=PathData.Path.Random.Local()*(2*boxSize[2])-boxSize[2];
#endif
  return myLoc;
}


bool OpenStageClass::Attempt (int &slice1, int &slice2, 
	      Array<int,1> &activeParticles, double &prevActionChange)
{
  
  cerr<<"Attempting worm move stage"<<endl;
  //    double actionChange = -log(Forw->CurrentCycle.P/Forw->Gamma[len-1]);
  //    Array<int,1> currentParticles=Forw->CurrentParticles();  
  //    double pi_ratio = exp(-actionChange+prevActionChange);
  if (activeParticles(0)==-1){
    //    double forwT=Sample(slice1,slice2,activeParticles);
    double Tratio = Sample(slice1,slice2,activeParticles); //forwT/revT;
    double acceptProb = min(1.0, 1/Tratio);
    //    prevActionChange = actionChange;

    double psi = PathData.Path.Random.Local();
    return (acceptProb > psi);
  }
  else 
    return true;
  
}




double OpenStageClass::Sample(int &slice1,int &slice2, 
			      Array<int,1> &activeParticles)
{
  cerr<<"Sampling"<<endl;
  if (activeParticles(0)==-1){
  bool shiftWasNecessary=false;
  if (PathData.Path.NowOpen){
    cerr<<"open"<<endl;
    PathData.Path(PathData.Path.OpenLink,PathData.Path.NumParticles())=
      PathData.Path(PathData.Path.OpenLink,PathData.Path.OpenPtcl);
    EndType Open;
    Open=TAIL;
    if (Open==HEAD){
      slice1=(int)PathData.Path.OpenLink;
      slice2=(1<<NumLevels)+slice1;
      if (PathData.Path.Communicator.NumProcs()>1)
	return 0.0;
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
	return 0.0;
      
      ///Shift the time slices so that from the tail there are
      ///2^numlevels slices available  
      while (slice1<0 || slice2==PathData.Path.NumTimeSlices()-1){
	cerr<<"In while"<<endl;
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
    PathData.Path.NowOpen=false;
    activeParticles(0)=PathData.Path.OpenPtcl;
    cerr<<"returning"<<endl;
    return 1.0;
    ///if you were trying to close and this gets accepted you know want
    ///to switch the openb to someone new.  or maybey ou jsut want to stay closed and have a move whose job it is to open..it's not really clear
  }
  else{
    cerr<<"close"<<endl;
    ///need to verify that the link that you choose is such that you can piock time slices to do the bisection over.  Also need to verify that the tail slice is the correct thing!
    PathData.Path.OpenLink=
      PathData.Path.Random.LocalInt(PathData.Path.NumTimeSlices());
    PathData.Path.OpenPtcl=
      PathData.Path.Random.LocalInt(PathData.Path.NumParticles());
    //    dVec newBoxLocation
    //    PathData.Path(PathData.Path.OpenLink,PathData.Path.OpenPtcl)=
    //      RandomBoxLocation();

    PathData.Path(PathData.Path.OpenLink,PathData.Path.NumParticles())=
      RandomBoxLocation();
    

    EndType Open;
    Open=TAIL;
    if (Open==HEAD){
      slice1=(int)PathData.Path.OpenLink;
      slice2=(1<<NumLevels)+slice1;
      if (PathData.Path.Communicator.NumProcs()>1)
	return 0.0;
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
      cerr<<"In while"<<endl;
      slice2=(int)PathData.Path.OpenLink;
      slice1=slice2-(1<<NumLevels);
      if (PathData.Path.Communicator.NumProcs()>1)
	return 0.0;
      
      ///Shift the time slices so that from the tail there are
      ///2^numlevels slices available  
      while (slice1<0 || slice2==PathData.Path.NumTimeSlices()-1){
	cerr<<"Another while"<<endl;
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
    PathData.Path.NowOpen=true;
    activeParticles(0)=PathData.Path.OpenPtcl;
    cerr<<"out"<<endl;

    return 1.0/((double)PathData.Path.NumParticles()*
		PathData.Path.NumTimeSlices());

  }
  cerr<<"You shouldn't actually be here"<<endl;
  return 1.0;
  }
  else{
    //    if (PathData.Path.NowClosed){
    return 1.0;
  }

}

