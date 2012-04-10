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

#include "VariationalDisplaceMove.h"

double VariationalDisplaceStageClass::Sample (int &slice1, int &slice2,
				   Array<int,1> &activeParticles)
{
  /// Now, choose a random displacement 
  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    dVec disp;
    PathData.Path.Random.CommonGaussianVec (Sigma, disp);

    // Actually displace the path
    SetMode(NEWMODE);
    for (int slice=0; slice<PathData.NumTimeSlices(); slice++)
      PathData.Path(slice, ptcl) = PathData.Path(slice, ptcl) + disp;
  }

  // And return sample probability ratio
  return 1.0;
}

void
VariationalDisplaceMoveClass::Read (IOSectionClass &in)
{
  cerr<<"Calling read of variational displace move"<<endl;
  // Construct action list
  assert(in.ReadVar("Sigma",VariationalDisplaceStage.Sigma));
  VariationalDisplaceStage.Actions.push_back(&PathData.Actions.ShortRangeOn);
  VariationalDisplaceStage.Actions.push_back(&PathData.Actions.VariationalPI);

  
  // Now construct stage list
  Stages.push_back(&VariationalDisplaceStage);

  ActiveParticles.resize(1);
}


bool 
VariationalDisplaceStageClass::Attempt(int &slice1, int &slice2, 
				       Array<int,1> &activeParticles,
				       double &prevActionChange)
{
  assert (slice1 == 0);
  assert (slice2 == PathData.NumTimeSlices()-1);

  SetMode (NEWMODE);
  double sampleRatio=Sample(slice1,slice2,activeParticles);
  SetMode(OLDMODE);
  double oldAction= PathData.Actions.ShortRangeOn.SingleAction(slice1,slice2,activeParticles,0);//GlobalStageAction(activeParticles);
  //  double oldAction=0.0;
  SetMode(NEWMODE);
  double newAction = PathData.Actions.ShortRangeOn.SingleAction(slice1,slice2,activeParticles,0); //GlobalStageAction(activeParticles);
  //  double newAction=0.0;


  double changeNodalAction;
  changeNodalAction=log(PathData.Actions.TruncatedInverse.SingleAction(slice1,slice2,activeParticles,0));
  double dummy=log(PathData.Actions.VariationalPI.SingleAction(slice1,slice2,activeParticles,0));  
//  changeNodalAction=0.0;
  //  perr << "oldAction = " << oldAction
  //       << "newAction = " << newAction << endl;
  //  double currActionChange=2*(newAction-oldAction-changeNodalAction);
  double currActionChange=(newAction-oldAction+changeNodalAction);
  double logAcceptProb=log(sampleRatio)-currActionChange+prevActionChange;
  bool toAccept = logAcceptProb>=log(PathData.Path.Random.Common()); /// Accept condition
  if (toAccept)
    NumAccepted++;
  NumAttempted++;
  prevActionChange=currActionChange;
  return toAccept;
}


void
VariationalDisplaceMoveClass::MakeMove ()
{
  
  // First, choose particle to move
  int chosenPtcl =PathData.Path.Random.CommonInt(PathData.Path.NumParticles());
  ActiveParticles(0) = chosenPtcl;
  assert(chosenPtcl<PathData.Path.NumParticles());
  // Next, set timeslices
  Slice1 = 0;
  Slice2 = PathData.Path.NumTimeSlices()-1;
  //  cerr<<"Move begun"<<endl;
  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;

  while (stageIter!=Stages.end() && toAccept){
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);
    stageIter++;
  }
  //  cerr<<"After levi flight"<<endl;
  //  PathData.Path.PrintRealSlices();

  TimesCalled++;
  if (toAccept)
    MultiStageClass::Accept();
  else 
    MultiStageClass::Reject();
  //  cerr<<"Move done"<<endl;
  //MoveClass::MakeMove();
}
