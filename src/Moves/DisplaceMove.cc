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

#include "DisplaceMove.h"

void DisplaceMoveClass::WriteRatio()
{
  NumAttemptedVar.Write(NumAttempted);
  NumAcceptedVar.Write(NumAccepted);
  MultiStageClass::WriteRatio();

  double AcceptRatio = (double)NumAccepted/(double)NumAttempted;
  Sigma *= 1.0 - DesiredAcceptRatio + AcceptRatio; // Recalculate step size
  dVec Box = PathData.Path.GetBox();
  if (Sigma > Box(0))
    Sigma = Box(0)/2.;
  cout << "Sigma: " << Sigma << " AcceptRatio: " << AcceptRatio << endl;
}

void DisplaceStageClass::Accept()
{
  CommonStageClass::Accept();
  Path.RefPath.AcceptCopy();
  //do nothing for now
}

void DisplaceStageClass::Reject()
{
  CommonStageClass::Reject();
  Path.RefPath.RejectCopy();
  //do nothing for now
}

double DisplaceStageClass::Sample (int &slice1, int &slice2, Array<int,1> &activeParticles)
{
  /// Now, choose a random displacement
  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    dVec disp;
    ///    PathData.Path.Random.CommonGaussianVec (Sigma, disp);
#if NDIM==3
    disp(0)=PathData.Path.Random.Common()-0.5;
    disp(1)=PathData.Path.Random.Common()-0.5;
    disp(2)=PathData.Path.Random.Common()-0.5;
#endif
#if NDIM==2
    disp(0)=PathData.Path.Random.Common()-0.5;
    disp(1)=PathData.Path.Random.Common()-0.5;
#endif
    disp=disp*Sigma;

    // Actually displace the path
    SetMode(NEWMODE);
    for (int slice=0; slice<PathData.Path.NumTimeSlices(); slice++)
      PathData.Path(slice, ptcl) = PathData.Path(slice, ptcl) + disp;
  }

  PathClass &Path=PathData.Path;
  // Broadcast the new reference path to all the other processors
  PathData.Path.BroadcastRefPath();

  // And return log sample probability ratio
  return log(1.0);
}

void DisplaceMoveClass::Read (IOSectionClass &in)
{
  string moveName = "Displace";
  CurrentPtcl=0;
  assert(in.ReadVar ("Sigma", Sigma));
  DisplaceStage.Sigma = Sigma;
  Array<string,1> activeSpeciesNames;

  int numToMove;
  assert(in.ReadVar("NumToMove", numToMove));
  SetNumParticlesToMove(numToMove);

  assert(in.ReadVar("DesiredAcceptRatio",DesiredAcceptRatio));

  // Read in the active species.
  assert(in.ReadVar ("ActiveSpecies", activeSpeciesNames));
  Array<int,1> activeSpecies(activeSpeciesNames.size());
  for (int i=0; i<activeSpecies.size(); i++)
    activeSpecies(i) = PathData.Path.SpeciesNum(activeSpeciesNames(i));
  SetActiveSpecies (activeSpecies);

  // // Move all particles at the same time.
  // int totalNum = 0;
  // for (int si=0; si<activeSpecies.size(); si++)
  //   totalNum += PathData.Path.Species(activeSpecies(si)).NumParticles;
  // SetNumParticlesToMove(totalNum);

  // Construct action list
  Array<string,1> samplingActions;
  assert(in.ReadVar("SamplingActions",samplingActions));
  for (int i=0;i<samplingActions.size();i++) {
    cerr<<PathData.Path.Communicator.MyProc()<<" "<<moveName<<" Adding "<<(*PathData.Actions.GetAction(samplingActions(i))).GetName()<<" Action"<<endl;
    DisplaceStage.Actions.push_back(PathData.Actions.GetAction(samplingActions(i)));
  }
  // if (PathData.Path.OrderN)
  //   DisplaceStage.Actions.push_back(&PathData.Actions.ShortRangeOn);
  // else
  //   //  DisplaceStage.Actions.push_back(&PathData.Actions.DiagonalAction);
  //   DisplaceStage.Actions.push_back(&PathData.Actions.ShortRange);
  if (PathData.Path.LongRange) {
    if (PathData.Actions.UseRPA) {
      cerr<<PathData.Path.Communicator.MyProc()<<" "<<moveName<<" Adding LongRangeRPA Action"<<endl;
      DisplaceStage.Actions.push_back(&PathData.Actions.LongRangeRPA);
    } else if (PathData.Path.DavidLongRange) {
      cerr<<PathData.Path.Communicator.MyProc()<<" "<<moveName<<" Adding DavidLongRange Action"<<endl;
      DisplaceStage.Actions.push_back(&PathData.Actions.DavidLongRange);
    } else {
      cerr<<PathData.Path.Communicator.MyProc()<<" "<<moveName<<" Adding LongRange Action"<<endl;
      DisplaceStage.Actions.push_back(&PathData.Actions.LongRange);
    }
  }
  for (int i=0; i<activeSpecies.size(); i++) {
    int speciesNum = activeSpecies(i);
    if ((PathData.Actions.NodalActions(speciesNum)!=NULL)) {
      cerr<<PathData.Path.Communicator.MyProc()<<" "<<moveName<<" "<<activeSpeciesNames(i)<<" Adding Node Action"<<endl;
      DisplaceStage.Actions.push_back(PathData.Actions.NodalActions(speciesNum));
    }
  }
  theSpecies = activeSpecies(0);
  // Now construct stage list
  Stages.push_back(&DisplaceStage);

  ActiveParticles.resize(NumParticlesToMove);
  NumAttempted = 0;
}

void DisplaceMoveClass::MakeMove()
{
  // Next, set timeslices
  Slice1 = 0;
  Slice2 = PathData.Path.NumTimeSlices()-1;

  for (int i=0;i<PathData.Path.NumParticles();i++){
    if (theSpecies == PathData.Path.ParticleSpeciesNum(i)) {
      ActiveParticles(0)=i;
      // Now call MultiStageClass' MakeMove
      MultiStageClass::MakeMove();
      NumAttempted++;
    }
  }

}
