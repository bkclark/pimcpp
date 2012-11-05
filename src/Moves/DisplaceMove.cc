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
  if (PathData.Path.Equilibrate && DesiredAcceptRatio>0)
    Sigma *= 1.0 - DesiredAcceptRatio + AcceptRatio; // Recalculate step size
  dVec Box = PathData.Path.GetBox();
  if (Sigma > Box(0))
    Sigma = Box(0)/2.;
}

void DisplaceStageClass::Accept()
{
  CommonStageClass::Accept();
  Path.RefPath.AcceptCopy();
  if (Path.StoreNodeDist)
    Path.NodeDist.AcceptCopy();
  if (Path.StoreNodeDet)
    Path.NodeDet.AcceptCopy();
}

void DisplaceStageClass::Reject()
{
  CommonStageClass::Reject();
  Path.RefPath.RejectCopy();
  if (Path.StoreNodeDist)
    Path.NodeDist.RejectCopy();
  if (Path.StoreNodeDet)
    Path.NodeDet.RejectCopy();
}

double DisplaceStageClass::Sample (int &slice1, int &slice2, Array<int,1> &activeParticles)
{

  // See if we have any with just identity permutation
  int N = PathData.Path.NumParticles();
  Array<int,1> TotalPerm(N), doDisplace(activeParticles.size());
  PathData.Path.TotalPermutation(TotalPerm);
  /// Only proc 0 gets TotalPerm
  doDisplace = 0;
  if (Path.Communicator.MyProc() == 0) {
    for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
      int ptcl = activeParticles(ptclIndex);
      if (TotalPerm(ptcl) == ptcl) // Only displace identity permutations
        doDisplace(ptclIndex) = 1;
    }
  }
  PathData.Path.Communicator.Broadcast(0, doDisplace);

  /// Now, choose a random displacement
  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    if (doDisplace(ptclIndex)) {
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
  }

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

  // Read in the active species.
  assert(in.ReadVar ("ActiveSpecies", activeSpeciesNames));
  activeSpecies.resize(activeSpeciesNames.size());
  for (int i=0; i<activeSpecies.size(); i++)
    activeSpecies(i) = PathData.Path.SpeciesNum(activeSpeciesNames(i));
  SetActiveSpecies (activeSpecies);

  // Determine number of particles to move
  int numToMove = 0;
  MoveAllParticles = false;
  in.ReadVar("MoveAll",MoveAllParticles);
  if (MoveAllParticles) {
    for (int i=0; i<activeSpecies.size(); i++) {
      int speciesNum = activeSpecies(i);
      numToMove += PathData.Path.Species(speciesNum).NumParticles;
    }
  } else
    assert(in.ReadVar("NumToMove", numToMove));
  SetNumParticlesToMove(numToMove);

  DesiredAcceptRatio=-1;
  in.ReadVar("DesiredAcceptRatio",DesiredAcceptRatio);

  // // Move all particles at the same time.
  // int totalNum = 0;
  // for (int si=0; si<activeSpecies.size(); si++)
  //   totalNum += PathData.Path.Species(activeSpecies(si)).NumParticles;
  // SetNumParticlesToMove(totalNum);

  // Construct action list
  Array<string,1> samplingActions;
  assert(in.ReadVar("SamplingActions",samplingActions));
  for (int i=0;i<samplingActions.size();i++) {
    cout<<PathData.Path.CloneStr<<" "<<moveName<<" Adding "<<(*PathData.Actions.GetAction(samplingActions(i))).GetName()<<" Action"<<endl;
    DisplaceStage.Actions.push_back(PathData.Actions.GetAction(samplingActions(i)));
  }
  // if (PathData.Path.OrderN)
  //   DisplaceStage.Actions.push_back(&PathData.Actions.ShortRangeOn);
  // else
  //   //  DisplaceStage.Actions.push_back(&PathData.Actions.DiagonalAction);
  //   DisplaceStage.Actions.push_back(&PathData.Actions.ShortRange);
  if (PathData.Path.LongRange) {
    if (PathData.Actions.UseRPA) {
      cout<<PathData.Path.CloneStr<<" "<<moveName<<" Adding LongRangeRPA Action"<<endl;
      DisplaceStage.Actions.push_back(&PathData.Actions.LongRangeRPA);
    } else if (PathData.Path.DavidLongRange) {
      cout<<PathData.Path.CloneStr<<" "<<moveName<<" Adding DavidLongRange Action"<<endl;
      DisplaceStage.Actions.push_back(&PathData.Actions.DavidLongRange);
    } else {
      cout<<PathData.Path.CloneStr<<" "<<moveName<<" Adding LongRange Action"<<endl;
      DisplaceStage.Actions.push_back(&PathData.Actions.LongRange);
    }
  }
  for (int i=0; i<activeSpecies.size(); i++) {
    int speciesNum = activeSpecies(i);
    if ((PathData.Actions.NodalActions(speciesNum)!=NULL)) {
      cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<activeSpeciesNames(i)<<" Adding Node Action"<<endl;
      DisplaceStage.Actions.push_back(PathData.Actions.NodalActions(speciesNum));
    }
  }
  // Now construct stage list
  Stages.push_back(&DisplaceStage);

  ActiveParticles.resize(numToMove);
  if (MoveAllParticles) {
    int k = 0;
    for (int i=0; i<activeSpecies.size(); i++) {
      int speciesNum = activeSpecies(i);
      for (int j=0; j<PathData.Path.Species(speciesNum).NumParticles; j++) {
        ActiveParticles(k) = PathData.Path.Species(speciesNum).FirstPtcl + j;
        k += 1;
      }
    }
  }

  NumAttempted = 0;
}



void DisplaceMoveClass::MakeMove()
{
  PathClass &Path = PathData.Path;

  // Next, set timeslices
  Slice1 = 0;
  Slice2 = Path.NumTimeSlices()-1;

  // PROBABLY A SMARTER WAY OF DOING THIS WITH VECTORS (WHAT TO DO IF SAME PARTICLE CHOSEN TWICE)
  Array<bool,1> ptclChecker(PathData.Path.NumParticles());
  ptclChecker = 0;

  if (!MoveAllParticles) {
    for (int j = 0; j < ActiveParticles.size(); j++) {
      while (1) {
        int speciesNum = activeSpecies(Path.Random.CommonInt(activeSpecies.size()));
        int ptclNum = Path.Species(speciesNum).FirstPtcl + Path.Random.CommonInt(Path.Species(speciesNum).NumParticles);
        if (!ptclChecker(ptclNum)) {
          ActiveParticles(j) = ptclNum;
          ptclChecker(ptclNum) = 1;
          //cout << j << " " << speciesNum << "  " << ActiveParticles(j) << endl;
          break;
        }
      }
    }
  }

  // See if we have any with just identity permutation
  int N = PathData.Path.NumParticles();
  Array<int,1> TotalPerm(N), doDisplace(ActiveParticles.size());
  PathData.Path.TotalPermutation(TotalPerm);
  /// Only proc 0 gets TotalPerm
  doDisplace = 0;
  if (Path.Communicator.MyProc() == 0) {
    for (int ptclIndex=0; ptclIndex<ActiveParticles.size(); ptclIndex++) {
      int ptcl = ActiveParticles(ptclIndex);
      if (TotalPerm(ptcl) == ptcl) // Only displace identity permutations
        doDisplace(ptclIndex) = 1;
    }
  }
  PathData.Path.Communicator.Broadcast(0, doDisplace);

  bool MakeTheMove = 0;
  for (int ptclIndex=0; ptclIndex<ActiveParticles.size(); ptclIndex++) {
    if (doDisplace(ptclIndex))
      MakeTheMove = 1;
  }

  // Now call MultiStageClass' MakeMove
  if (MakeTheMove)
    MultiStageClass::MakeMove();
  NumAttempted++;

}
