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

#include "RefSliceShift.h"

void RefSliceShiftClass::Read(IOSectionClass &in)
{
  string moveName = "RefSliceShift";
  int myProc = PathData.Path.Communicator.MyProc();

  // Read in the active species.
  Array<string,1> activeSpeciesNames;
  assert(in.ReadVar ("ActiveSpecies", activeSpeciesNames));
  activeSpecies.resize(activeSpeciesNames.size());
  for (int i=0; i<activeSpecies.size(); i++)
    activeSpecies(i) = PathData.Path.SpeciesNum(activeSpeciesNames(i));
  SetActiveSpecies (activeSpecies);

  // Set number of particles to move
  int numToMove = 0;
  for (int i=0; i<activeSpecies.size(); i++) {
    int speciesNum = activeSpecies(i);
    numToMove += PathData.Path.Species(speciesNum).NumParticles;
  }
  SetNumParticlesToMove(numToMove);
  int k = 0;
  for (int i=0; i<activeSpecies.size(); i++) {
    int speciesNum = activeSpecies(i);
    for (int j=0; j<PathData.Path.Species(speciesNum).NumParticles; j++) {
      ActiveParticles(k) = PathData.Path.Species(speciesNum).FirstPtcl + j;
      k += 1;
    }
  }

  // Construct action list
  for (int i=0; i<activeSpecies.size(); i++) {
    int speciesNum = activeSpecies(i);
    if ((PathData.Actions.NodalActions(speciesNum)!=NULL)) {
      if (myProc == 0)
        cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<activeSpeciesNames(i)<<" Adding Node Action"<<endl;
      RefSliceShiftStage.Actions.push_back(PathData.Actions.NodalActions(speciesNum));
    }
  }
  // Now construct stage list
  Stages.push_back(&RefSliceShiftStage);

  // Reset counter
  NumAttempted = 0;

}

void RefSliceShiftStageClass::Accept()
{
  CommonStageClass::Accept();
  Path.RefPath.AcceptCopy();
  if (Path.StoreNodeDist)
    Path.NodeDist.AcceptCopy();
  if (Path.StoreNodeDet)
    Path.NodeDet.AcceptCopy();
}

void RefSliceShiftStageClass::Reject()
{
  CommonStageClass::Reject();
  Path.RefSlice = oldRefSlice;
  Path.BroadcastRefPath();
  Path.RefPath.RejectCopy();
  if (Path.StoreNodeDist)
    Path.NodeDist.RejectCopy();
  if (Path.StoreNodeDet)
    Path.NodeDet.RejectCopy();
}

double RefSliceShiftStageClass::Sample (int &slice1, int &slice2, Array<int,1> &activeParticles)
{
  PathClass &Path = PathData.Path;

  // Pick new reference slice and broadcast it to all the other processors
  oldRefSlice = Path.GetRefSlice();
  int newRefSlice = Path.Random.CommonInt(Path.TotalNumSlices/2);
  Path.RefSlice = newRefSlice;
  Path.BroadcastRefPath();

  // And return log sample probability ratio
  return log(1.0);
}

void RefSliceShiftClass::MakeMove()
{
  PathClass &Path = PathData.Path;

  // Next, set timeslices
  Slice1 = 0;
  Slice2 = Path.NumTimeSlices()-1;

  MultiStageClass::MakeMove();
}
