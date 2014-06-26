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
  SetActiveSpecies (activeSpeciesNames);

  // Whether or not to always accept the move
  bool alwaysAccept;
  if(!in.ReadVar ("AlwaysAccept", alwaysAccept))
    alwaysAccept = false;
  RefSliceShiftStage.alwaysAccept = alwaysAccept;

  // Distance to shift reference slice
  int shiftDistance;
  if(!in.ReadVar ("ShiftDistance", shiftDistance))
    shiftDistance = -1;
  RefSliceShiftStage.shiftDistance = shiftDistance;

  // Construct action list
  for (int i=0; i<ActiveSpecies.size(); i++) {
    int speciesNum = ActiveSpecies(i);
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

bool RefSliceShiftStageClass::Attempt(int &slice1, int &slice2, Array<int,1> &activeParticles, double &prevActionChange)
{
  if(alwaysAccept)
    return true;
  else
    return CommonStageClass::Attempt(slice1,slice2,activeParticles,prevActionChange);
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
  int newRefSlice;
  if(shiftDistance == -1)
    newRefSlice = Path.Random.CommonInt(Path.TotalNumSlices/2);
  else
    newRefSlice = (oldRefSlice + shiftDistance) % Path.TotalNumSlices;
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
