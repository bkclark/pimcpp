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

#include "NodalModelMove.h"
#include "sys/time.h"

void NodalModelMoveClass::WriteRatio()
{

  NumAttemptedVar.Write(NumAttempted);
  NumAcceptedVar.Write(NumAccepted);
  MultiStageClass::WriteRatio();
}

void NodalModelStageClass::WriteRatio()
{
  AttemptArrayVar.Write(AttemptArray);
  AttemptArray = 0;
  AttemptArrayVar.Flush();
  AcceptArrayVar.Write(AcceptArray);
  AcceptArray = 0;
  AcceptArrayVar.Flush();
}

void NodalModelStageClass::Accept()
{
  CommonStageClass::Accept();

  for (int i=0; i<activeSpecies.size(); i++) {
    int SpeciesNum = activeSpecies(i);
    int newModel = PathData.Actions.NodalActions(SpeciesNum) -> GetModel();
    AcceptArray(oldModel,newModel) += 1;
    AttemptArray(oldModel,newModel) += 1;
  }

  Path.RefPath.AcceptCopy();
  if (Path.StoreNodeDist)
    Path.NodeDist.AcceptCopy();
  if (Path.StoreNodeDet)
    Path.NodeDet.AcceptCopy();
}


void NodalModelStageClass::Reject()
{

  for (int i=0; i<activeSpecies.size(); i++) {
    int SpeciesNum = activeSpecies(i);
    int newModel = PathData.Actions.NodalActions(SpeciesNum) -> GetModel();
    AttemptArray(oldModel,newModel) += 1;
    if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
      PathData.Actions.NodalActions(SpeciesNum) -> ChangeModel(oldModel);
    }
  }

  CommonStageClass::Reject();
  Path.RefPath.RejectCopy();
  if (Path.StoreNodeDist)
    Path.NodeDist.RejectCopy();
  if (Path.StoreNodeDet)
    Path.NodeDet.RejectCopy();
}


double NodalModelStageClass::Sample (int &slice1, int &slice2, Array<int,1> &activeParticles)
{
  for (int i=0; i<activeSpecies.size(); i++) {
    int SpeciesNum = activeSpecies(i);
    if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
      oldModel = PathData.Actions.NodalActions(SpeciesNum) -> GetModel();
      int NumModels = PathData.Actions.NodalActions(SpeciesNum) -> GetNumModels();
      int newModel = PathData.Path.Random.CommonInt(NumModels); /// get model
      PathData.Actions.NodalActions(SpeciesNum) -> ChangeModel(newModel);
    }
  }

  // And return log sample probability ratio
  return log(1.0);
}


void NodalModelMoveClass::Read (IOSectionClass &in)
{
  int myProc = PathData.Path.Communicator.MyProc();

  string moveName = "NodalModel";

  // Read in the active species.
  Array<string,1> activeSpeciesNames;
  assert(in.ReadVar ("ActiveSpecies", activeSpeciesNames));
  SetActiveSpecies(activeSpeciesNames);

  int NumModels;
  for (int i=0; i<ActiveSpecies.size(); i++) {
    int SpeciesNum = ActiveSpecies(i);
    if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
      if (myProc == 0)
        cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<activeSpeciesNames(i)<<" Adding Node Action"<<endl;
      NodalModelStage.Actions.push_back(PathData.Actions.NodalActions(SpeciesNum));
    }
    NumModels = PathData.Actions.NodalActions(SpeciesNum) -> GetNumModels();
  }

  NodalModelStage.activeSpecies.resize(ActiveSpecies.size());
  for (int i=0; i<ActiveSpecies.size(); i++) {
    NodalModelStage.activeSpecies(i) = ActiveSpecies(i);
    int SpeciesNum = ActiveSpecies(i);
    if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
      NodalModelStage.oldModel = PathData.Actions.NodalActions(SpeciesNum) -> GetModel();
    }
  }

  NodalModelStage.AcceptArray.resize(NumModels,NumModels);
  NodalModelStage.AttemptArray.resize(NumModels,NumModels);

  // Now construct stage list
  Stages.push_back(&NodalModelStage);

  ActiveParticles.resize(0);
}


void NodalModelMoveClass::MakeMove()
{
  PathClass &Path = PathData.Path;

  // Next, set timeslices
  Slice1 = 0;
  Slice2 = Path.NumTimeSlices()-1;

  // Now call MultiStageClass' MakeMove
  MultiStageClass::MakeMove();
  NumAttempted++;

}
