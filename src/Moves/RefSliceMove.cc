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

#include "RefSliceMove.h"
#include "NoPermuteStage.h"

void RefSliceMoveClass::WriteRatio()
{
  MultiStageClass::WriteRatio();
}

void RefSliceMoveClass::Read(IOSectionClass &in)
{
  string moveName = "RefSliceMove";
  string permuteType, speciesName;

  assert (in.ReadVar ("NumLevels", NumLevels));
  assert (in.ReadVar ("Species", speciesName));
  SpeciesNum = PathData.Path.SpeciesNum (speciesName);

  if(!in.ReadVar ("DoSlaveMoves", DoSlaveMoves))
    DoSlaveMoves = false;

  /// Set up permutation
  assert (in.ReadVar ("PermuteType", permuteType));
  if (permuteType == "TABLE")
    PermuteStage = new TablePermuteStageClass (PathData, SpeciesNum, NumLevels, IOSection);
  else if (permuteType == "NONE")
    PermuteStage = new NoPermuteStageClass(PathData, SpeciesNum, NumLevels, IOSection);
  else {
    cerr << "Unrecognized PermuteType, """ << permuteType << """\n";
    exit(EXIT_FAILURE);
  }
  PermuteStage->Read (in);
  Stages.push_back (PermuteStage);

  for (int level=NumLevels-1; level>=0; level--) {
    BisectionStageClass *newStage = new BisectionStageClass (PathData, level, IOSection);
    newStage -> TotalLevels = NumLevels;
    newStage -> BisectionLevel = level;
    newStage -> UseCorrelatedSampling=false;
    cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding Kinetic Action"<<endl;
    newStage -> Actions.push_back(&PathData.Actions.Kinetic);
    if (level == 0) {
      cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding ShortRange Action"<<endl;
      newStage -> Actions.push_back(&PathData.Actions.ShortRange);
      if (PathData.Path.LongRange) {
        if (PathData.Actions.UseRPA){
          cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding LongRangeRPA Action"<<endl;
          newStage -> Actions.push_back(&PathData.Actions.LongRangeRPA);
        } else if (PathData.Path.DavidLongRange) {
          cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding DavidLongRange Action"<<endl;
          newStage -> Actions.push_back(&PathData.Actions.DavidLongRange);
        } else {
          cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding LongRangeAction"<<endl;
          newStage -> Actions.push_back(&PathData.Actions.LongRange);
        }
      }
      /// No need to calculate this twice!
      //
      //if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
      //  cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding Node Action"<<endl;
      //  newStage -> Actions.push_back(PathData.Actions.NodalActions(SpeciesNum));
      //}
    }
    Stages.push_back (newStage);
  }
  // Add the second stage of the permutation step
  Stages.push_back (PermuteStage);
}


bool RefSliceMoveClass::NodeCheck()
{
  if (PathData.Actions.NodalActions(SpeciesNum) != NULL) {
    PathClass &Path=PathData.Path;
    // Broadcast the new reference path to all the other processors
    PathData.Path.BroadcastRefPath();

    // Calculate local nodal action
    SetMode(OLDMODE);
    double oldLocalNode = PathData.Actions.NodalActions(SpeciesNum) -> Action(0, Path.NumTimeSlices()-1, ActiveParticles, 0);
    SetMode(NEWMODE);
    double newLocalNode = PathData.Actions.NodalActions(SpeciesNum) -> Action(0, Path.NumTimeSlices()-1, ActiveParticles, 0);

    // Do global sum over processors
    double localChange = newLocalNode - oldLocalNode;
    double globalChange = PathData.Path.Communicator.AllSum (localChange);
    bool toAccept = (-globalChange)>=log(PathData.Path.Random.Common());

    // Check if Broken
    if ((abs(newLocalNode) > 1e50 || abs(oldLocalNode) > 1e50) && toAccept) {
      cerr << Path.CloneStr << " Broken NodeCheck!!! " <<Path.NumTimeSlices()-1<< " " <<  toAccept << " " << PathData.Path.Species(SpeciesNum).Name << " " << PathData.Path.GetRefSlice() << " " << PathData.Path.SliceOwner(PathData.Path.GetRefSlice()) << " " << oldLocalNode << " " << newLocalNode << " " << localChange << " " << globalChange << endl;
      if (abs(newLocalNode > 1e50)) {
        toAccept = 0;
        assert(1==2);
      } else {
        toAccept = 1;
      }
    }
    //else {
    //  cout << Path.CloneStr << " NodeCheck " <<Path.NumTimeSlices()-1<< " " <<  toAccept << " " << PathData.Path.Species(SpeciesNum).Name << " " << PathData.Path.GetRefSlice() << " " << PathData.Path.SliceOwner(PathData.Path.GetRefSlice()) << " " << oldLocalNode << " " << newLocalNode << " " << localChange << " " << globalChange << endl;
    //}

    return toAccept;
  }
  else
    return true;
}


/// This version is for the processor with the reference slice
void RefSliceMoveClass::MakeMoveMaster()
{
  // static bool firstTime=true;
  // cerr<<"MAKING REF SLICE MOVE"<<endl;
  // cerr<<"NODE CHECK IS "<<NodeCheck()<<endl;
  // cerr<<"DONE NODECHECK"<<endl;
  // if (!firstTime)
  //   return;
  PathClass &Path = PathData.Path;
  int myProc = PathData.Path.Communicator.MyProc();
  int firstSlice, lastSlice;
  Path.SliceRange (myProc, firstSlice, lastSlice);
  // Choose time slices
  int localRef = Path.GetRefSlice() - firstSlice;
  int bisectSlices = (1 << NumLevels);
  int minStart = max(0, localRef - bisectSlices);
  int maxStart = min(Path.NumTimeSlices()-1-bisectSlices, localRef);
  maxStart = max(0, maxStart);
  Slice1 = Path.Random.LocalInt(maxStart-minStart) + minStart;
  Slice2 = Slice1 + bisectSlices;
  assert (Slice1 >= 0);
  assert ((localRef >= Slice1) && (localRef <= Slice2));
  // Move the join out of the way.
  PathData.MoveJoin(Slice2);
  ActiveParticles.resize(1);
  ActiveParticles(0) = -1;
  // Go through local stages
  bool toAccept = true;
  list<StageClass*>::iterator stageIter = Stages.begin();
  double prevActionChange = 0.0;
  int stageCounter = 0;
  ((PermuteStageClass*)PermuteStage) -> InitBlock(Slice1, Slice2);
  while (stageIter != Stages.end() && toAccept){
    //cerr<<"Attempting: level("<<(*stageIter)->BisectionLevel<<")"endl;
    toAccept = (*stageIter) -> Attempt(Slice1, Slice2, ActiveParticles, prevActionChange);
    stageCounter++;
    stageIter++;
  }
  // Broadcast acceptance or rejection
  SetMode (NEWMODE);
  int accept = toAccept ? 1 : 0;
  PathData.Path.Communicator.Broadcast (myProc, accept);

  // Now, if we accept local stages, move on to global nodal decision
  if (toAccept) {
    if (NodeCheck()) {
      NodeAccept++;
      Accept();
      Path.RefPath.AcceptCopy();
      if (Path.UseNodeDist)
        Path.NodeDist.AcceptCopy();
      if (Path.UseNodeDet)
        Path.NodeDet.AcceptCopy();
      //cerr<<Path.CloneStr << " ACCEPTING REF SLICE MOVE"<<endl;
    }
    else {
      NodeReject++;
      Reject();
      Path.RefPath.RejectCopy();
      if (Path.UseNodeDist)
        Path.NodeDist.RejectCopy();
      if (Path.UseNodeDet)
        Path.NodeDet.RejectCopy();
      //cerr<<Path.CloneStr << " REJECTING REF SLICE MOVE BY NODECHECK"<<endl;
    }
  }
  // Otherwise, reject the whole move
  else {
    NodeReject++;
    Reject();
    Path.RefPath.RejectCopy();
    if (Path.UseNodeDist)
      Path.NodeDist.RejectCopy();
    if (Path.UseNodeDet)
      Path.NodeDet.RejectCopy();
    //cerr<<Path.CloneStr << " REJECTING REF SLICE MOVE BY LOCAL REJECT"<<endl;
  }

}


/// This version is for processors that do not own the reference slice 
void RefSliceMoveClass::MakeMoveSlave()
{
  PathClass &Path=PathData.Path;
  int myProc = PathData.Path.Communicator.MyProc();
  int master = Path.SliceOwner (Path.GetRefSlice());

  if (DoSlaveMoves) {
    /// Choose time slices for local bisections
    /// We don't have the reference slice, so we can be cavalier
    int firstSlice, lastSlice;
    Path.SliceRange (myProc, firstSlice, lastSlice);
    // Choose time slices
    int bisectSlices = (1 << NumLevels);
    int maxStart = Path.NumTimeSlices()-bisectSlices;
    int slice1 = Path.Random.LocalInt(maxStart);
    int slice2 = slice1 + bisectSlices;
    ActiveParticles.resize(1);
    ActiveParticles(0) = (PathData.Path.Random.LocalInt(PathData.Path.Species(SpeciesNum).NumParticles)) + PathData.Path.Species(SpeciesNum).FirstPtcl;
    bool toAccept = true;
    list<StageClass*>::iterator stageIter = Stages.begin();
    double prevActionChange = 0.0;
    int stageCounter = 0;
    while (stageIter != Stages.end() && toAccept){
      toAccept = (*stageIter) -> Attempt(slice1, slice2, ActiveParticles, prevActionChange);
      stageCounter++;
      stageIter++;
    }
  }

  int accept;
  /// Receive broadcast from Master.
  SetMode (NEWMODE);
  PathData.Path.Communicator.Broadcast (master, accept);
  if (accept==1) {
    if (NodeCheck()) {
      NodeAccept++;
      //Accept();
      Path.RefPath.AcceptCopy();
      if (Path.UseNodeDist)
        Path.NodeDist.RejectCopy();
      if (Path.UseNodeDet)
        Path.NodeDet.RejectCopy();
      //cerr<<Path.CloneStr << " ACCEPTING REF SLICE MOVE"<<endl;
    } else {
      NodeReject++;
      //Reject();
      Path.RefPath.RejectCopy();
      if (Path.UseNodeDist)
        Path.NodeDist.RejectCopy();
      if (Path.UseNodeDet)
        Path.NodeDet.RejectCopy();
      //cerr<<Path.CloneStr << " REJECTING REF SLICE MOVE BY NODECHECK"<<endl;
    }
  } else {
    NodeReject++;
    //Reject();
    Path.RefPath.RejectCopy();
    if (Path.UseNodeDist)
      Path.NodeDist.RejectCopy();
    if (Path.UseNodeDet)
      Path.NodeDet.RejectCopy();
    //cerr<<Path.CloneStr << " REJECTING REF SLICE MOVE BY LOCAL REJECT"<<endl;
  }

}


void RefSliceMoveClass::MakeMove()
{
  //cerr << "RefSlice Move" << endl;
  PathClass &Path = PathData.Path;
  MasterProc = Path.SliceOwner (Path.GetRefSlice());
  //cerr<<"Starting RefSlice move."<<endl;
  if (PathData.Path.Communicator.MyProc() == MasterProc){
    //cout<<"MakeMoveMaster();"<<endl;
    MakeMoveMaster();
  } else {
    //cout<<"MakeMoveSlave();"<<endl;
    MakeMoveSlave();
  }
  if ((NodeAccept+NodeReject) % 10000 == 9999)
    cout << PathData.Path.CloneStr <<  " Node accept ratio = " << (double)NodeAccept/(double)(NodeAccept+NodeReject);
  //cerr<<"Finished RefSlice move."<<endl;
}
