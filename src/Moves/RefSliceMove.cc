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
  double ratio = (double)NodeAccept/(double)(NodeAccept+NodeReject);
  RatioVar.Write(ratio);  
}

void RefSliceMoveClass::Read(IOSectionClass &in)
{
  string permuteType, speciesName;

  assert (in.ReadVar ("NumLevels", NumLevels));
  assert (in.ReadVar ("Species", speciesName));
  SpeciesNum = PathData.Path.SpeciesNum (speciesName);

  /// Set up permutation
  assert (in.ReadVar ("PermuteType", permuteType));
  if (permuteType == "TABLE") 
    PermuteStage = new TablePermuteStageClass (PathData,SpeciesNum,NumLevels,
					       IOSection);
  else if (permuteType == "NONE") 
    PermuteStage = new NoPermuteStageClass(PathData, SpeciesNum, NumLevels,
					   IOSection);
  else {
    cerr << "Unrecognized PermuteType, """ << permuteType << """\n";
    exit(EXIT_FAILURE);
  }
  PermuteStage->Read (in);
  Stages.push_back (PermuteStage);
  
  for (int level=NumLevels-1; level>=0; level--) {
    BisectionStageClass *newStage = new BisectionStageClass (PathData, level,
							     IOSection);
    newStage->Actions.push_back(&PathData.Actions.Kinetic);
    newStage->UseCorrelatedSampling=false;
    ///HACK!    newStage->Actions.push_back(&PathData.Actions.ShortRange);
    if (level == 0) {
      if (PathData.Path.DavidLongRange){
	cerr<<"Pushing david long range at lowetst level"<<endl;
	newStage->Actions.push_back(&PathData.Actions.DavidLongRange);
      }
      else  if (PathData.Path.LongRange) 
	if (PathData.Actions.UseRPA)
	  newStage->Actions.push_back(&PathData.Actions.LongRangeRPA);
	else
	  newStage->Actions.push_back(&PathData.Actions.LongRange);
      
//       for (int i=0; i<PathData.Actions.NodalActions.size(); i++)
// 	newStage->Actions.push_back(PathData.Actions.NodalActions(i));

      if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
	cerr << "Adding fermion node action for species " 
	     << speciesName << endl;
	newStage->Actions.push_back(PathData.Actions.NodalActions(SpeciesNum));
      }
    }
    newStage->BisectionLevel = level;
    newStage->TotalLevels=NumLevels;
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
    double oldLocalNode = 
      PathData.Actions.NodalActions(SpeciesNum)->Action
      (0, Path.NumTimeSlices()-1, ActiveParticles, 0);
    SetMode(NEWMODE);
    double newLocalNode = 
      PathData.Actions.NodalActions(SpeciesNum)->Action
      (0, Path.NumTimeSlices()-1, ActiveParticles, 0);
    
//     if (newLocalNode > 1.0e20) 
//       cerr << "Local node crossing in NodeCheck().\n";

    // Do global sum over processors
    double localChange = newLocalNode - oldLocalNode;
    double globalChange = PathData.Path.Communicator.AllSum (localChange);
    bool toAccept = (-globalChange)>=log(PathData.Path.Random.Common()); 

//     fprintf (stderr, "old = %1.12e\n", oldLocalNode);
//     fprintf (stderr, "new = %1.12e\n", newLocalNode);
//     fprintf (stderr, "globalChange = %1.12e\n", globalChange);
//     fprintf (stderr, "localchange = %1.12e\n", localChange);
//     fprintf (stderr, "globalchange = %1.12e\n", globalChange);
//     cerr << "toAccept = " << (toAccept ? "true\n" : "false\n");

    return toAccept;
  }
  else
    return true;
}


/// This version is for the processor with the reference slice
void RefSliceMoveClass::MakeMoveMaster()
{
//    static bool firstTime=true;
//  cerr<<"MAKING REF SLICE MOVE"<<endl;
//  cerr<<"NODE CHECK IS "<<NodeCheck()<<endl;
//  cerr<<"DONE NODECHECK"<<endl;
//  if (!firstTime)
//    return;
  PathClass &Path=PathData.Path;
  int myProc = PathData.Path.Communicator.MyProc();
  int firstSlice, lastSlice;
  Path.SliceRange (myProc, firstSlice, lastSlice);
  // Choose time slices
  int localRef = Path.GetRefSlice() - firstSlice;
  int bisectSlices = (1 << NumLevels);
  int minStart = max(0, localRef - bisectSlices);
  int maxStart = 
    min (Path.NumTimeSlices()-1-bisectSlices, localRef);
  maxStart = max (0, maxStart);
  Slice1 = Path.Random.LocalInt(maxStart-minStart) + minStart;
  Slice2 = Slice1 + bisectSlices;
  assert (Slice1 >= 0);
  assert ((localRef >= Slice1) && (localRef <= Slice2));
  // Move the join out of the way.
  PathData.MoveJoin(Slice2);
  ActiveParticles.resize(1);
  ActiveParticles(0) = -1;
  // Go through local stages
  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;
  int stageCounter = 0;
  ((PermuteStageClass*)PermuteStage)->InitBlock(Slice1, Slice2);
  while (stageIter!=Stages.end() && toAccept){
    cerr<<"Attempting level"<<(*stageIter)->BisectionLevel<<endl;
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);
    stageCounter++;
    stageIter++;
  }
  // Broadcast acceptance or rejection 
  int accept = toAccept ? 1 : 0;
  PathData.Path.Communicator.Broadcast (myProc, accept);

  // Now, if we accept local stages, move on to global nodal
  // decision. 
  if (toAccept) {
//     if ((NodeAccept+NodeReject) % 1000 == 999)
//       fprintf (stderr, "Node accept ratio = %5.3f\n",
// 	       (double)NodeAccept/(double)(NodeAccept+NodeReject));

    if (NodeCheck()) {

      NodeAccept++;
      Accept();
      Path.RefPath.AcceptCopy();
      cerr<<"ACCEPTNG REF SLICE MOVE"<<endl;
    }
//    else if (firstTime){
//      firstTime=false;
//      NodeAccept++;
//      Accept();
//      Path.RefPath.AcceptCopy();
//      cerr<<"ACCEPTNG REF SLICE MOVE FOR FIRSTTIME"<<endl;
//    }
    else {
      NodeReject++;
      Reject();
      Path.RefPath.RejectCopy();
      cerr<<"REJECTING REF SLICE MOVE BY NODECHECK"<<endl;
    }
  }
  // Otherwise, reject the whole move
  else {
    Reject();
    cerr<<"REJECTING REF SLICE MOVE BY LOCAL REJECt"<<endl;
    Path.RefPath.RejectCopy();
  }

}


/// This version is for processors that do no own the reference slice 
void RefSliceMoveClass::MakeMoveSlave()
{
  PathClass &Path=PathData.Path;
  int myProc = PathData.Path.Communicator.MyProc();
  int master = Path.SliceOwner (Path.GetRefSlice());

//   /// Choose time slices for local bisections
//   /// We don't have the reference slice, so we can be cavalier
//   int firstSlice, lastSlice;
//   Path.SliceRange (myProc, firstSlice, lastSlice);
//   // Choose time slices
//   int bisectSlices = (1 << NumLevels);
//   int maxStart = Path.NumTimeSlices()-bisectSlices;
//   int slice1 = Path.Random.LocalInt(maxSlice);
//   int slice2 = slice1 + bisectSlices;

  int accept;
  /// Receive broadcast from Master.
  SetMode (NEWMODE);
  PathData.Path.Communicator.Broadcast (master, accept);
  if (accept==1) {
    if (NodeCheck()) {
      Path.RefPath.AcceptCopy();
      NodeAccept++;
    }
    else {
      Path.RefPath.RejectCopy();
      NodeReject++;
    }
  }
}


void RefSliceMoveClass::MakeMove()
{
  PathClass &Path=PathData.Path;
  MasterProc = Path.SliceOwner (Path.GetRefSlice());
  if (PathData.Path.Communicator.MyProc() == MasterProc)
    MakeMoveMaster();
  else
    MakeMoveSlave();
  if ((NodeAccept+NodeReject) % 1000 == 999)
    fprintf (stderr, "Node accept ratio = %5.3f\n",
	     (double)NodeAccept/(double)(NodeAccept+NodeReject));
}
