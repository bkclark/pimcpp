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

#include "OpenEndMove.h"
#include "BisectionStage.h"
#include "EndStage.h"

void OpenEndMoveClass::WriteRatio()
{
  MultiStageClass::WriteRatio();
  //Do nothing for now
}

void OpenEndMoveClass::MakeMove()
{
  //  cerr<<"GRR!!  I'm HERE!"<<endl;
  //  sleep(1000);
  //  cerr<<"Making open end move"<<endl;
  ///WARNING! HACK! HACK! HACK! 
//   if (!PathData.Path.NowOpen){
//     MultiStageClass::Reject();
//     return;
//   }
  ///END HACK!
  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;
  while (stageIter!=Stages.end() && toAccept){
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);
    stageIter++;
  }
  if (toAccept)
    MultiStageClass::Accept();
  else 
    MultiStageClass::Reject();
  //GRR!  MultiStageClass::MoveClass::MakeMove();
}

void OpenEndMoveClass::Read(IOSectionClass &in)
{
  string speciesName;
  StageClass *permuteStage;
  assert (in.ReadVar ("NumLevels", NumLevels));
  SpeciesNum = PathData.Path.OpenSpeciesNum;
  StageClass* endStage;
  endStage=new EndStageClass(PathData,NumLevels,IOSection);
  endStage->Read(in);
  endStage->Actions.push_back(&PathData.Actions.Kinetic);
  endStage->BisectionLevel = NumLevels;
  Stages.push_back (endStage);
  for (int level=NumLevels-1; level>=0; level--) {
    BisectionStageClass *newStage = 
      new BisectionStageClass (PathData, level, IOSection);

    newStage->UseCorrelatedSampling=false;
    newStage->TotalLevels=NumLevels;

    newStage->Actions.push_back(&PathData.Actions.Kinetic);
    // May need later
    //if (level==0)
      //newStage->Actions.push_back(&PathData.Actions.OpenLoopImportance);
    if (level>0) 
      newStage->Actions.push_back(&PathData.Actions.ShortRangeApproximate);
    else if (PathData.Path.OrderN)
      newStage->Actions.push_back(&PathData.Actions.ShortRangeOn);
    else 
      newStage->Actions.push_back(&PathData.Actions.ShortRange);
    if (level == 0) {
      if (PathData.Path.LongRange){
         if (PathData.Actions.UseRPA)
           newStage->Actions.push_back(&PathData.Actions.LongRangeRPA);
         ///If it's David's long range class then do this
         else if (PathData.Path.DavidLongRange){
           newStage->Actions.push_back(&PathData.Actions.DavidLongRange);
         }
         ////
         else
           newStage->Actions.push_back(&PathData.Actions.LongRange);
      }
      if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
	      cerr << "Adding fermion node action for species " 
	     << speciesName << endl;

	      newStage->Actions.push_back(PathData.Actions.NodalActions(SpeciesNum));
      }
    }
    newStage->BisectionLevel = level;
    Stages.push_back (newStage);
  }
}



