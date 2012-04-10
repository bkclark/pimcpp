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

#include "CouplingMove.h"
#include "BisectionStage.h"
#include "CouplingStage.h"

void CouplingMoveClass::WriteRatio()
{
  MultiStageClass::WriteRatio();
  //Do nothing for now
}

void CouplingMoveClass::MakeMove()
{

  bool toAccept=true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange=0.0;
  int stageCount=0;
  while (stageIter!=Stages.end() && toAccept){
    //    cerr<<"This is stage "<<stageCount<<endl;
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,
				     ActiveParticles,prevActionChange);
    stageIter++;
    //    cerr<<endl;
  }
  if (toAccept)
    MultiStageClass::Accept();
  else 
    MultiStageClass::Reject();
}

void CouplingMoveClass::Read(IOSectionClass &in)
{
  StageClass* coupleStage;
  coupleStage=new CouplingStageClass(PathData,NumLevels,IOSection);
  coupleStage->Read(in);
  coupleStage->Actions.push_back(&PathData.Actions.Kinetic);
  coupleStage->Actions.push_back(&PathData.Actions.ShortRange);
  Stages.push_back (coupleStage);
}



