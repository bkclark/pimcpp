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


#include "MultiStage.h"

void MultiStageClass::WriteRatio()
{
   list<StageClass*>::iterator stageIter=Stages.begin();
   double prevActionChange = 0.0;
   while (stageIter!=Stages.end()) {
     (*stageIter)->WriteRatio();
     stageIter++;
   }
   MoveClass::WriteRatio();

   // Center of mass movement
   CenterOfMassVar.Write(cm2);
   cm2 = 0.;
}


void MultiStageClass::MakeMove()
{
  bool toAccept = true;
  list<StageClass*>::iterator stageIter=Stages.begin();
  double prevActionChange = 0.0;
  while (stageIter!=Stages.end() && toAccept) {
    toAccept = (*stageIter)->Attempt(Slice1,Slice2,ActiveParticles,prevActionChange);
    stageIter++;
  }

  NumAttempted++;
  if (toAccept)
    Accept();
  else
    Reject();
}


void MultiStageClass::Accept()
{
  PathData.AcceptMove(Slice1,Slice2,ActiveParticles);
  for (list<StageClass*>::iterator stageIter=Stages.begin();stageIter!=Stages.end();stageIter++)
    (*stageIter)->Accept();
  NumAccepted++;

  // Center of mass
  cm2 += PathData.Path.cm2;
}


void MultiStageClass::Reject()
{
  PathData.RejectMove(Slice1,Slice2,ActiveParticles);
  for (list<StageClass*>::iterator stageIter=Stages.begin();stageIter!=Stages.end();stageIter++)
    (*stageIter)->Reject();
}

