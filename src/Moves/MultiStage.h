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

#ifndef MULTI_STAGE_H
#define MULTI_STAGE_H


#include <list>
#include "MoveBase.h"
#include "../Observables/ObservableVar.h"
#include "../Actions/ActionBase.h"
#include "StageClass.h"

///One method of making a move is to build it out of a set of
///stages. This class allows for such a construct. There is a list of
///stages (in the variable Stages).  By default, this moves iterates
///over the stages, running each one sequentially. Each stage must
///return true or false indicating
// whether or not that stage is accepted. If the stage is
///accepted then the next stage is run. 
class MultiStageClass : public ParticleMoveClass
{
protected:
  list<StageClass*> Stages;
  int NumSteps;
  int NumAttempted;
  int Slice1,Slice2;
  double TimeSpent2;
  ObservableDouble CenterOfMassVar;

  inline double AcceptanceRatio() 
  {
    return (double)(NumAccepted)/(double)NumAttempted;
  }

public:
  double cm2;
  void Read(IOSectionClass &io);
  void Accept();
  void Reject();
  virtual void WriteRatio();
  double NewMoveProb;
  double OldMoveProb;
  ///Why was this MakeMove()=0 and virtual?
  void MakeMove();
  MultiStageClass(PathDataClass &pathData, IOSectionClass &outSection) : 
    CenterOfMassVar("CenterOfMassDrift",IOSection,pathData.Path.Communicator),
    ParticleMoveClass(pathData,outSection) , TimeSpent2(0)
  {
    //do nothing for now
  }
};


#endif
