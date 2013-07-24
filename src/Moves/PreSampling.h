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

#ifndef PRE_SAMPLING_MOVE_H
#define PRE_SAMPLING_MOVE_H

#include <list>
#include "MoleculeMoveManager.h"
#include "../Observables/ObservableVar.h"
#include "../Actions/ActionBase.h"
#include "StageClass.h"

/// This move is intended in situations where you
/// have a ``cheap'' potential that correlates
/// with a costly potential whose distribution you
/// want to sample.  A markov chain is generated
/// using the cheap potential.  This series of moves
/// is followed by an accept/reject based on the
/// expensive potential.
/// See L.D. Gelb, JCP 118, 7747 (2003) for a
/// derivation
class PreSamplingClass : public MoleculeMoveStageManagerClass
{
protected:
  list<StageClass*> PreStages;
  list<ActionBaseClass*> PreActions;
  list<ActionBaseClass*> FinalActions;
  //StageClass* FinalStage;
  int TotalNumPreSteps;
  double PreDeltaAction;
  bool FirstMove;
  int NumSteps;
  int NumPreAccept, NumFinalAccept;
  int Slice1,Slice2;
  int NumPtcls, NumSlices;
  Array<dVec,2> InitialPath;
  Array<int,1> NumPreSampleSteps;
  Array<int,1> NumPreSampleAccept;

public:
  void Read(IOSectionClass &io);
  void Accept();
  void Reject();
  double StageAction(std::list<ActionBaseClass*> ActionList, int startSlice,int endSlice, const Array<int,1> &changedParticles);
  virtual void WriteRatio();
  double NewMoveProb;
  double OldMoveProb;
  void MakeMove();
  void StoreInitialPath();
  void AssignInitialPath();
  PreSamplingClass(PathDataClass &pathData, IOSectionClass &outSection) : 
    MoleculeMoveStageManagerClass(pathData,outSection) 
  {
    //do nothing for now
    NumSteps = 0;
    NumFinalAccept = 0;
  }
};

class PreSampleDummy : public DummyEvaluate 
{
  list<ActionBaseClass*> PreActions;
  int numPre, numFinal;
  int toRead, startI;
  int numMol, molIndex;

  public:
  double PreSampleAction(int startSlice,int endSlice,
		     const Array<int,1> &changedParticles);
  bool Attempt(int &slice1, int &slice2, 
			      Array<int,1> &activeParticles,
			      double &prevActionChange);
  void Read (IOSectionClass &in);

  PreSampleDummy(PathDataClass& PathData, IOSectionClass& IO, int actionsToRead, int startIndex);

};


#endif
