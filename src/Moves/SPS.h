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

#ifndef STOCHASTIC_SWITCHING_POTENTIAL_MOVE_H
#define STOCHASTIC_SWITCHING_POTENTIAL_MOVE_H

#include <list>
#include "MoleculeMoveManager.h"
#include "../Observables/ObservableVar.h"
#include "../Actions/ActionBase.h"
#include "StageClass.h"

/// This move class implements the idea
/// described in C.H. Mak, JCP 122, 214110 (2005)
/// The concept is similar to parallel tempering:
/// A Metropolis-like decision is made to switch
/// to an alternate potential, which is presumably
/// cheaper to evaluate than the "correct" one.

enum ActionMode{FULL, SWITCH};

class SPSClass : public MoleculeMoveStageManagerClass
{
protected:
  list<MolMoveClass*> FullMoveStages;
  list<MolMoveClass*> SwitchMoveStages;
  list<ActionBaseClass*> FullActions;
  list<ActionBaseClass*> SwitchActions;
  // not used yet; might be useful for adjusting actions ptcl-wise
  Array<int,1> PtclKey;
  ActionMode mode;
  double Ts2f, Tf2s;
  int NumSteps;
  int Slice1,Slice2;
  int NumPtcls, NumSlices;
  int NumMakeMove, NumFMakeMove, NumSMakeMove;
  int NumF2SProp, NumS2FProp, NumF2SAcc, NumS2FAcc;
  int NumFAttempted, NumFAccepted, NumSAttempted, NumSAccepted;

public:
  double ComputeEnergy(ActionMode mode);
  void Read(IOSectionClass &io);
  void Accept();
  void Reject();
  virtual void WriteRatio();
  void MakeMove();
  void GetMode(ActionMode& setMode);
  void GetMode(string& setMode);
  bool Attempt(int &slice1, int &slice2, Array<int,1> &activeParticles, double &prevActionChange);
  SPSClass(PathDataClass &pathData, IOSectionClass &outSection) : 
    MoleculeMoveStageManagerClass(pathData,outSection) 
  {
    //do nothing for now
    NumSteps = 0;
  }
};

#endif
