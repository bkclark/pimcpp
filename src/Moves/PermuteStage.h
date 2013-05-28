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

#ifndef PERMUTE_STAGE_CLASS_H
#define PERMUTE_STAGE_CLASS_H

#include "MultiStage.h"
#include "PermuteTableClass.h"
#include "PermuteTableOnClass.h"
#include "../Observables/ObservableVar.h"

class PermuteStageClass : public LocalStageClass
{
protected:
  int SpeciesNum, NumLevels;
public:
  /// This function is responsible for setting the activeParticles if
  /// activeParticles(0) == -1; It is called a second time after the
  /// bisection stages.  We impose the condtion that the product of
  /// the two return values must be equal to true transition
  /// probablity ratio.  This allows avoiding computing the reverse
  /// probability if the move is rejected before this stage.

  virtual double Sample (int &slice1, int &slice2, 
			 Array<int,1> &activeParticles) = 0;
  virtual bool Attempt (int &slice1, int &slice2, 
			Array<int,1> &activeParticles, double &prevActionChange) = 0;
  virtual void InitBlock(int &slice1,int &slice2);
  virtual void Read (IOSectionClass &in);
  virtual void Accept();
  virtual void Reject();
  PermuteStageClass(PathDataClass &pathData, int speciesNum,
		    int numLevels,IOSectionClass &outSection) : 
    LocalStageClass (pathData,outSection), 
    SpeciesNum (speciesNum), 
    NumLevels(numLevels)

  {
    // do nothing for now 
  }
};



class TablePermuteStageClass : public PermuteStageClass
{
private:
  PermuteTableClass Table1, Table2;
  PermuteTableClass *Forw, *Rev;
  ObservableVecDouble1 AcceptanceRatioVar;
  ObservableVecInt1 AcceptanceTotalVar;
  Array<int,1> NumAccepted;
  Array<int,1> NumAttempted;
  bool NeedToRebuildTable;
  bool HaveBeenAcceptedOrRejected;
  double forwT;
  bool zFocus;
  bool NeedToUpdateHTableOnReject;
public:
  /// This function will construct a new permutation if
  /// activeParticles is set to the array, [ -1 ];  In this case,
  /// it will set the activeParticles.  It returns 1.0
  /// as the transition probability ratio at this time.
  /// If called with a valid set of particles in activeParticles,
  /// it changes nothing, but returns the transition probability
  /// ratio for the sampling.  This is so we can avoid calculating
  /// that ratio if the move is rejected, saving time.  Thus, this
  /// function is called twice during a successful multistage move.
 
  void InitBlock(int &slice1,int &slice2);
  void Read (IOSectionClass &in);
  double Sample (int &slice1, int &slice2,
		 Array<int,1> &activeParticles);
  bool Attempt (int &slice1, int &slice2, 
		   Array<int,1> &activeParticles, double &prevActionChange);
  void Accept ();
  void Reject();
  void WriteRatio();
  TablePermuteStageClass (PathDataClass &pathData, int speciesNum, int numLevels,IOSectionClass &outSection) : 
    PermuteStageClass(pathData, speciesNum, numLevels,outSection),
    Table1(pathData), Table2(pathData),
    AcceptanceRatioVar("Acceptance Ratio",OutSection,PathData.Path.Communicator),    AcceptanceTotalVar("Perms Tried",OutSection,PathData.Path.Communicator)

  {
    Forw = &Table1;
    Rev  = &Table2;
    NumAccepted.resize(4);
    NumAttempted.resize(4);
    NumAccepted=0;
    NumAttempted=0;
    NeedToRebuildTable=true;
  }
};


#endif
