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

#ifndef NO_PERMUTE_STAGE_CLASS_H
#define NO_PERMUTE_STAGE_CLASS_H

#include "PermuteStage.h"
#include "MultiStage.h"
#include "PermuteTableClass.h"
#include "../Observables/ObservableVar.h"

class NoPermuteStageClass : public PermuteStageClass
{
private:
  int ChooseParticle();
  int NumLevels;
  int NumMoves, NumAccepted;
public:
  double Sample (int &slice1, int &slice2,
		 Array<int,1> &activeParticles);
  bool Attempt (int &slice1, int &slice2, 
		 Array<int,1> &activeParticles, double &prevActionChange);

  void InitBlock(int &slice1,int &slice2);
  NoPermuteStageClass (PathDataClass &pathData, int speciesNum, int numLevels,
		       IOSectionClass &outSection) 
    : PermuteStageClass(pathData, speciesNum, numLevels,outSection)
  {
    // do nothing for now
  }
};




#endif
