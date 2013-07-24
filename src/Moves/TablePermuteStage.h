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

#ifndef TABLE_PERMUTE_STAGE_CLASS_H
#define TABLE_PERMUTE_STAGE_CLASS_H

#include "PermuteStageClass.h"

class TablePermuteStageClass : public PermuteStageClass
{
private:
  int SpeciesNum;
  bool HaveBeenAcceptedOrRejected;
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
  double forwT;
  bool zFocus;
  double Sample (int &slice1,int &slice2, Array<int,1> &changedParticles);
  TablePermuteStageClass (PathDataClass &pathData, int speciesNum, int numLevels)
    : PermuteStageClass(pathData, speciesNum, numLevels)
  {
    // do nothing for now
  }

};

#endif
