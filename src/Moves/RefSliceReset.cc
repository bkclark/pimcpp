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

#include "RefSliceReset.h"

void RefSliceResetClass::Read(IOSectionClass &in)
{
  // Read in the active species.
  Array<string,1> activeSpeciesNames;
  assert(in.ReadVar ("ActiveSpecies", activeSpeciesNames));
  SetActiveSpecies (activeSpeciesNames);
}

void RefSliceResetClass::MakeMove()
{
  // Reset to reference point position
  SetMode(NEWMODE);
  for (int ptclIndex=0; ptclIndex<ActiveParticles.size(); ptclIndex++) {
    int ptcl = ActiveParticles(ptclIndex);
    for (int slice=0; slice<Path.TotalNumSlices; slice++)
      Path(slice, ptcl) = Path.RefPath(ptcl);
    Path.Permutation(ptcl) = ptcl;
  }

  // Always accepted
  Path.AcceptCopy(0,Path.TotalNumSlices,ActiveParticles);
  NumAccepted++;
  NumAttempted++;
}
