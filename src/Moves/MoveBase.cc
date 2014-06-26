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

#include "MoveBase.h"
#include "time.h"

void MoveClass::DoEvent()
{
  TimesCalled++;
  MakeMove();
}

void MoveClass::WriteRatio()
{
  RatioVar.Write(AcceptanceRatio());
  NumAccepted = NumAttempted = 0;
}

void ParticleMoveClass::SetActiveSpecies (Array<string,1> activeSpeciesNames)
{
  ActiveSpecies.resize(activeSpeciesNames.size());
  for (int i=0; i<activeSpeciesNames.size(); i++)
    ActiveSpecies(i) = Path.SpeciesNum(activeSpeciesNames(i));

  // Set number of particles to move
  int numToMove = 0;
  for (int i=0; i<ActiveSpecies.size(); i++) {
    int speciesNum = ActiveSpecies(i);
    numToMove += Path.Species(speciesNum).NumParticles;
  }
  ActiveParticles.resize(numToMove);
  int k = 0;
  for (int i=0; i<ActiveSpecies.size(); i++) {
    int speciesNum = ActiveSpecies(i);
    for (int j=0; j<Path.Species(speciesNum).NumParticles; j++) {
      ActiveParticles(k) = Path.Species(speciesNum).FirstPtcl + j;
      k += 1;
    }
  }

}
