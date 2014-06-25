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
  PathData.moveClock++;
  TimesCalled++;
  MakeMove();
}

void MoveClass::WriteRatio()
{
  RatioVar.Write(AcceptanceRatio());
  NumAccepted = TimesCalled = 0;
}

void ParticleMoveClass::SetActiveSpecies (Array<int,1> ActSpecies)
{
  ActiveSpecies.resize(ActSpecies.size());
  ActiveSpecies = ActSpecies;
  ///This calculates the total number of particles
  TotalParticles = 0;
  for (int i=0; i<ActSpecies.size(); i++) {
    int CurrentNumPtcls = PathData.Path.Species(i).NumParticles; 
    TotalParticles += CurrentNumPtcls;
  }
}

/// So do we still want to choose particles by dumping everything
/// into some mapping array from the active particles and dealing 
// with it that way? I think this is doing duplicate stuff in here.
void ParticleMoveClass::ChooseParticles()
{
  for (int i=0; i<NumParticlesToMove; i++) { 
    bool Redundant;
    do {
      MyParticleIndices(i) = PathData.Path.Random.LocalInt(TotalParticles);
      while (PathData.Path.OpenPaths && 
	     MyParticleIndices(i)==(int)(PathData.Path.OpenPtcl)){ 
	MyParticleIndices(i) = PathData.Path.Random.LocalInt(TotalParticles);
      } 
      Redundant = false;
      for (int j=0; j<i; j++){
	if (MyParticleIndices(i) == MyParticleIndices(j)){
	  Redundant = true;
	  break;
	}
      }      
    } while (Redundant); 
  }
  for (int i=0; i<NumParticlesToMove; i++) 
    ActiveParticles(i) = MyParticleIndices(i);
}

/// So do we still want to choose particles by dumping everything
/// into some mapping array from teh active particles and dealing 
// with it that way? I think this is doing duplicate stuff in here.
void ParticleMoveClass::ChooseParticlesOpen()
{
  for (int i=0; i<NumParticlesToMove; i++) { 
    bool Redundant;
    do {
      MyParticleIndices(i) = PathData.Path.Random.LocalInt(TotalParticles);
      Redundant = false;
      for (int j=0; j<i; j++){
	if (MyParticleIndices(i) == MyParticleIndices(j)){
	  Redundant = true;
	  break;
	}
      }      
    } while (Redundant); 
  }
  for (int i=0; i<NumParticlesToMove; i++) 
    ActiveParticles(i) = MyParticleIndices(i);
}
