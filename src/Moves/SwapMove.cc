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

#include "SwapMove.h"

double SwapStageClass::Sample (int &slice1, int &slice2,
				   Array<int,1> &activeParticles)
{
  int helium3Particle=PathData.Path.Species(0).FirstPtcl;
  assert(PathData.Path.Permutation(helium3Particle)==helium3Particle);  
  
  SetMode(NEWMODE);
  dVec temp;
  for (int slice=0; slice<PathData.NumTimeSlices(); slice++){
    temp=PathData.Path(slice, activeParticles(0));
    PathData.Path(slice,activeParticles(0))=PathData.Path(slice, helium3Particle);
    PathData.Path(slice, helium3Particle) =temp;
  }
  // And return sample probability ratio
    return 1.0/(PathData.Path.Species(1).LastPtcl-PathData.Path.Species(1).FirstPtcl+1);
 
}

void
SwapMoveClass::Read (IOSectionClass &in)
{
  assert(in.ReadVar ("Sigma", Sigma));
  SwapStage.Sigma = Sigma;
  Array<string,1> activeSpeciesNames;

  int numToMove;
  assert(in.ReadVar("NumToMove", numToMove));
  SetNumParticlesToMove(numToMove);

  // Read in the active species.
  assert(in.ReadVar ("ActiveSpecies", activeSpeciesNames));
  Array<int,1> activeSpecies(activeSpeciesNames.size());
  for (int i=0; i<activeSpecies.size(); i++)
    activeSpecies(i) = PathData.Path.SpeciesNum(activeSpeciesNames(i));
  SetActiveSpecies (activeSpecies);

//   // Move all particles at the same time.
//   int totalNum = 0;
//   for (int si=0; si<activeSpecies.size(); si++)
//     totalNum += PathData.Path.Species(activeSpecies(si)).NumParticles;
//   SetNumParticlesToMove(totalNum);

  // Construct action list
  SwapStage.Actions.push_back(&PathData.Actions.ShortRange);
  SwapStage.Actions.push_back(&PathData.Actions.Kinetic);
  if (PathData.Path.LongRange) 
    if (PathData.Actions.UseRPA)
      SwapStage.Actions.push_back(&PathData.Actions.LongRangeRPA);
    else
      SwapStage.Actions.push_back(&PathData.Actions.LongRange);

  for (int i=0; i<activeSpecies.size(); i++) {
    int speciesNum = activeSpecies(i);
    if ((PathData.Actions.NodalActions(speciesNum)!=NULL)) {
      cerr << "SwapMove adding fermion node action for species " 
	   << activeSpeciesNames(i) << endl;
      SwapStage.Actions.push_back
	(PathData.Actions.NodalActions(speciesNum));
    }
  }
    
  // Now construct stage list
  Stages.push_back(&SwapStage);

  ActiveParticles.resize(NumParticlesToMove);
}

void
SwapMoveClass::MakeMove ()
{

  if (PathData.Path.Random.Common()>0.5)
    return;
  // Move the Join out of the way.
  PathData.MoveJoin (PathData.Path.NumTimeSlices()-1);
  Array<int,1> permVec;
  PathData.Path.TotalPermutation(permVec);
  PathData.Path.Communicator.Broadcast(0,permVec);
  vector<int> ptclList;
  int numLeft=0;
  for (int i=0; i<ActiveSpecies.size(); i++) {
    SpeciesClass &species = PathData.Path.Species(ActiveSpecies(i));
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      if (permVec(ptcl)==ptcl){
	ptclList.push_back(ptcl);
	numLeft++;
      }
    }
  }
  
  // First, choose particle to move
  for (int i=0; i<NumParticlesToMove; i++) {
    int index = PathData.Path.Random.CommonInt(numLeft);
    vector<int>::iterator iter = ptclList.begin();
    ActiveParticles(i) = ptclList[index];
    for (int j=0; j<index; j++)
      iter++;
    ptclList.erase(iter);
    numLeft--;
  }

//   int ptclIndex = 0;
//   for (int si=0; si < ActiveSpecies.size(); si++) {
//     SpeciesClass &species = PathData.Path.Species(ActiveSpecies(si));
//     for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
//       ActiveParticles(ptclIndex) = ptcl;
//       ptclIndex++;
//     }
//   }
  // Next, set timeslices
  Slice1 = 0;
  Slice2 = PathData.Path.NumTimeSlices()-1;

  // Now call MultiStageClass' MakeMove
  MultiStageClass::MakeMove();
}
