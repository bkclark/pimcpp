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

#include "DisplaceMoveFast.h"

double DisplaceFastStageClass::Sample (int &slice1, int &slice2,
				   Array<int,1> &activeParticles)
{
  /// Now, choose a random displacement 
  for (int ptclIndex=0; ptclIndex<activeParticles.size(); ptclIndex++) {
    int ptcl = activeParticles(ptclIndex);
    dVec disp;
    ///    PathData.Path.Random.CommonGaussianVec (Sigma, disp);
    disp(0)=PathData.Path.Random.Common()-0.5;disp(1)=PathData.Path.Random.Common()-0.5;disp(2)=PathData.Path.Random.Common()-0.5;
    disp=disp*Sigma;
    

    // Actually displace the path
    SetMode(NEWMODE);
    for (int slice=0; slice<PathData.NumTimeSlices(); slice++)
      PathData.Path(slice, ptcl) = PathData.Path(slice, ptcl) + disp;
  }

  // And return sample probability ratio
  return 1.0;
}

void
DisplaceFastMoveClass::Read (IOSectionClass &in)
{
  CurrentPtcl=0;
  assert(in.ReadVar ("Sigma", Sigma));
  DisplaceFastStage.Sigma = Sigma;
  DisplaceFastStage.BisectionLevel=3;
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
  if (PathData.Path.OrderN){
    DisplaceFastStage.Actions.push_back(&PathData.Actions.ShortRangeOn);
    EmptyStage.Actions.push_back(&PathData.Actions.ShortRangeOn);
  }
  else
    DisplaceFastStage.Actions.push_back(&PathData.Actions.ShortRange);
  if (PathData.Path.LongRange) 
    if (PathData.Actions.UseRPA)
      DisplaceFastStage.Actions.push_back(&PathData.Actions.LongRangeRPA);
    else if (PathData.Path.DavidLongRange){
      EmptyStage.Actions.push_back(&PathData.Actions.DavidLongRange);
    }
    else
      DisplaceFastStage.Actions.push_back(&PathData.Actions.LongRange);

  for (int i=0; i<activeSpecies.size(); i++) {
    int speciesNum = activeSpecies(i);
    if ((PathData.Actions.NodalActions(speciesNum)!=NULL)) {
      cerr << "DisplaceFastMove adding fermion node action for species " 
	   << activeSpeciesNames(i) << endl;
      DisplaceFastStage.Actions.push_back
	(PathData.Actions.NodalActions(speciesNum));
    }
  }
    
  // Now construct stage list
  Stages.push_back(&DisplaceFastStage);
  Stages.push_back(&EmptyStage);
  ActiveParticles.resize(NumParticlesToMove);
}

void 
DisplaceFastMoveClass::MakeMove()
{
  // Next, set timeslices
  Slice1 = 0;
  Slice2 = PathData.Path.NumTimeSlices()-1;

  for (int i=0;i<PathData.Path.NumParticles();i++){
    ActiveParticles(0)=i;
    // Now call MultiStageClass' MakeMove
    if (PathData.Path.Random.Common()>0.1)
      MultiStageClass::MakeMove();
  }


}

// void
// DisplaceFastMoveClass::MakeMove ()
// {
//   // Move the Join out of the way.
//   PathData.MoveJoin (PathData.Path.NumTimeSlices()-1);
//   Array<int,1> permVec;
//   PathData.Path.TotalPermutation(permVec);
//   PathData.Path.Communicator.Broadcast(0,permVec);
//   vector<int> ptclList;
//   int numLeft=0;
//   for (int i=0; i<ActiveSpecies.size(); i++) {
//     SpeciesClass &species = PathData.Path.Species(ActiveSpecies(i));
//     for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
//       if (permVec(ptcl)==ptcl){
// 	ptclList.push_back(ptcl);
// 	numLeft++;
//       }
//     }
//   }
//   if (numLeft==0) //if there are no non-permuted particles, abort move
//     return;
//   // First, choose particle to move
//   for (int i=0; i<NumParticlesToMove; i++) {
//     int index = PathData.Path.Random.CommonInt(numLeft);
//     //    cerr<<"Index is "<<index<<endl;
//     //    cerr<<"numLeft is  "<<numLeft<<endl;
//     vector<int>::iterator iter = ptclList.begin();
//     ActiveParticles(i) = ptclList[index];
//     for (int j=0; j<index; j++)
//       iter++;
//     ptclList.erase(iter);
//     numLeft--;
//   }

// //   int ptclIndex = 0;
// //   for (int si=0; si < ActiveSpecies.size(); si++) {
// //     SpeciesClass &species = PathData.Path.Species(ActiveSpecies(si));
// //     for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
// //       ActiveParticles(ptclIndex) = ptcl;
// //       ptclIndex++;
// //     }
// //   }
//   // Next, set timeslices
//   Slice1 = 0;
//   Slice2 = PathData.Path.NumTimeSlices()-1;

//   // Now call MultiStageClass' MakeMove
//   MultiStageClass::MakeMove();
// }
