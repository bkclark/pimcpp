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

#include "BisectionSphereBlock.h"
#include "NoPermuteStage.h"
#include "StructureRejectStage.h"

//void BisectionSphereBlockClass::WriteRatio()
//{
//  double AcceptanceRatio=(double)NumAccepted/(double)NumAttempted;
//  double SumAcceptance=PathData.Path.Communicator.Sum(AcceptanceRatio);
//  AcceptanceRatioVar.Write(SumAcceptance);
//  MultiStageClass::WriteRatio();


//}

void BisectionSphereBlockClass::Read(IOSectionClass &in)
{
  string permuteType, speciesName;
  StageClass *permuteStage;
  assert (in.ReadVar ("NumLevels", NumLevels));
  assert (in.ReadVar ("Species", speciesName));
  assert (in.ReadVar ("StepsPerBlock", StepsPerBlock));
  cerr << "speciesName = " << speciesName << ".\n";
  SpeciesNum = PathData.Path.SpeciesNum (speciesName);
  HaveRefslice = 
    ((PathData.Path.Species(SpeciesNum).GetParticleType() == FERMION) &&
     (PathData.Actions.NodalActions(SpeciesNum) != NULL) &&
     (!PathData.Actions.NodalActions(SpeciesNum)->IsGroundState()));

  /// Set up permutation
  assert (in.ReadVar ("PermuteType", permuteType));
  if (permuteType == "TABLE") 
    permuteStage = new TablePermuteStageClass(PathData, SpeciesNum, NumLevels,
					      IOSection);
  else if (permuteType=="COUPLE")
    permuteStage= new CoupledPermuteStageClass(PathData,SpeciesNum,NumLevels,
					       IOSection);
  else if (permuteType == "NONE") 
    permuteStage = new NoPermuteStageClass(PathData, SpeciesNum, NumLevels,
					   IOSection);
  else {
    cerr << "Unrecognized PermuteType, """ << permuteType << """\n";
    exit(EXIT_FAILURE);
  }
  permuteStage->Read (in);
  Stages.push_back (permuteStage);
  
  for (int level=NumLevels-1; level>=0; level--) {
    BisectionStageSphereClass *newStage = 
      new BisectionStageSphereClass (PathData, level, IOSection);
    newStage->Actions.push_back(&PathData.Actions.Kinetic);
    if (PathData.Path.OpenPaths && level==0)
      newStage->Actions.push_back(&PathData.Actions.OpenLoopImportance);
    if (PathData.Path.OpenPaths && level>0)
      newStage->Actions.push_back(&PathData.Actions.ShortRangeApproximate);
    else
      newStage->Actions.push_back(&PathData.Actions.ShortRange);
    if (level == 0) {
      ////ADDED JUST FOR THE MINUTE HACK! FOR THE hehp project
      ////      newStage->Actions.push_back(&PathData.Actions.StructureReject);
      ///If it's David's long range class then do this
      if (PathData.Path.DavidLongRange)
	newStage->Actions.push_back(&PathData.Actions.DavidLongRange);
      else if (PathData.Path.LongRange) {
	if (PathData.Actions.UseRPA)
	  newStage->Actions.push_back(&PathData.Actions.LongRangeRPA);
	else
	  newStage->Actions.push_back(&PathData.Actions.LongRange);
      }
      if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
	cerr << "Adding fermion node action for species " 
	     << speciesName << endl;
	
	newStage->Actions.push_back(PathData.Actions.NodalActions(SpeciesNum));
      }
//       for (int i=0; i<PathData.Actions.NodalActions.size(); i++)
// 	newStage->Actions.push_back(PathData.Actions.NodalActions(i));
    }

    newStage->BisectionLevel = level;
    Stages.push_back (newStage);
  }

  // Add the second stage of the permutation step
  Stages.push_back (permuteStage);

//   ///HACK! Addding a stage that will reject the move if the structure
//   //factor gets too large
//   StructureRejectStageClass* structureReject =
//     new StructureRejectStageClass(PathData,in);
//   structureReject->Read(in);
//   Stages.push_back(structureReject);
  // Add the nodal action stage, if necessary
  
}


void BisectionSphereBlockClass::ChooseTimeSlices()
{
  PathClass &Path = PathData.Path;
  int myProc = PathData.Path.Communicator.MyProc();
  // do something special to avoid moving reference slice
  if (HaveRefslice &&
      Path.SliceOwner(Path.GetRefSlice()) == myProc) {
    
    cerr << "Avoid reference slice.\n";
    int bSlices = 1<<NumLevels;
    int myStart, myEnd;
    Path.SliceRange (myProc, myStart, myEnd);
    // refSlice is relative to my first slice
    int refSlice = Path.GetRefSlice() - myStart;
    int numSlices = Path.NumTimeSlices();
    assert(bSlices*2<numSlices);
    if (refSlice < bSlices) {  
      int numStarts = numSlices - bSlices - refSlice;
      Slice1 = refSlice + Path.Random.LocalInt(numStarts);
      Slice2 = Slice1 + bSlices;
      if (Slice2>=numSlices){
	cerr << "(Slice 2, numSlices): (" << Slice2 << ", " 
	     << numSlices << ")" << " "<< bSlices << " "
	     << numStarts << " " << refSlice << endl;
      }
      assert (Slice2 < numSlices);
    }
    else if (refSlice > (numSlices-1-bSlices) && refSlice<numSlices) {
      int numStarts = refSlice - bSlices + 1;
      Slice1 = Path.Random.LocalInt(numStarts);
      Slice2 = Slice1+bSlices;
      assert (Slice2 <= refSlice);
    }
    else {
      int numBefore = refSlice - bSlices +1;
      int numAfter = numSlices -bSlices - refSlice;
      int numStarts = numBefore+numAfter;
      int start = Path.Random.LocalInt (numStarts);
      if (start < numBefore)
	Slice1 = start;
      else
	Slice1 = start-numBefore+refSlice;
      Slice2 = Slice1+bSlices;
      assert ((refSlice <= Slice1) || (refSlice >= Slice2));
      assert (Slice2 < numSlices);
    }
  }
  // Bryan should fill this part in.
  //  else if (PathData.Path.OpenPtcl) {
    

  //  }
  else {
    int sliceSep = 1<<NumLevels;
    assert (sliceSep < PathData.Path.NumTimeSlices());
    int numLeft = PathData.Path.NumTimeSlices()-sliceSep;
    Slice1 = PathData.Path.Random.LocalInt (numLeft);
    Slice2 = Slice1+sliceSep;
  }
}

void BisectionSphereBlockClass::MakeMove()
{
  //  perr << "BisectionSphereBlock MakeMove.\n";
  //  cerr<<"Starting my bisection block"<<endl;
  ChooseTimeSlices();
  //  cerr<<"Choosing Time slices"<<endl;
  PathData.MoveJoin(Slice2);
  //  cerr<<"Moving Join"<<endl;
  //  sleep(10);
  
  ActiveParticles.resize(1);
  for (int step=0; step<StepsPerBlock; step++) {
    //    cerr<<"Step number "<<step<<endl;
    //    sleep(10);
    ActiveParticles(0)=-1;
    MultiStageClass::MakeMove();
    //    cerr<<"Step number "<<step<<" done"<<endl;
    //    sleep(10);
  }
  //  cerr<<"Ending my bisection block"<<endl;
  //  sleep(10);
  
}
