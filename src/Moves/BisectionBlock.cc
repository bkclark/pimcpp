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

#include "BisectionBlock.h"
#include "EmptyStage.h"
#include "StructureRejectStage.h"
#include "CouplingStage.h"
//#include "WormPermuteStage.h"
#include "OpenStage.h"
#include "NoPermuteStage.h"
#include "sys/time.h"


void BisectionBlockClass::Read(IOSectionClass &in)
{
  int myProc = PathData.Path.Communicator.MyProc();
  string moveName = "BisectionBlock";
  bool useCorrelatedSampling;
  if (!in.ReadVar("UseCorrelatedSampling",useCorrelatedSampling))
    useCorrelatedSampling=false;
  if (useCorrelatedSampling)
    cout<<"Using correlated sampling"<<endl;
  string permuteType, speciesName;
  assert (in.ReadVar ("NumLevels", NumLevels));
  LowestLevel = 0;
  if (!in.ReadVar ("LowestLevel", LowestLevel))
    LowestLevel = 0;
  assert (LowestLevel < NumLevels);
  assert (in.ReadVar ("Species", speciesName));
  assert (in.ReadVar ("StepsPerBlock", StepsPerBlock));
  SpeciesNum = PathData.Path.SpeciesNum (speciesName);
  if (PathData.Path.Species(SpeciesNum).GetParticleType() == FERMION)
    HaveRefslice=true;
  HaveRefslice = ((PathData.Path.Species(SpeciesNum).GetParticleType() == FERMION) &&
                  (PathData.Actions.NodalActions(SpeciesNum) != NULL) &&
                  (!PathData.Actions.NodalActions(SpeciesNum)->IsGroundState()));
  /// Set up permutation
  assert (in.ReadVar ("PermuteType", permuteType));
  if (permuteType == "TABLE")
    PermuteStage = new TablePermuteStageClass(PathData, SpeciesNum, NumLevels, IOSection);
  else if (permuteType=="COUPLE")
    PermuteStage = new CoupledPermuteStageClass(PathData, SpeciesNum, NumLevels, IOSection);
  else if (permuteType=="WORMMOVE")
    PermuteStage = new OpenStageClass(PathData, SpeciesNum, NumLevels, IOSection);
  else if (permuteType == "NONE")
    PermuteStage = new NoPermuteStageClass(PathData, SpeciesNum, NumLevels, IOSection);
  else {
    cerr << "Unrecognized PermuteType, """ << permuteType << """\n";
    exit(EXIT_FAILURE);
  }
  PermuteStage->Read (in);
  Stages.push_back (PermuteStage);

  // IsFree flag for faster runtimes for free particles
  bool IsFree;
  if (!in.ReadVar ("IsFree",IsFree))
    IsFree = false;
  if (IsFree)
    cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" Free particle sampling " << IsFree << endl;

  for (int level=NumLevels-1; level>=LowestLevel; level--) {
    BisectionStageClass *newStage = new BisectionStageClass (PathData, level, IOSection);
    newStage->TotalLevels = NumLevels;
    newStage->BisectionLevel = level;
    newStage->UseCorrelatedSampling = useCorrelatedSampling;
    newStage->IsFree = IsFree;
    if (!IsFree) {
      if (myProc == 0)
        cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding Kinetic Action"<<endl;
      newStage->Actions.push_back(&PathData.Actions.Kinetic);
    }
    if (level!=LowestLevel){
      Array<string,1> higherLevelActions;
      if (in.ReadVar("HigherLevelActions",higherLevelActions)){
        for (int i=0;i<higherLevelActions.size();i++) {
          if (myProc == 0)
            cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding "<<(*PathData.Actions.GetAction(higherLevelActions(i))).GetName()<<" Action"<<endl;
          newStage->Actions.push_back(PathData.Actions.GetAction(higherLevelActions(i)));
        }
      }
    }
    else if (level == LowestLevel) {
      Array<string,1> samplingActions;
      if(in.ReadVar("SamplingActions",samplingActions)) {
        for (int i=0;i<samplingActions.size();i++) {
          if (myProc == 0)
            cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding "<<(*PathData.Actions.GetAction(samplingActions(i))).GetName()<<" Action"<<endl;
          newStage -> Actions.push_back(PathData.Actions.GetAction(samplingActions(i)));
        }
      } else {
        if (myProc == 0)
          cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" WARNING: No sampling actions found! Treating as free particles."<<endl;
      }
      if ((PathData.Actions.NodalActions(SpeciesNum)!=NULL)) {
        if (myProc == 0)
          cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding Node Action"<<endl;
        newStage->Actions.push_back(PathData.Actions.NodalActions(SpeciesNum));
      }
    }
    //HACK!!! These used to be only pushed on at the lowest level
    //    if (PathData.Actions.UseNonlocal) {
    //      cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding Nonlocal Action"<<endl;
    //      newStage->Actions.push_back(&PathData.Actions.Nonlocal);
    //    }
    Stages.push_back (newStage);
  }
  // Add the second stage of the permutation step
  /// HACK!!!  Pushing onto the stack twice causes the stage
  /// to be accepted twice, which causes swapping the forward and
  // reverse tables twice!
  Stages.push_back (PermuteStage);

}


void BisectionBlockClass::ChooseTimeSlices()
{
  //  if (PathData.Path.Communicator.MyProc()==0)
    //    cerr<<"Choosing time slices"<<endl;
  PathClass &Path = PathData.Path;
  int myProc = PathData.Path.Communicator.MyProc();
  // do something special to avoid moving reference slice
  if (HaveRefslice &&
      Path.SliceOwner(Path.GetRefSlice()) == myProc) {
    int bSlices = 1<<NumLevels;
    int myStart, myEnd;
    Path.SliceRange (myProc, myStart, myEnd);
    // refSlice is relative to my first slice
    int refSlice = Path.GetRefSlice() - myStart;
    int numSlices = Path.NumTimeSlices();
    //cerr << "bSlices: " << bSlices << " numSlices:" << numSlices << " myProc:" << myProc << endl;
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
    if (numLeft < 0) {
      cerr << "Not enough slices to bisect with " << NumLevels 
	   << " levels in BisectionBlock.\n";
      abort();
    }
    Slice1 = PathData.Path.Random.LocalInt (numLeft);
    Slice2 = Slice1+sliceSep;
  }
  //  if (PathData.Path.Communicator.MyProc()==0)
  //    cerr<<"Ending Choosing time slices"<<endl;
  ////  cerr<<"Slices: "<<Slice1<<" "<<Slice2<<endl;
}


void BisectionBlockClass::MakeMove()
{
  struct timeval start, end;
  struct timezone tz;

  //HACK!
  ChooseTimeSlices();
  PathData.MoveJoin(Slice2);

  ((PermuteStageClass*)PermuteStage)->InitBlock(Slice1,Slice2);
  ActiveParticles.resize(1);
  for (int step=0; step<StepsPerBlock; step++) {
    NumAttempted++;
    ActiveParticles(0)=-1;
    gettimeofday(&start, &tz);
    MultiStageClass::MakeMove();
    gettimeofday(&end,   &tz);
    TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);
  }

  if (LowestLevel != 0)
    MakeStraightPaths();

}


void BisectionBlockClass::WriteRatio()
{
  //  PrintTimeSpent();
  MultiStageClass::WriteRatio();

}


void BisectionBlockClass::PrintTimeSpent()
{
  
  list<StageClass*>::iterator stageIter=Stages.begin();
  while (stageIter!=Stages.end()){
    //    cerr<<"LEVEL A IS "<<(*stageIter)->TimeSpent<<endl;
    cout<<"START"<<endl;
    for (list<ActionBaseClass*>::iterator actionIter=(*stageIter)->Actions.begin();actionIter!=(*stageIter)->Actions.end();actionIter++){
      cout<<"    Action value was "<<(*actionIter)->TimeSpent<<" "<<(*actionIter)->GetName()<<endl;
    }
    cout<<"Time spent is "<<TimeSpent<<endl;
    cout<<"END"<<endl;
    stageIter++;
  }
}


void BisectionBlockClass::MakeStraightPaths()
{
  PathClass &Path = PathData.Path;
  SetMode(NEWMODE);
  int skip = 1<<LowestLevel;
  int first = Path.Species(SpeciesNum).FirstPtcl;
  int last = Path.Species(SpeciesNum).LastPtcl;
  double inc = 1.0/(double)skip;
  for (int slice=Slice1; slice < Slice2; slice += skip) 
    for (int ptcl=first; ptcl<=last; ptcl++) {
      dVec delta = Path.Velocity(slice, slice+skip, ptcl);
      double frac = inc;
      for (int s=1; s<skip; s++) {
	Path(slice+s,ptcl) = Path(slice,ptcl) + frac*delta;
	frac += inc;
      }
    }
  Array<int,1> ptcls(last-first+1);
  for (int ptcl=first; ptcl<=last; ptcl++)
    ptcls(ptcl-first) = ptcl;
  PathData.AcceptMove(Slice1, Slice2, ptcls);
}
