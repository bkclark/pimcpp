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

#include "CorrelatedBisectionBlock.h"
#include "StructureRejectStage.h"
#include "NoPermuteStage.h"

void CorrelatedBisectionBlockClass::Read(IOSectionClass &in)
{

  string permuteType, speciesName;
  //  StageClass *permuteStage;
  assert (in.ReadVar ("NumLevels", NumLevels));
  assert (in.ReadVar ("Species", speciesName));
  assert (in.ReadVar ("StepsPerBlock", StepsPerBlock));

  BisectionStages.resize(NumLevels);
  for (int level=0;level<NumLevels;level++)
    BisectionStages(NumLevels-level-1)=new BisectionStageClass(PathData,level,IOSection);
  
  SpeciesNum = PathData.Path.SpeciesNum (speciesName);
  
  HaveRefslice = 
    ((PathData.Path.Species(SpeciesNum).GetParticleType() == FERMION) &&
     (PathData.Actions.NodalActions(SpeciesNum) != NULL) &&
     (!PathData.Actions.NodalActions(SpeciesNum)->IsGroundState()));
  /// Set up permutation
  assert (in.ReadVar ("PermuteType", permuteType));
  if (permuteType == "TABLE") 
    PermuteStage = 
      new TablePermuteStageClass(PathData, SpeciesNum, NumLevels,IOSection);
  else if (permuteType=="COUPLE")
    PermuteStage = 
      new CoupledPermuteStageClass(PathData,SpeciesNum,NumLevels,IOSection);
  else if (permuteType == "NONE") 
    PermuteStage = 
      new NoPermuteStageClass(PathData, SpeciesNum, NumLevels, IOSection);
  else {
    cerr << "Unrecognized PermuteType, """ << permuteType << """\n";
    exit(EXIT_FAILURE);
  }
  PermuteStage->Read (in);
  
  BosonActions.push_back(&PathData.Actions.Kinetic);
  BosonActions.push_back(&PathData.Actions.ShortRange);
  if (PathData.Path.LongRange){
    if (PathData.Actions.UseRPA)
      BosonActions.push_back(&PathData.Actions.LongRangeRPA);
    else
      BosonActions.push_back(&PathData.Actions.LongRange);
  }
  for (int species=0;species<PathData.NumSpecies();species++)
    if ((PathData.Actions.NodalActions(species)!=NULL)) {
      cerr << "Adding fermion node action for species " 
	   << speciesName << endl;
      NodalActions.push_back(PathData.Actions.NodalActions(species));
    }

  /// Initialize values for total action
  double actionA = 0.0; double actionB = 0.0;
  int s1 = 0; 
  int s2 = PathData.Path.NumTimeSlices()-1;
  Array<int,1> allParticles(PathData.Path.NumParticles());
  for (int i=0; i<allParticles.size(); i++)
    allParticles(i) = i;
  list<ActionBaseClass*>::iterator iter;
  PathData.Path.SetIonConfig(0);
  for (iter=BosonActions.begin(); iter!=BosonActions.end(); iter++)
    actionA += (*iter)->Action(s1, s2, allParticles, 0);
  for (iter=NodalActions.begin(); iter!=NodalActions.end(); iter++)
    actionA += (*iter)->Action(s1, s2, allParticles, 0);
  PathData.Path.SetIonConfig(1);
  for (iter=BosonActions.begin(); iter!=BosonActions.end(); iter++)
    actionB += (*iter)->Action(s1, s2, allParticles, 0);
  for (iter=NodalActions.begin(); iter!=NodalActions.end(); iter++)
    actionB += (*iter)->Action(s1, s2, allParticles, 0);
  PathData.Actions.TotalA = PathData.Path.Communicator.AllSum(actionA);
  PathData.Actions.TotalB = PathData.Path.Communicator.AllSum(actionB);
}


void CorrelatedBisectionBlockClass::ChooseTimeSlices()
{
  //  if (PathData.Path.Communicator.MyProc()==0)
    //    cerr<<"Choosing time slices"<<endl;
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
  //  if (PathData.Path.Communicator.MyProc()==0)
    //    cerr<<"Ending Choosing time slices"<<endl;
}

void CorrelatedBisectionBlockClass::MakeMove()
{
  perr << "In CorrelatedBisectionBlockClass::MakeMove().\n";
  ChooseTimeSlices();
  PathData.MoveJoin(Slice2);
  PermuteStage->InitBlock(Slice1,Slice2);
  ActiveParticles.resize(1);
  Array<int,1> allParticles(PathData.Path.NumParticles());
  for (int i=0; i<allParticles.size(); i++)
    allParticles(i) = i;
  for (int step=0; step<StepsPerBlock; step++) {
    ActiveParticles(0) = -1;
    double logSampleProb=0.0;
    ///This is where we call the permtuations
    double prevActionChange=0.0;
    PermuteStage->Attempt(Slice1,Slice2,ActiveParticles,prevActionChange);
    for (int bIndex=0;bIndex<BisectionStages.size();bIndex++)
      logSampleProb += 
	log(BisectionStages(bIndex)->Sample(Slice1,Slice2,ActiveParticles));
    
    list<ActionBaseClass*>::iterator iter;
    PathData.Path.SetIonConfig(0);
    double bosonSAOld, bosonSANew, bosonSBOld, bosonSBNew;
    bosonSAOld = 0.0; bosonSANew = 0.0; bosonSBOld = 0.0; bosonSBNew =0.0;
    double nodalSAOld, nodalSANew, nodalSBOld, nodalSBNew;
    nodalSAOld = 0.0; nodalSANew = 0.0; nodalSBOld = 0.0; nodalSBNew =0.0;
    for (iter=BosonActions.begin();iter!=BosonActions.end();iter++) {
      SetMode(NEWMODE);
      bosonSANew += (*iter)->Action(Slice1,Slice2,ActiveParticles,0);
      SetMode(OLDMODE);
      bosonSAOld += (*iter)->Action(Slice1,Slice2,ActiveParticles,0);
    }
    for (iter=NodalActions.begin();iter!=NodalActions.end();iter++) {
      SetMode(NEWMODE);
      nodalSANew += (*iter)->Action(Slice1,Slice2,ActiveParticles,0);
      SetMode(OLDMODE);
      nodalSAOld += (*iter)->Action(Slice1,Slice2,ActiveParticles,0);
    }

    PathData.Path.SetIonConfig(1);
    for (iter=BosonActions.begin();iter!=BosonActions.end();iter++) {
      SetMode(NEWMODE);
      bosonSBNew += (*iter)->Action(Slice1,Slice2,ActiveParticles,0);
      SetMode(OLDMODE);
      bosonSBOld += (*iter)->Action(Slice1,Slice2,ActiveParticles,0);
    }
    for (iter=NodalActions.begin();iter!=NodalActions.end();iter++) {
      SetMode(NEWMODE);
      nodalSBNew += (*iter)->Action(Slice1,Slice2,ActiveParticles,0);
      SetMode(OLDMODE);
      nodalSBOld += (*iter)->Action(Slice1,Slice2,ActiveParticles,0);
    }

//     double localSBarOld=0.5*(bosonSAOld+bosonSBOld);
//     double localSBarNew=0.5*(bosonSANew+bosonSBNew);
    double localSBarOld=0.5*(bosonSAOld+bosonSBOld+nodalSAOld+nodalSBOld);
    double localSBarNew=0.5*(bosonSANew+bosonSBNew+nodalSANew+nodalSBNew);
    double logAcceptProb= logSampleProb - localSBarNew + localSBarOld;
    bool toAcceptLocal = logAcceptProb>=log(PathData.Path.Random.Local());

    /// HACK
    toAcceptLocal = true;
    ///////////////////////////////////////////
    /// New total calculation using updates ///
    ///////////////////////////////////////////
    double Achange = 
      toAcceptLocal ? (bosonSANew+nodalSANew - (bosonSAOld+nodalSAOld)) : 0.0;
    double Bchange = 
      toAcceptLocal ? (bosonSBNew+nodalSBNew - (bosonSBOld+nodalSBOld)) : 0.0;

    AllSumVecIn(0) = Achange;   AllSumVecIn(1) = Bchange;
    PathData.Path.Communicator.AllSum (AllSumVecIn, AllSumVecOut);
    Achange = AllSumVecOut(0);  Bchange = AllSumVecOut(1);

    double totalSANew   = PathData.Actions.TotalA + Achange;
    double totalSBNew   = PathData.Actions.TotalB + Bchange;
    double deltaSNew    = 0.5*(totalSANew - totalSBNew);
    double deltaSOld    = 0.5*(PathData.Actions.TotalA-PathData.Actions.TotalB);
    double totalSbarNew = 0.5*(totalSANew + totalSBNew);
    double totalSbarOld = 0.5*(PathData.Actions.TotalA+PathData.Actions.TotalB);

    /// HACK HACK HACK
    //    commonAcceptProb = cosh(deltaSNew)/cosh(deltaSOld);
 //    double commonAcceptProb = 
//       cosh(deltaSNew)/cosh(deltaSOld) * exp(totalSbarOld - totalSbarNew);
    logSampleProb = PathData.Path.Communicator.AllSum(logSampleProb);
    double logCoshNew, logCoshOld;
    if (deltaSNew > 100.0)
      logCoshNew = log(0.5) + deltaSNew;
    else if (deltaSNew < -100.0)
      logCoshNew = log(0.5) - deltaSNew;
    else
      logCoshNew = log (cosh(deltaSNew));

    if (deltaSOld > 100.0)
      logCoshOld = log(0.5) + deltaSOld;
    else if (deltaSOld < -100.0)
      logCoshOld = log(0.5) - deltaSOld;
    else
      logCoshOld = log (cosh(deltaSOld));


    double logCommonAcceptProb = 
      logCoshNew - logCoshOld + totalSbarOld-totalSbarNew+logSampleProb;

    //    cerr << "logCommonAcceptProb = " << logCommonAcceptProb << endl;

    bool toAcceptCommon = 
      logCommonAcceptProb>=log(PathData.Path.Random.Common());
    
    double wA, wB;
    if (toAcceptCommon) {      /// Accept
      /// HACK HACK HACK
      wA = exp(-deltaSNew)/(2.0*cosh(deltaSNew));
      wB = 1.0-wA;
//      wA = exp(-deltaSNew);
//      wB = exp(+deltaSNew);
      PathData.Actions.TotalA = totalSANew;
      PathData.Actions.TotalB = totalSBNew;
      if (toAcceptLocal)
	Accept();
      else
	Reject();
      perr << "Accepted. SA = " << totalSANew 
	   << "    SB = " << totalSBNew << endl;
    }
    else {                      /// Reject
      /// HACK HACK HACK
      wA = exp(-deltaSOld)/(2.0*cosh(deltaSOld));
      wB = 1.0-wA;
//      wA = exp(-deltaSNew);
//      wB = exp(+deltaSNew);
      Reject();
      perr << "Rejected. SA = " << PathData.Actions.TotalA 
	   << "    SB = " << PathData.Actions.TotalB << endl;

    }

    int s1 = 0;
    int s2 = PathData.Path.NumTimeSlices()-1;
    double energyA, energyB;
    energyA = energyB = 0.0;

//     PathData.Path.SetIonConfig(0);
//     for (iter=BosonActions.begin(); iter!=BosonActions.end(); iter++) 
//       energyA += (*iter)->d_dBeta(s1,s2,0);
//     for (iter=NodalActions.begin(); iter!=NodalActions.end(); iter++) 
//       energyA += (*iter)->d_dBeta(s1,s2,0);


//     PathData.Path.SetIonConfig(1);
//     for (iter=BosonActions.begin(); iter!=BosonActions.end(); iter++) 
//       energyB += (*iter)->d_dBeta(s1,s2,0);
//     for (iter=NodalActions.begin(); iter!=NodalActions.end(); iter++) 
//       energyB += (*iter)->d_dBeta(s1,s2,0);

//     AllSumVecIn(0) = energyA;   AllSumVecIn(1) = energyB;
//     PathData.Path.Communicator.AllSum(AllSumVecIn, AllSumVecOut);
//     energyA = AllSumVecOut(0) / PathData.Path.TotalNumSlices;
//     energyB = AllSumVecOut(1) / PathData.Path.TotalNumSlices;

    PathData.Path.WeightA.push_back(wA);
    PathData.Path.EnergyA.push_back(energyA);
    PathData.Path.WeightB.push_back(wB);
    PathData.Path.EnergyB.push_back(energyB);

    StepNum++;
  }

  if (((StepNum/StepsPerBlock) % 20) == 19) {
    double wASum   = 0.0;      double wBSum   = 0.0;
    double wAEASum   = 0.0;    double wBEBSum   = 0.0;

    for (int i=0; i<PathData.Path.WeightA.size(); i++) {
      wAEASum += (PathData.Path.EnergyA[i]*PathData.Path.WeightA[i]);
      wASum   += PathData.Path.WeightA[i];
      wBEBSum += (PathData.Path.EnergyB[i]*PathData.Path.WeightB[i]);
      wBSum   += PathData.Path.WeightB[i];
    }

    wAEAvar.Write(wAEASum);    wBEBvar.Write(wBEBSum);
    wAvar.Write(wASum);        wBvar.Write(wBSum);
    wAvar.Flush();

    PathData.Path.WeightA.clear();
    PathData.Path.EnergyA.clear();
    PathData.Path.WeightB.clear();
    PathData.Path.EnergyB.clear();
  }
}
