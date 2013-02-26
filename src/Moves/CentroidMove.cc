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

#include "CentroidMove.h"
#include "NoPermuteStage.h"
#include "sys/time.h"


void CentroidMoveClass::WriteRatio()
{
  MultiStageClass::WriteRatio();
}


void CentroidMoveClass::Read(IOSectionClass &in)
{
  int myProc = PathData.Path.Communicator.MyProc();
  string moveName = "CentroidMove";
  string permuteType, speciesName;
  assert (in.ReadVar ("NumLevels", NumLevels));
  int LowestLevel = 0;
  assert (in.ReadVar ("Species", speciesName));
  SpeciesNum = PathData.Path.SpeciesNum (speciesName);

  /// Set up permutation
  // Right now it only makes sense to have no permutation
  //assert (in.ReadVar ("PermuteType", permuteType));
  PermuteStage = new NoPermuteStageClass(PathData, SpeciesNum, NumLevels, IOSection);
  PermuteStage->Read (in);
  Stages.push_back (PermuteStage);

  /// Sampling actions
  Array<string,1> samplingActions;
  samplingActions.resize(0);
  if(!in.ReadVar("SamplingActions",samplingActions))
    cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" WARNING: No sampling actions found! Treating as free particles."<<endl;

  /// Add Bisection stages
  for (int level=NumLevels-1; level>=LowestLevel; level--) {
    BisectionStageClass *newStage = new BisectionStageClass (PathData, level, IOSection);
    newStage->TotalLevels = NumLevels;
    newStage->BisectionLevel = level;
    newStage->UseCorrelatedSampling = false;
    if (myProc == 0)
      cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding Kinetic Action"<<endl;
    newStage->Actions.push_back(&PathData.Actions.Kinetic);
    if (level!=LowestLevel) {
      Array<string,1> higherLevelActions;
      if (in.ReadVar("HigherLevelActions",higherLevelActions)){
        for (int i=0;i<higherLevelActions.size();i++) {
          if (myProc == 0)
            cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding "<<(*PathData.Actions.GetAction(higherLevelActions(i))).GetName()<<" Action"<<endl;
          newStage->Actions.push_back(PathData.Actions.GetAction(higherLevelActions(i)));
        }
      }
    } else if (level == LowestLevel) {
      for (int i=0;i<samplingActions.size();i++) {
        if (myProc == 0)
          cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" "<<level<<" Adding "<<(*PathData.Actions.GetAction(samplingActions(i))).GetName()<<" Action"<<endl;
        newStage -> Actions.push_back(PathData.Actions.GetAction(samplingActions(i)));
      }
    }
    Stages.push_back (newStage);
  }

  /// Add centroid displace stage
  for (int i=0;i<samplingActions.size();i++) {
    if (myProc == 0)
      cout<<PathData.Path.CloneStr<<" "<<moveName<<" Adding "<<(*PathData.Actions.GetAction(samplingActions(i))).GetName()<<" Action"<<endl;
    DisplaceStage.Actions.push_back(PathData.Actions.GetAction(samplingActions(i)));
  }
  Stages.push_back(&DisplaceStage);

  // Add the second stage of the permutation step
  /// HACK!!!  Pushing onto the stack twice causes the stage
  /// to be accepted twice, which causes swapping the forward and
  // reverse tables twice!
  Stages.push_back (PermuteStage);

}


void CentroidMoveClass::ChooseTimeSlices()
{
  PathClass &Path = PathData.Path;
  int myProc = PathData.Path.Communicator.MyProc();
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


void CentroidMoveClass::MakeMove()
{
  struct timeval start, end;
  struct timezone tz;

  //HACK!
  ChooseTimeSlices();
  PathData.MoveJoin(Slice2);

  ((PermuteStageClass*)PermuteStage)->InitBlock(Slice1,Slice2);
  ActiveParticles.resize(1);
  NumAttempted++;
  ActiveParticles(0)=-1;
  gettimeofday(&start, &tz);
  MultiStageClass::MakeMove();
  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);
}


void CentroidStageClass::Accept()
{
  CommonStageClass::Accept();
}


void CentroidStageClass::Reject()
{
  CommonStageClass::Reject();
}


double CentroidStageClass::Sample (int &slice1, int &slice2, Array<int,1> &activeParticles)
{
  int N = activeParticles.size();

  /// Get old centroids
  Array<TinyVector<double,NDIM>,1> oldCentPos(N);
  SetMode(OLDMODE);
  PathData.GetCentroids(oldCentPos, activeParticles);

  /// Get new centroids
  Array<TinyVector<double,NDIM>,1> newCentPos(N);
  SetMode(NEWMODE);
  PathData.GetCentroids(newCentPos, activeParticles);

  // Actually displace the path
  SetMode(NEWMODE);
  for (int i=0; i<N; i++) {
    int ptcl = activeParticles(i);
    dVec disp;
    disp = newCentPos(i) - oldCentPos(i);
    Path.PutInBox(disp); // HACK!!! Do I need this?
    for (int slice=0; slice < PathData.Path.NumTimeSlices(); slice++)
      PathData.Path(slice, ptcl) = PathData.Path(slice, ptcl) - disp;
  }

  // And return log sample probability ratio
  return log(1.0);
}


