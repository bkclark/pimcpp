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
    Stages.push_back (newStage);
  }

  /// Add centroid displace stage
  /// HACK! Not adding Kinetic to save time.
  //if (myProc == 0)
  //  cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" Adding Kinetic Action"<<endl;
  //DisplaceStage.Actions.push_back(&PathData.Actions.Kinetic);
  for (int i=0;i<samplingActions.size();i++) {
    if (myProc == 0)
      cout<<PathData.Path.CloneStr<<" "<<moveName<<" "<<speciesName<<" Adding "<<(*PathData.Actions.GetAction(samplingActions(i))).GetName()<<" Action"<<endl;
    DisplaceStage.Actions.push_back(PathData.Actions.GetAction(samplingActions(i)));
  }
  Stages.push_back(&DisplaceStage);

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

  ActiveParticles.resize(1);
  NumAttempted++;
  int myPtcl = (PathData.Path.Random.LocalInt(PathData.Path.Species(SpeciesNum).NumParticles)) + PathData.Path.Species(SpeciesNum).FirstPtcl;
  ActiveParticles(0) = myPtcl;
  gettimeofday(&start, &tz);
  MultiStageClass::MakeMove();
  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);
}


void CentroidStageClass::Accept()
{
  NumAccepted++;
  NumAttempted++;
  CommonStageClass::Accept();
}


void CentroidStageClass::Reject()
{
  NumAttempted++;
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

  PathData.GetCentroids(newCentPos, activeParticles);
  for (int i=0; i<N; i++) {
    for (int d=0; d<NDIM; d++) {
      if(abs(float(oldCentPos(i)(d))-float(newCentPos(i)(d))) > 1e-8) {
        cout << oldCentPos(i)(d) << " " << newCentPos(i)(d) << endl;
        cout << oldCentPos << endl;
        cout << newCentPos << endl;
        abort();
      }
    }
  }

  return log(1.0);
}


bool CentroidStageClass::Attempt (int &slice1, int &slice2, Array<int,1> &activeParticles, double &prevActionChange)
{
  // HACK: We are assuming previous step was accepted by setting
  // prevActionChange = 0. This really only works for Boltzmannons.
  prevActionChange = 0;

  bool toAccept = CommonStageClass::Attempt(slice1, slice2, activeParticles, prevActionChange);
  return toAccept;
}

