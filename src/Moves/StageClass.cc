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

#include "StageClass.h"
#include "sys/time.h"

void StageClass::Read(IOSectionClass &in)
{
  ///Do nothing for now
}


void StageClass::WriteRatio()
{
  AcceptRatioVar.Write((double)NumAccepted/(double)NumAttempted);
  AcceptRatioVar.Flush();
}


//BUG: DOES NOT HAVE CORRECT SLICES!
void StageClass::Accept()
{
  NumAccepted++;
  NumAttempted++;
  for (list<ActionBaseClass*>::iterator actionIter = Actions.begin(); actionIter != Actions.end(); actionIter++)
    (*actionIter) -> AcceptCopy(0,0);
}


//BUG: DOES NOT HAVE CORRECT SLICES
void StageClass::Reject()
{
  NumAttempted++;
  for (list<ActionBaseClass*>::iterator actionIter = Actions.begin(); actionIter != Actions.end(); actionIter++)
    (*actionIter) -> RejectCopy(0,0);
}


bool LocalStageClass::Attempt(int &slice1, int &slice2, Array<int,1> &activeParticles, double &prevActionChange)
{
  struct timeval start, end;
  struct timezone tz;

  SetMode (NEWMODE);
  double logSampleRatio = Sample(slice1,slice2,activeParticles);

  bool toAccept;
  if (!IsFree) {
    SetMode (OLDMODE);
    gettimeofday(&start, &tz);
    double oldAction = StageAction(slice1,slice2,activeParticles);

    SetMode(NEWMODE);
    double newAction = StageAction(slice1,slice2,activeParticles);
    double currActionChange = newAction - oldAction;
    double logAcceptProb = logSampleRatio - currActionChange + prevActionChange;
    toAccept = logAcceptProb >= log(PathData.Path.Random.Local()); // Accept condition

    //cout << "Local Staging: " << toAccept << " " << BisectionLevel << " " << slice1 << " " << slice2 << " " << PathData.Path.GetRefSlice() << " " << PathData.Path.SliceOwner(PathData.Path.GetRefSlice()) << " " << PathData.Path.Communicator.MyProc() << " " << oldAction << " " << newAction << " " << logSampleRatio << " " << currActionChange << " " << prevActionChange << endl;
    if (abs(newAction) > 1e50 || abs(oldAction) > 1e50) {
      if (toAccept) {
        if (abs(newAction) > 1e50 && abs(oldAction) < 1e50) {
          cerr << PathData.Path.CloneStr <<" Broken Local Staging (new): " << BisectionLevel << " " << slice1 << " " << slice2 << " " << PathData.Path.GetRefSlice() << " " << PathData.Path.SliceOwner(PathData.Path.GetRefSlice()) << " " << oldAction << " " << newAction << " " << logSampleRatio << " " << currActionChange << " " << prevActionChange << endl;
          toAccept = 0;
          assert(1==2);
        } else if (abs(oldAction) > 1e50 && abs(newAction) < 1e50) {
          cerr << PathData.Path.CloneStr << " Broken Local Staging (old): " << BisectionLevel << " " << slice1 << " " << slice2 << " " << PathData.Path.GetRefSlice() << " " << PathData.Path.SliceOwner(PathData.Path.GetRefSlice()) << " " <<  oldAction << " " << newAction << " " << logSampleRatio << " " << currActionChange << " " << prevActionChange << endl;
        } else {
          cerr << PathData.Path.CloneStr <<" Broken Local Staging (both): " << BisectionLevel << " " << slice1 << " " << slice2 << " " << PathData.Path.GetRefSlice() << " " << PathData.Path.SliceOwner(PathData.Path.GetRefSlice()) << " " << oldAction << " " << newAction << " " << logSampleRatio << " " << currActionChange << " " << prevActionChange << endl;
          toAccept = 0;
        }
      }
    }

    prevActionChange = currActionChange;
  } else {
    toAccept = 0.0 >= log(PathData.Path.Random.Local()); // Accept condition
  }

  if (toAccept)
    NumAccepted++;
  NumAttempted++;

  gettimeofday(&end, &tz);
  TimeSpent += (double)(end.tv_sec-start.tv_sec) + 1.0e-6*(double)(end.tv_usec-start.tv_usec);

  return toAccept;
}


bool CommonStageClass::Attempt (int &slice1, int &slice2, Array<int,1> &activeParticles, double &prevActionChange)
{
  slice1 = 0;
  slice2 = PathData.NumTimeSlices()-1;
  assert (slice1 == 0);
  assert (slice2 == PathData.NumTimeSlices()-1);

  SetMode(OLDMODE);
  double oldAction = GlobalStageAction(activeParticles);

  SetMode(NEWMODE);
  double logSampleRatio = Sample(slice1,slice2,activeParticles);
  logSampleRatio = PathData.Path.Communicator.AllSum (logSampleRatio);

  SetMode(NEWMODE);
  double newAction = GlobalStageAction(activeParticles);
  double currActionChange = newAction - oldAction;
  double logAcceptProb = logSampleRatio - currActionChange + prevActionChange;
  bool toAccept = logAcceptProb >= log(PathData.Path.Random.Common()); /// Accept condition

  //cout << "Common Staging: " << toAccept << " " << BisectionLevel << " " << slice1 << " " << slice2 << " " << PathData.Path.GetRefSlice() << " " << PathData.Path.SliceOwner(PathData.Path.GetRefSlice()) << " " << PathData.Path.Communicator.MyProc() << " " << oldAction << " " << newAction << " " << logSampleRatio << " " << currActionChange << " " << prevActionChange << endl;
  if (abs(newAction) > 1e50 || abs(oldAction) > 1e50) {
    if (toAccept) {
      if (abs(newAction) > 1e50 && abs(oldAction) < 1e50) {
        cerr << PathData.Path.CloneStr <<" Broken Common Staging (new): " << BisectionLevel << " " << slice1 << " " << slice2 << " " << PathData.Path.GetRefSlice() << " " << PathData.Path.SliceOwner(PathData.Path.GetRefSlice()) << " " << oldAction << " " << newAction << " " << logSampleRatio << " " << currActionChange << " " << prevActionChange << endl;
        toAccept = 0;
        assert(1==2);
      } else if (abs(oldAction) > 1e50 && abs(newAction) < 1e50) {
        cerr << PathData.Path.CloneStr <<" Broken Common Staging (old): " << BisectionLevel << " " << slice1 << " " << slice2 << " " << PathData.Path.GetRefSlice() << " " << PathData.Path.SliceOwner(PathData.Path.GetRefSlice()) << " " << oldAction << " " << newAction << " " << logSampleRatio << " " << currActionChange << " " << prevActionChange << endl;
      } else {
        cerr << PathData.Path.CloneStr <<" Broken Common Staging (both): " << BisectionLevel << " " << slice1 << " " << slice2 << " " << PathData.Path.GetRefSlice() << " " << PathData.Path.SliceOwner(PathData.Path.GetRefSlice()) << " " << oldAction << " " << newAction << " " << logSampleRatio << " " << currActionChange << " " << prevActionChange << endl;
        toAccept = 0;
      }
    }
  }

  if (toAccept) {
    NumAccepted++;
  }
  NumAttempted++;
  prevActionChange = currActionChange;

  return toAccept;
}


double StageClass::StageAction (int startSlice, int endSlice, const Array<int,1> &changedParticles)
{
  double TotalAction = 0.0;
  list<ActionBaseClass*>::iterator actionIter = Actions.begin();
  while (actionIter != Actions.end()) {
    double TempAction = ((*actionIter) -> Action(startSlice, endSlice, changedParticles, BisectionLevel));
    TotalAction += TempAction;
    if (GetMode()==OLDMODE && (abs(TempAction) > 1e50))
      cerr << PathData.Path.CloneStr << " WARNING: " << (*actionIter) -> GetName() << " " << TempAction << " " << startSlice << " " << endSlice << " " << BisectionLevel << endl;
    actionIter++;
  }
  //cerr << "Local Stage Action " << TotalAction << endl;

  return TotalAction;
}


inline double StageClass::GlobalStageAction (const Array<int,1> &changedParticles)
{
  int slice1 = 0;
  int slice2 = PathData.Path.NumTimeSlices()-1;
  double localAction = 0.0;
  list<ActionBaseClass*>::iterator actionIter=Actions.begin();
  while (actionIter!=Actions.end()) {
    double TempAction = ((*actionIter) -> Action(slice1, slice2, changedParticles, BisectionLevel));
    localAction += TempAction;
    //cerr << (*actionIter) -> GetName() << " " << TempAction << " " << slice1 << " " << slice2 << " " << BisectionLevel <<  endl;
    actionIter++;
  }
  double globalAction = PathData.Path.Communicator.AllSum (localAction);
  //cerr << "Global Stage Action " << localAction << " " << globalAction << endl;

  return globalAction;
}

