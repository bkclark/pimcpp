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

#ifndef STAGE_H
#define STAGE_H

// #include "../MPI/Communication.h"
#include "../Actions/ActionBase.h"
#include <list>
#include "../Observables/ObservableVar.h"

// Some moves are built out of stages.  This class is the parent class
// for a stage.  If it is the first stage it will typically need to
// set the slices adn the activeparticles for the other stages to
// perform correctly. Each stage has a list of actions that lets it
// calculate its total action by summing the value of each action in
// its list.
class StageClass
{
protected:
  PathDataClass &PathData;
  PathClass &Path;
  IOSectionClass OutSection;
  ObservableDouble AcceptRatioVar;
public:
  int NumAccepted, NumAttempted;
  double TimeSpent;
  int BisectionLevel;
  list<ActionBaseClass*> Actions;

  // Stores the transition probability that you have made this move. 
  // This is not implemented in all stages.
  double NewSample;
  double OldSample;
  double AcceptProb;
  double OldAcceptProb;

  // Flag for faster simulation for free particles
  bool IsFree;

  // The first stage will set the slices and activeParticles
  // This returns transition probability ratio T(new->old)/T(old->new)
  virtual double Sample (int &slice1, int &slice2, Array<int,1> &activeParticles) = 0;
  virtual void Read (IOSectionClass &in);
  virtual void Accept();
  virtual void Reject();
  virtual void WriteRatio();
  double StageAction(int startSlice, int endSlice, const Array<int,1> &changedParticles);
  inline double GlobalStageAction (const Array<int,1> &changeParticles);
  inline double AcceptRatio () {
    return (double)NumAccepted / (double) NumAttempted; 
  }

  virtual bool Attempt(int &slice1, int &slice2, Array<int,1> &activeParticles, double &prevActionChange) = 0;


  StageClass(PathDataClass &pathData, IOSectionClass outSection) :
    PathData(pathData), Path(pathData.Path), NumAccepted(0), NumAttempted(0), BisectionLevel(0), OutSection(outSection), AcceptRatioVar("AcceptRatio",OutSection,pathData.Path.Communicator)
  {
    if (PathData.Path.Communicator.MyProc()==0)
      OutSection.NewSection("Stage");
    IsFree = false;
  }
};

class CommonStageClass : public StageClass
{
public:
  virtual bool Attempt(int &slice1, int &slice2, Array<int,1> &activeParticles, double &prevActionChange);
  CommonStageClass(PathDataClass &pathData,IOSectionClass &outSection) :
    StageClass(pathData,outSection)
  {
    // Do nothing for now
  }
};

class LocalStageClass : public StageClass
{
public:
  double UAction; // jg just for testing
  bool Attempt(int &slice1, int &slice2, Array<int,1> &activeParticles, double &prevActionChange);
  LocalStageClass(PathDataClass &pathData, IOSectionClass &outSection) :
    StageClass(pathData, outSection)
  {
    UAction = 0.0;
  }
};

#endif
