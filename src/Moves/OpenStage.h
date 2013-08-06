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

#ifndef OPEN_STAGE_CLASS_H
#define OPEN_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableBase.h"
#include "PermuteStage.h"

class OpenStageClass : public PermuteStageClass
{
private:
  IOSectionClass OutSection;
  Array<int,1> AcceptRatio;
  ObservableVecDouble1 AcceptRatioVar;
  int EndAttempts;
public:
  void WriteRatio();
  void Accept();
  void Reject();
  void Read(IOSectionClass  &in);
  bool Attempt (int &slice1, int &slice2, 
		Array<int,1> &activeParticles, double &prevActionChange);
  dVec RandomBoxLocation();
  //  int NumLevels;
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  OpenStageClass(PathDataClass &pathData, int speciesNum,
		 int numLevels,
		 IOSectionClass &outSection) : 
    //    LocalStageClass(pathData,outSection),
    PermuteStageClass(pathData, speciesNum, numLevels,outSection),
    OutSection(outSection),
    AcceptRatioVar("Acceptance Ratio",OutSection,PathData.Path.Communicator)
  { 
    EndAttempts=0;
    AcceptRatio.resize(2);
    AcceptRatio=0;
    ///    BisectionLevel=level;
    //do nothing for now

  }
};

#endif

