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

#ifndef COUPLING_STAGE_CLASS_H
#define COUPLING_STAGE_CLASS_H

#include "MultiStage.h"
#include "../Observables/ObservableVar.h"


class CouplingStageClass : public LocalStageClass
{
private:
  void ChooseTimeSlices(int &slice1,int &slice2);
  IOSectionClass OutSection;
  Array<int,1> AcceptRatio;
  ObservableVecDouble1 AcceptRatioVar;
public:
  void WriteRatio();
  void Accept();
  void Reject();
  void Read(IOSectionClass  &in);
  double Sample(int &slice1,int &slice2, 
		Array<int,1> &activeParticles);
  CouplingStageClass(PathDataClass &pathData, int numLevels,
		IOSectionClass &outSection) : 
    LocalStageClass(pathData,outSection),
    OutSection(outSection),
    AcceptRatioVar("Acceptance Ratio",OutSection,PathData.Path.Communicator)

  { 
    AcceptRatio.resize(2);
    AcceptRatio=0;

    //do nothing for now

  }
};

#endif

