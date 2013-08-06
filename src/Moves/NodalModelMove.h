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

#ifndef NODALMODEL_MOVE_H
#define NODALMODEL_MOVE_H

#include "../PathDataClass.h"
#include "MultiStage.h"


class NodalModelStageClass : public CommonStageClass
{
private:
  ObservableVecInt2 AcceptArrayVar;
  ObservableVecInt2 AttemptArrayVar;
public:
  Array<int,2> AcceptArray, AttemptArray;
  void Accept();
  void Reject();
  void WriteRatio();

  int oldModel;
  Array<int,1> activeSpecies;

  double Sample (int &slice1, int &slice2, Array <int,1> &activeParticles);
  NodalModelStageClass (PathDataClass &pathData,IOSectionClass &outSection) :
    CommonStageClass (pathData,outSection),
    AttemptArrayVar("AttemptRatios", OutSection, pathData.Path.Communicator),
    AcceptArrayVar("AcceptRatios", OutSection, pathData.Path.Communicator)
  {
    // Do nothing for now.
  }
};


class NodalModelMoveClass : public MultiStageClass
{
private:
  NodalModelStageClass NodalModelStage;
  ObservableInt NumAttemptedVar;
  ObservableInt NumAcceptedVar;
  int NumAttempted;
  Array<int,1> activeSpecies;
public:

  inline double AcceptanceRatio()
  {
    double accept_ratio=(double)(NumAccepted)/(double)NumAttempted;
    return accept_ratio;
  }
  void WriteRatio();
  void Read (IOSectionClass &in);
  // Actually attempts the move and accepts or rejects
  void MakeMove();
  NodalModelMoveClass (PathDataClass &pathData, IOSectionClass &outSection) :
    NumAttempted(0),
    MultiStageClass(pathData, outSection),
    NodalModelStage(pathData, outSection),
    NumAttemptedVar("NumAttempted", IOSection, pathData.Path.Communicator),
    NumAcceptedVar("NumAccepted", IOSection, pathData.Path.Communicator)
  {
    NumAccepted = NumAttempted = 0;
  }
};

#endif
