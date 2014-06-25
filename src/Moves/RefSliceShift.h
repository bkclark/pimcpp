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

#ifndef REF_SLICE_SHIFT_H
#define REF_SLICE_SHIFT_H

#include "../PathDataClass.h"
#include "MultiStage.h"

class RefSliceShiftStageClass : public CommonStageClass
{
private:
  int oldRefSlice;

public:
  int maxShift;
  void Accept();
  void Reject();

  double Sample (int &slice1, int &slice2, Array <int,1> &activeParticles);
  RefSliceShiftStageClass (PathDataClass &pathData,IOSectionClass &outSection) :
    CommonStageClass (pathData,outSection)
  {}
};


class RefSliceShiftClass : public MultiStageClass
{
private:
  int maxShift;
  RefSliceShiftStageClass RefSliceShiftStage;
  Array<int,1> activeSpecies;
public:
  // Read the parameters from the input file
  void Read (IOSectionClass &in);
  // Actually attempts the move and accepts or rejects
  void MakeMove();
  RefSliceShiftClass (PathDataClass &pathData, IOSectionClass &outSection) :
    MultiStageClass(pathData, outSection), RefSliceShiftStage(pathData, outSection)
  {}
};

#endif
