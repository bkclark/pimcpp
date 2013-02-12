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

#ifndef WINDING_NUMBER__H
#define WINDING_NUMBER__H

#include "ObservableBase.h"



class WindingNumberClass : public ObservableClass
{
 protected:
  ObservableVecDouble1 WNVar;
  ObservableVecDouble1 WNVarLowVariance;
  /// This stores a list of integers corresponding to the species that
  /// are included in the winding number calculation.
  Array<int,1> SpeciesList;
  /// This stores the block sum of the winding numbers
  dVec W2Sum;
  /// Vector which stores all of the local winding numbers since the
  /// last writeblock.
  std::vector<dVec> WNVec;
  Array<double,1> WN2Array;
  Array<double,1> WN2ArrayLowVariance;
  int SamplesInBlock;
  void CalcWN2();
 public:
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& IO);
  WindingNumberClass(PathDataClass &myPathData,IOSectionClass &ioSection) :
    ObservableClass(myPathData,ioSection),
    WNVar("W2", IOSection, myPathData.Path.Communicator),
    WNVarLowVariance("W2LowVariance", IOSection, myPathData.Path.Communicator)
  {
    W2Sum = 0.0;
    SamplesInBlock = 0;
  }

};



#endif
