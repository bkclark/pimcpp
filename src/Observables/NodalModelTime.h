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

#ifndef NODALMODELTIME_H
#define NODALMODELTIME_H

#include "ObservableBase.h"
#include "../MatrixOps/MatrixOps.h"

class NodalModelTimeClass : public ObservableClass
{
private:
  ObservableVecDouble1 NodalModelTimeVar;
  Array<double,1> counts;
  int NumSamples;
  string Species;
public:
  /// This grid defines the bins.  Bin 0 is bounded by 0 on the 
  /// bottom and grid(0) on the top.
  void Accumulate();
  void Read(IOSectionClass& in);
  void WriteBlock();
  void WriteInfo();
  NodalModelTimeClass(PathDataClass &myPathData, IOSectionClass &ioSection) :
    ObservableClass(myPathData,ioSection),
    NodalModelTimeVar("Time", IOSection, myPathData.Path.Communicator)
  {
    NumSamples = 0;
  }

};


#endif
