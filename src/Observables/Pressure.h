/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
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
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef PRESSURE_H
#define PRESSURE_H

#include "ObservableBase.h"

class PressureClass : public ObservableClass
{
protected:
  ObservableDouble PressureVar, KineticVar, ShortRangeVar, LongRangeVar,
    NodeVar;
  ObservableVecDouble1 ShortRangeVecVar;
  ObservableVecDouble1 KineticVecVar;
  double Psum;
  double KineticSum, ShortRangeSum, LongRangeSum, NodeSum;
  int NumSamples;
  double PartitionPressure();
  double KineticPressure();
  double ShortRangePressure();
  double LongRangePressure();
  double NodePressure();
  Array<double,1> ShortRangePressureVec;
  Array<double,1> KineticPressureVec;
public:
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass &in);
  PressureClass(PathDataClass &pathData, IOSectionClass &ioSection)
    : ObservableClass (pathData, ioSection),
      PressureVar  ("Total",      IOSection, pathData.Path.Communicator),
      KineticVar   ("Kinetic",    IOSection, pathData.Path.Communicator),
      ShortRangeVar("ShortRange", IOSection, pathData.Path.Communicator),
      ShortRangeVecVar("ShortRangeComponents",IOSection,pathData.Path.Communicator),
      KineticVecVar("KineticComponents",IOSection,pathData.Path.Communicator),
      LongRangeVar ("LongRange",  IOSection, pathData.Path.Communicator),
      NodeVar      ("Node",       IOSection, pathData.Path.Communicator),
      NumSamples(0), Psum(0.0), KineticSum(0.0), ShortRangeSum(0.0),
      LongRangeSum(0.0), NodeSum(0.0)
  {

  }
};

#endif
