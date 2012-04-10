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

#ifndef WEIGHT_H
#define WEIGHT_H



#include "ObservableBase.h"




class WeightClass : public ObservableClass
{

private:
  Array<double,1> Weight;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
  ObservableVecDouble1 Tot;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);
  WeightClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection),
      Tot("Total",IOSection,PathData.Path.Communicator)
  {
    Weight.resize(4);
    Weight = 0.0;
    NumSamples = 0;
    TimesCalled=0;
  }
};


#endif 
