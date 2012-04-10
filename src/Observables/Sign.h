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

#ifndef SIGN_H
#define SIGN_H



#include "ObservableBase.h"




class SignClass : public ObservableClass
{

private:
  double Weight;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
  ObservableDouble Tot;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);
  SignClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection),
      Tot("Total",IOSection,PathData.Path.Communicator)
  {

    Weight = 0.0;
    NumSamples = 0;
    TimesCalled=0;
  }
};


#endif 
