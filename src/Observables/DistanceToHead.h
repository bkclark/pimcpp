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

#ifndef DISTANCE_TO_HEAD_H
#define DISTANCE_TO_HEAD_H

#include "ObservableBase.h"

class HeadLocClass : public ObservableClass
{

private:
  //Must be initialized
  Array<int,1> HeadLoc;
  Array<int,1> TailLoc;
  ///This is the set of locations you should compare against to decide
  ///the location of the head and the tail
  Array<dVec,1> FixedLoc;
  ObservableVecDouble1 HeadLocVar,TailLocVar;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  void WriteInfo();
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& in);
  HeadLocClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      HeadLocVar("HeadLocation",IOSection,myPathData.Path.Communicator),
      TailLocVar("TailLocation",IOSection,myPathData.Path.Communicator)
  {
    NumSamples = 0;
    TimesCalled=0;
  }
};


#endif 
