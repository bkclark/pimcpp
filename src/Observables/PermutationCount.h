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

#ifndef PERMUTATIONCOUNT_H
#define PERMUTATIONCOUNT_H

#include "ObservableBase.h"

class PermutationCountClass : public ObservableClass
{
 private:
  Array<double,1> PermutationNumber;
  ObservableVecDouble1 PermutationNumberVar;
  int Freq,DumpFreq;
  int TotalCounts;
  int TimesCalled;
  ObservableVecDouble1 CycleCountVar;
  Array<bool,1> CountedAlready;
  Array<int,1> TotalPerm;

  ObservableVecDouble1 SectorCountVar;
  Array<double,1> SectorCount;
  int NumSamples;
public:
  Array<double,1> CycleCount;
  int Species;
  void Accumulate();
  void Read(IOSectionClass& in);
  void WriteBlock();
  PermutationCountClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    PermutationNumberVar("Partition Function",IOSection,myPathData.Path.Communicator),
    CycleCountVar("y", IOSection, myPathData.Path.Communicator),
    SectorCountVar  ("SectorCount",IOSection,myPathData.Path.Communicator)
  {
    TimesCalled=0;
    SectorCount = 0.0;
  }

};


#endif
