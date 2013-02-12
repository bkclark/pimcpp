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

#ifndef PERMUTATIONCOUNT_H
#define PERMUTATIONCOUNT_H

#include "ObservableBase.h"

class PermutationCountClass : public ObservableClass
{
 private:
  ObservableVecDouble1 SectorCountVar;
  ObservableVecDouble1 CycleCountVar;
  ObservableVecInt2 PossPermsVar;
  vector<int> SectorCount;
  Array<double,1> CycleCount;
  int NumSamples;
  int Species;
public:
  void Accumulate();
  void Read(IOSectionClass& in);
  void WriteBlock();
  PermutationCountClass(PathDataClass &myPathData, IOSectionClass &ioSection) :
    ObservableClass(myPathData,ioSection),
    PossPermsVar("PossPerms", IOSection, myPathData.Path.Communicator),
    CycleCountVar("y", IOSection, myPathData.Path.Communicator),
    SectorCountVar("SectorCount",IOSection,myPathData.Path.Communicator)
  {
    NumSamples = 0;
    //SectorCount = 0.0;
    CycleCount = 0.0;
  }

};


#endif
