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

#ifndef COUPLING_H
#define COUPLING_H

#include "ObservableBase.h"

/// A pair correlation function observable.
class CouplingClass : public ObservableClass
{
private:
  ObservableVecDouble1 couplingVar;
  /// Stores number of counts in each bin
  Array<int,1> Coupling;
  /// Stores the total number of counts
  int TotalCounts;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  /// My specialization of the virtual function.
  void Accumulate();
  /// My specialization of the virtual function.
  void Initialize();
  /// My specialization of the virtual function.
  void WriteBlock();
  void WriteInfo();
  void Read(IOSectionClass& IO);
  CouplingClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    couplingVar("y", IOSection, myPathData.Path.Communicator)
  {
    Initialize();
  }


};


#endif
