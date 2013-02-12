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

#ifndef FORCES_H
#define FORCES_H

#include "ObservableBase.h"

/// A pair correlation function observable.
class ForcesClass : public ObservableClass
{
private:
  ObservableVecDouble2 ForcesVar;
  /// Stores number of counts in each bin
  Array<dVec,1> Forces, SumTmp;
  Array<double,2> ForcesArray;
  int SpeciesNum;
  Array<int,1> Ptcls;
  /// Stores the total number of counts
  int TimesCalled, Counts;
  int Freq;
  int DumpFreq;
public:
  /// My specialization of the virtual function.
  void Accumulate();
  /// My specialization of the virtual function.
  void Initialize();
  /// My specialization of the virtual function.
  void Print();
  void WriteBlock();
  void WriteInfo();
  void Read(IOSectionClass& IO);
  void SetSpecies (int speciesNum);
  ForcesClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), Counts(0),
    ForcesVar("F", IOSection, myPathData.Path.Communicator)
  {
    TimesCalled=0;
  }
  ForcesClass(PathDataClass &myPathData, IOSectionClass &ioSection,
	      int speciesNum) : 
    ObservableClass(myPathData, ioSection), Counts(0),
    ForcesVar("y", IOSection, myPathData.Path.Communicator)
  { 
    SetSpecies (speciesNum);
  }
};

#endif
