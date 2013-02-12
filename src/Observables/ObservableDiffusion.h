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

#ifndef MS_DIFFUSION_H
#define MS_DIFFUSION_H

#include "ObservableBase.h"

/// Measures the mean squared diffusion of the ptcl/COM
/// as a function of MC timestep.
class ObsDiffusionClass : public ObservableClass
{
  ObservableDouble HistVar;
  ObservableDouble TimeStepVar;
  Array<dVec, 2> R0;
	double TotalMSD;
  int TotalCounts;
  int TimesCalled;
	string COMSpecies;
	int totalSlices;
	int totalMol;
  int dumpFrequency;

public:
  /// My specialization of the virtual function.
  void Accumulate();
  /// My specialization of the virtual function.
  void Initialize();
  /// My specialization of the virtual function.
  void Print();
  void WriteBlock();
  void LocalWriteBlock();
  void WriteInfo();
  void Read(IOSectionClass& IO);
  ObsDiffusionClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    HistVar ("MeanSqDisp",  IOSection, myPathData.Path.Communicator),
    TimeStepVar ("time",  IOSection, myPathData.Path.Communicator)
  {
    TimesCalled=0;
    Initialize();
  }


};
#endif
