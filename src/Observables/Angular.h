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

#ifndef ANGULAR_H
#define ANGULAR_H

#include "ObservableBase.h"

class AngularClass : public ObservableClass
{
 private:
  int Freq,DumpFreq;
  int TotalCounts;
  int TimesCalled;
  ObservableVecDouble2 CorVar;
  
public:
  ///l x t
  Array<double,2> Correlation;
  int Species;
  void Accumulate();
  void Read(IOSectionClass& in);
  void WriteBlock();
  AngularClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    CorVar("y", IOSection, myPathData.Path.Communicator)
  {
    TimesCalled=0;
  }

};


#endif
