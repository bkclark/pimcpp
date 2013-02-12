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

#ifndef SUPERFLUIDROP_H
#define SUPERFLUIDROP_H

#include "ObservableBase.h" // what are we inheriting?

class SuperfluiDrop : public ObservableClass
{
 private:
  int Freq,DumpFreq;
  int TotalCounts;
  int TimesCalled;
  //double superfluidity;
  ObservableDouble areaSquared,momInertia;
public:
  double area,anorm;
  double momi,mominorm;
  int Species;
  void Accumulate();
  void Read(IOSectionClass& in);
  void WriteBlock();
  SuperfluiDrop(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    areaSquared("area", IOSection, myPathData.Path.Communicator),
    momInertia("mominertia", IOSection, myPathData.Path.Communicator)
  {
    TimesCalled=0;
  }
};


#endif
