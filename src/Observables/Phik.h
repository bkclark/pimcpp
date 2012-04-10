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

#ifndef PHI_K_H
#define PHI_K_H

#include "ObservableBase.h"

class PhiKClass : public ObservableClass
{

private:
  Array<double,1> PhiK;
  ObservableVecDouble1 PhiKVar;
  
  

  int NumSamples;
  int TimesCalled;
public:
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& in);
  PhiKClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      PhiKVar("PhiK",IOSection,myPathData.Path.Communicator)

  {
    PhiK.resize(PathData.Path.NumTimeSlices());
    PhiK=0;
    TimesCalled=0;
    
  }
};


#endif 
