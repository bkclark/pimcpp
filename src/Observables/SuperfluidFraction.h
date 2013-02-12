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

#ifndef SUPERFLUID_FRACTION__H
#define SUPERFLUID_FRACTION__H

#include "ObservableBase.h"
#include "WindingNumber.h"



class SuperfluidFractionClass : public WindingNumberClass
{
 protected:
  ObservableVecDouble1 SFVar;
 public:
  void Read(IOSectionClass& IO);
  void WriteBlock();
  SuperfluidFractionClass(PathDataClass &myPathData,IOSectionClass &ioSection) :
    WindingNumberClass(myPathData,ioSection),
    SFVar("SuperfluidFraction", IOSection, myPathData.Path.Communicator)
  {
    W2Sum = 0.0;
    SamplesInBlock = 0;
  }

};



#endif
