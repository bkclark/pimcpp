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
#include "../Communication/Communication.h"
#include "Sign.h"


// Fix to include final link between link M and 0
void SignClass::Accumulate()
{

  NumSamples++;
  double FullSign = 0.0;
  double currSign=PathData.Path.SignWeight;
  PathData.Path.Communicator.GatherProd(currSign,FullSign,0);
  Sign=Sign+FullSign;

}

void SignClass::ShiftData (int NumTimeSlices)
{
  // Do nothing
}

void SignClass::WriteBlock()
{
  Sign=Sign/(double)NumSamples;
  Tot.Write(Sign);
  if (PathData.Path.Communicator.MyProc()==0)
    IOSection.FlushFile();
  Sign = 0.0;
  NumSamples = 0;
}

void SignClass::Read(IOSectionClass &in)
{
  ObservableClass::Read(in);
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
}
