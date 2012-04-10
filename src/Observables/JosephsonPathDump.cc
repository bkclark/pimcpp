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

#include "JosephsonPathDump.h"



void JosephsonPathDumpClass::Accumulate()
{

}


void JosephsonPathDumpClass::WriteBlock()
{
  int slice1 = 0;
  int slice2 = PathData.Path.NumTimeSlices();
  int numTimeSlices=PathData.Path.NumTimeSlices();
  PathClass &Path= PathData.Path;
  for (int slice=slice1;slice<slice2;slice++)
    Phase(slice)=Path(slice,0)[0];
  PhaseVar.Write(Phase);  
}

void JosephsonPathDumpClass::Read(IOSectionClass &in)
{  
  ObservableClass::Read(in);
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
}

