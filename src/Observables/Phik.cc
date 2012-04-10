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

#include "Phik.h"



void PhiKClass::Accumulate()
{
  //Move the join to the end so we don't have to worry about permutations
  PathData.MoveJoin(PathData.NumTimeSlices()-1);
  
  NumSamples++;

  int slice1 = 0;
  int slice2 = PathData.Path.NumTimeSlices()-1;
  int numTimeSlices=PathData.Path.NumTimeSlices();
  PathClass &Path= PathData.Path;
  for (int tau=0;tau<numTimeSlices;tau++)
    for (int slice=slice1;slice<slice2;slice++)
      PhiK(tau) += Path(slice,0)[0]*Path((slice+tau) % PathData.Path.TotalNumSlices,0)[0];
}


void PhiKClass::WriteBlock()
{
  int nslices=PathData.Path.TotalNumSlices-1;
  double norm = 1.0/((double)NumSamples*(double)nslices);
  PhiK=PhiK*norm;
  PhiKVar.Write(PhiK);  
  PhiK=0;
  NumSamples = 0;
}

void PhiKClass::Read(IOSectionClass &in)
{  
  ObservableClass::Read(in);
  if (PathData.Path.Communicator.MyProc()==0){
    WriteInfo();
    IOSection.WriteVar("Type","Scalar");
  }
}

