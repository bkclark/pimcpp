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

#include "MetaMoves.h"



void PrintMoveClass::Read(IOSectionClass &IO)
{
  string typeCheck;
  assert(IO.ReadVar("Type",typeCheck));
  assert(typeCheck=="PrintMove");
  assert(IO.ReadVar("toprint",MyString));
}

//////////////////////////////////////////////////////////////

void ShiftMoveClass::Read(IOSectionClass &theInput)
{
  string typeCheck;
  assert(theInput.ReadVar("Type",typeCheck));
  assert(typeCheck=="ShiftMove");
}


void ShiftMoveClass::MakeMove()
{
  int slice1, slice2;
  // The last processor will have the least number of slices
  // possible.  Use that for the maximum shift.
  PathData.Path.SliceRange(PathData.Path.Communicator.NumProcs()-1,
			   slice1, slice2);

//   /// Shift between 7/16 and 1/2 of the maximum slices
//   int maxSlices = (slice2-slice1)>>1;
//   int minSlices = (7*(slice2-slice1))>>4;
//   if ((maxSlices-minSlices)<5)
//     minSlices = max(0, maxSlices - 5);

//   int numTimeSlicesToShift = 
//     minSlices + PathData.Path.Random.CommonInt(maxSlices-minSlices+1);
  
  int maxSlices=slice2-slice1;
  int numTimeSlicesToShift = PathData.Path.Random.CommonInt(maxSlices);
  //  There is no point in shifting more than maxSlices/2
  if (numTimeSlicesToShift > (maxSlices>>1))
    numTimeSlicesToShift -= (maxSlices>>1);


  PathData.MoveJoin(0);
  PathData.ShiftData(numTimeSlicesToShift);
  PathData.Join=numTimeSlicesToShift;
}
