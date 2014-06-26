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

#ifndef SHIFT_MOVE_H
#define SHIFT_MOVE_H

#include "MoveBase.h"

///This is a "psuedo-move" that is inherited from MoveClass.
///"Pseudo-moves do not increment the monte-carlo time slice. It
///shifts the data in memory by numTimeSlicesToShift. We do this to
///shift data between processors, etc.
class ShiftMoveClass : public MoveClass
{
 public:
  /// Contains the number of time slices to shift at a time. We have this change itself randomly from other objects.
  int numTimeSlicesToShift;

  /// Function to actually make a shift move.
  void MakeMove();

  //Currently we don't read anything for the shift move class.
  void Read(IOSectionClass &in);
  ShiftMoveClass (PathDataClass &myPathData, IOSectionClass outSection) : 
    MoveClass(myPathData, outSection)
  {}

};


#endif
