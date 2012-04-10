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

#ifndef META_MOVES_H
#define META_MOVES_H


#include "MoveBase.h"



class PrintMoveClass : public MoveClass
{
  string MyString;
 public:
  void Read(IOSectionClass &IO);
  double AcceptanceRatio() {return 1.0;}
  void MakeMove() {cerr<<"This is printing  "<<MyString<<endl;}
  PrintMoveClass(PathDataClass &myPathData, IOSectionClass outSection) : 
    MoveClass(myPathData, outSection)
    {MyString="Hi";}
		    
};


////////////////////////////////////////////////////////////////////////////////

class JoinMoveClass : public MoveClass
{
  int JoinLocation;
  void MoveJoin()
  {
    //    PathData.MoveJoin();
  }
  void Read(IOSectionClass &input){};

  void MakeMove() {};
  JoinMoveClass (PathDataClass &myPathData, IOSectionClass outSection) : 
    MoveClass(myPathData, outSection)
  {/* Do nothing for now. */ }

};



////////////////////////////////////////////////////////////////////////////////



///This is a "psuedo-move" that is inherited from MoveClass.
///"Pseudo-moves do not increment the monte-carlo time slice. It
///shifts the data in memory by numTimeSlicesToShift. We do this to
///shift data between processors, etc.
class ShiftMoveClass : public MoveClass
{
 public:
  /// Contains the number of time slices to shift at a time. We have
///this change itself randomly from other objects.
  int numTimeSlicesToShift;
  /// Function to actually make a shift move. 
  void MakeMove();
  //Currently we don't read anything for the shift move class.
  void Read(IOSectionClass &theInput);
  ShiftMoveClass (PathDataClass &myPathData, IOSectionClass outSection) : 
    MoveClass(myPathData, outSection)
  { /* Do nothing for now. */ }
  
};




#endif
