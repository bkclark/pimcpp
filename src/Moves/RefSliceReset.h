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

#ifndef REF_SLICE_RESET_H
#define REF_SLICE_RESET_H

#include "MoveBase.h"

// Puts all particle positions at the reference slice particle position
class RefSliceResetClass : public ParticleMoveClass
{
public:
  void MakeMove();
  void Read(IOSectionClass &in);

  RefSliceResetClass (PathDataClass &myPathData, IOSectionClass outSection) :
    ParticleMoveClass(myPathData, outSection)
  {}

};


#endif
